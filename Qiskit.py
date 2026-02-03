from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
from qiskit import QuantumCircuit, transpile
from qiskit_experiments.framework import ParallelExperiment
from qiskit_experiments.library import StateTomography
from qiskit_aer import AerSimulator
from qiskit.qasm2 import dumps
from typing import Tuple, Dict, List, Optional
from qiskit.quantum_info import DensityMatrix, partial_trace
from qiskit.circuit.library import UnitaryGate
import numpy as np

# Try importing qiskit-experiments (recommended). If not available, we'll fallback.
try:
    from qiskit_experiments.library import StateTomography
    HAS_QISKIT_EXPERIMENTS = True
except Exception:
    HAS_QISKIT_EXPERIMENTS = False


app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class ComplexNumber(BaseModel):
    re: float
    im: float

class Gate(BaseModel):
    type: str
    name: Optional[str] = None           # "m", "cx", etc.
    customType: Optional[str] = None     # "CUSTOM_MATRIX", "CUSTOM_CIRCUIT", "CUSTOM_CONTROL"
    params: list[int] = []
    angle: float | None = None
    matrix: Optional[List[List[ComplexNumber]]] = None
    subGates: Optional[List["Gate"]] = None   # recursive definition for CUSTOM_CIRCUIT

class CircuitRequest(BaseModel):
    numQubits: int
    gates: list[Gate]
    initialStates: list[int] | None = None
    targetQubit : int | None = None

def decode_complex_matrix(mat: List[List["ComplexNumber"]]) -> np.ndarray:
    rows, cols = len(mat), len(mat[0])
    out = np.zeros((rows, cols), dtype=complex)
    for i in range(rows):
        for j in range(cols):
            out[i, j] = complex(mat[i][j].re, mat[i][j].im)
    return out


def serialize_rho(rho: np.ndarray):
    """Convert a density matrix with complex numbers to JSON-safe format."""
    return [[{"re": float(np.real(val)), "im": float(np.imag(val))} for val in row] for row in rho]

def _ensure_backend(backend=None):
    """Return AerSimulator if backend is not provided."""
    return backend or AerSimulator()


def _counts_bit_for_qubit(counts: Dict[str, int], n_qubits: int, target: int) -> Tuple[float, float]:
    """Compute probabilities P(0), P(1) for target qubit given counts dict."""
    total, c0 = 0, 0
    for bitstr, cnt in counts.items():
        total += cnt
        rev = bitstr[::-1]  # reverse bitstring so index matches qubit
        if rev[target] == '0':
            c0 += cnt
    if total == 0:
        return 0.0, 0.0
    p0 = c0 / total
    p1 = 1.0 - p0
    return p0, p1


def _expectation_from_probs(p0: float, p1: float) -> float:
    return float(p0 - p1)


def _reconstruct_rho_from_xyz(x: float, y: float, z: float) -> np.ndarray:
    """Reconstruct density matrix rho = 1/2 (I + xX + yY + zZ)."""
    rho = 0.5 * np.array([[1 + z, x - 1j * y],
                          [x + 1j * y, 1 - z]], dtype=complex)
    return rho


def _make_meas_circuit_variant(base_qc: QuantumCircuit, target: int, basis: str) -> QuantumCircuit:
    """Return a copy of base_qc with measurement in chosen basis (Z, X, Y)."""
    qc = base_qc.copy()
    if basis == "X":
        qc.h(target)
    elif basis == "Y":
        qc.sdg(target)
        qc.h(target)
    elif basis == "Z":
        pass
    else:
        raise ValueError("basis must be 'X','Y' or 'Z'")
    qc.measure_all()
    return qc


# -------------------- Tomography functions --------------------

def reconstruct_single_qubit_rho_manual(base_qc: QuantumCircuit,
                                 target: int,
                                 shots: int = 2048,
                                 backend=None):
    """Reconstruct reduced density matrix for a single qubit."""
    n = base_qc.num_qubits
    if not (0 <= target < n):
        raise ValueError("Target qubit out of range.")

    backend = _ensure_backend(backend)
    circuits, bases = [], []

    for basis in ("Z", "X", "Y"):
        circuits.append(_make_meas_circuit_variant(base_qc, target, basis))
        bases.append(basis)

    tcirc = transpile(circuits, backend=backend)
    result = backend.run(tcirc, shots=shots).result()

    exps = {}
    for idx, basis in enumerate(bases):
        counts = result.get_counts(idx)
        p0, p1 = _counts_bit_for_qubit(counts, n, target)
        exps[basis] = _expectation_from_probs(p0, p1)

    x, y, z = exps["X"], exps["Y"], exps["Z"]
    rho = _reconstruct_rho_from_xyz(x, y, z)
    return x, y, z, rho
def reconstruct_single_qubit_rho_experiments(base_qc: QuantumCircuit,
                                             target: int,
                                             shots: int = 8192,
                                             backend=None) -> Tuple[float, float, float, np.ndarray]:
    """
    Use qiskit-experiments StateTomography when available.
    Returns (x,y,z,rho_numpy_array).
    """
    if not HAS_QISKIT_EXPERIMENTS:
        raise RuntimeError("qiskit-experiments not installed")

    if target is None:
        raise ValueError("target must be specified for tomography")

    n = base_qc.num_qubits
    if not (0 <= target < n):
        raise ValueError("Target qubit out of range.")

    backend = _ensure_backend(backend)

    # Build the experiment for tomography on a single qubit (subsystem)
    tomo = StateTomography(base_qc, [target])

    # Run experiment — ExperimentData object is returned
    exp_data = tomo.run(backend, shots=shots).block_for_results()

    # Try to extract a density matrix from the experiment analysis results.
    # Different versions of qiskit-experiments expose results differently; try a few reasonable access patterns.
    rho_np = None
    try:
        # pattern 1: query specifically for "state" analysis (works in some releases)
        ar = exp_data.analysis_results("state")
        # ar may be a single AnalysisResult or a list; normalize
        if isinstance(ar, (list, tuple)):
            ar = ar[0] if len(ar) > 0 else None
        if ar is not None:
            # Many AnalysisResult objects store the reconstructed state in `.value`
            if hasattr(ar, "value") and isinstance(ar.value, (DensityMatrix, np.ndarray)):
                rho_obj = ar.value
                rho_np = rho_obj.data if isinstance(rho_obj, DensityMatrix) else np.array(rho_obj)
    except Exception:
        rho_np = None

    # fallback attempt: iterate over all analysis results and try to find a DensityMatrix-like value
    if rho_np is None:
        try:
            for ar in exp_data.analysis_results():
                # each ar may have .value or .metadata etc.
                if hasattr(ar, "value") and ar.value is not None:
                    val = ar.value
                    if isinstance(val, DensityMatrix):
                        rho_np = val.data
                        break
                    # sometimes value might be a dict with a 'state' field or similar:
                    if isinstance(val, dict):
                        # try common keys
                        for key in ("density_matrix", "state", "rho", "reconstructed_state"):
                            if key in val:
                                v = val[key]
                                if isinstance(v, DensityMatrix):
                                    rho_np = v.data
                                    break
                                elif isinstance(v, (np.ndarray, list)):
                                    rho_np = np.array(v)
                                    break
                        if rho_np is not None:
                            break
        except Exception:
            rho_np = None

    # last fallback: if we still didn't get a matrix, raise up so caller can fallback to manual
    if rho_np is None:
        raise RuntimeError("Could not extract reconstructed density matrix from qiskit-experiments result")

    # compute Bloch coordinates from density matrix
    rho_mat = np.array(rho_np, dtype=complex)
    x = 2 * np.real(rho_mat[0, 1])
    y = -2 * np.imag(rho_mat[0, 1])
    z = np.real(rho_mat[0, 0] - rho_mat[1, 1])
    return float(x), float(y), float(z), rho_mat

def apply_gate(qc: QuantumCircuit, gate: Gate):
    g, p, a = gate.type, gate.params, gate.angle

    if g == "X": qc.x(p[0])
    elif g == "Y": qc.y(p[0])
    elif g == "Z": qc.z(p[0])
    elif g == "H": qc.h(p[0])
    elif g == "S": qc.s(p[0])
    elif g == "Sdg": qc.sdg(p[0])
    elif g == "T": qc.t(p[0])
    elif g == "Tdg": qc.tdg(p[0])
    elif g == "Rx": qc.rx(a, p[0])
    elif g == "Ry": qc.ry(a, p[0])
    elif g == "Rz": qc.rz(a, p[0])
    elif g == "Phase": qc.p(a, p[0])
    elif g == "CNOT": qc.cx(p[0], p[1])
    elif g == "CZ": qc.cz(p[0], p[1])
    elif g == "SWAP": qc.swap(p[0], p[1])
    elif g == "CCNOT": qc.ccx(p[0], p[1], p[2])

    else:
        raise ValueError(f"Unsupported gate: {gate}")


def build_circuit(data: CircuitRequest):
    n = data.numQubits
    if n<6:
        qc = QuantumCircuit(n, n)
    else:
        qc = QuantumCircuit(n)
    if len(data.initialStates) < n:
        # Pad with 0s
        raise ValueError("less states")
    elif len(data.initialStates) > n:
        raise ValueError("Too many initialStates for given numQubits")
    if data.initialStates:
        for idx, state in enumerate(data.initialStates):
            if int(state) == 1:
                qc.x(idx)
    for gate in data.gates:
        g, p, a = gate.type, gate.params, gate.angle
        if g == "X": qc.x(p[0])
        elif g == "Y": qc.y(p[0])
        elif g == "Z": qc.z(p[0])
        elif g == "H": qc.h(p[0])
        elif g == "S": qc.s(p[0])
        elif g == "Sdg": qc.sdg(p[0])
        elif g == "T": qc.t(p[0])
        elif g == "Tdg": qc.tdg(p[0])
        elif g == "Rx": qc.rx(a, p[0])
        elif g == "Ry": qc.ry(a, p[0])
        elif g == "Rz": qc.rz(a, p[0])
        elif g == "Phase": qc.p(a, p[0])
        elif g == "CNOT": qc.cx(p[0], p[1])
        elif g == "CZ": qc.cz(p[0], p[1])
        elif g == "SWAP": qc.swap(p[0], p[1])
        elif g == "CCNOT": qc.ccx(p[0], p[1], p[2])
        elif g == "CUSTOM" and gate.customType == "CUSTOM_MATRIX":
            matrix_np = decode_complex_matrix(gate.matrix)
            unitary = UnitaryGate(matrix_np, label=gate.name)
            qc.append(unitary, p)
        elif g == "CUSTOM" and gate.customType == "CUSTOM_CIRCUIT":
            # Find how many qubits the subcircuit needs
            max_qubit = -1
            for sg in gate.subGates or []:
                if sg.params:
                    max_qubit = max(max_qubit, max(sg.params))
            num_sub_qubits = max_qubit + 1 if max_qubit >= 0 else len(p)

            sub_qc = QuantumCircuit(num_sub_qubits)

            for sg in gate.subGates or []:
                apply_gate(sub_qc, sg)  # ✅ now sub_qc has correct size

            qc.compose(sub_qc, p, inplace = True)
        elif g == "CUSTOM" and gate.customType == "CUSTOM_CONTROL":
            if not gate.subGates or len(gate.subGates) != 1:
                raise ValueError("CUSTOM_CONTROL must wrap exactly one subGate")
            
            sg = gate.subGates[0]

            # Frontend supplies exact qubit indices: last = target, rest = controls
            ctrl_qubits = p[:-1]
            target_qubit = p[-1]

            # Validate qubits are within circuit range
            all_qubits = ctrl_qubits + [target_qubit]
            if any(q >= data.numQubits or q < 0 for q in all_qubits):
                raise ValueError(f"Control or target qubit out of range: {all_qubits}")

            # Apply native multi-controlled gates
            if sg.type == "X":
                qc.mcx(ctrl_qubits, target_qubit)
            elif sg.type == "Y":
                qc.mcy(ctrl_qubits, target_qubit)
            elif sg.type == "Z":
                qc.mcz(ctrl_qubits, target_qubit)
            elif sg.type == "H":
                qc.mch(ctrl_qubits, target_qubit)
            elif sg.type == "Rx":
                qc.mcrx(sg.angle, ctrl_qubits, target_qubit)
            elif sg.type == "Ry":
                qc.mcry(sg.angle, ctrl_qubits, target_qubit)
            elif sg.type == "Rz":
                qc.mcrz(sg.angle, ctrl_qubits, target_qubit)
            elif sg.type == "Phase":
                qc.mcp(sg.angle, ctrl_qubits, target_qubit)
            elif sg.type == "S":
                # S gate controlled: apply as phase pi/2
                qc.mcp(np.pi/2, ctrl_qubits, target_qubit)
            elif sg.type == "T":
                # T gate controlled: apply as phase pi/4
                qc.mcp(np.pi/4, ctrl_qubits, target_qubit)
        else:
            raise ValueError(f"Unsupported controlled single-qubit gate: {sg.type}")


    # ✅ Always measure all qubits at the end
    if(n<6):
        qc.measure(range(n), range(n))

    return qc

def clean_rho_and_bloch(rho: np.ndarray, tol: float = 1e-12):
    """
    Take a possibly noisy single-qubit density matrix and return:
      - cleaned (physical) rho
      - Bloch vector coordinates (x,y,z)
    """
    # --- Step 1: enforce Hermiticity
    rho = 0.5 * (rho + rho.conj().T)

    # --- Step 2: eigen-decompose
    vals, vecs = np.linalg.eigh(rho)

    # --- Step 3: clip small negatives, renormalize
    vals = np.where(vals < tol, 0.0, vals)
    rho_clean = vecs @ np.diag(vals) @ vecs.conj().T
    rho_clean /= np.trace(rho_clean)

    # --- Step 4: compute Bloch coordinates
    # rho = 1/2 (I + xX + yY + zZ)
    x = 2.0 * np.real(rho_clean[0, 1])
    y = -2.0 * np.imag(rho_clean[0, 1])
    z = np.real(rho_clean[0, 0] - rho_clean[1, 1])

    return rho_clean, {"x": float(x), "y": float(y), "z": float(z)}

@app.get("/")
def home():
    return {"message": "Qiskit backend is running!"}

@app.post("/run")
def run_circuit(request: CircuitRequest):
    if(request.numQubits<6):
        qc = build_circuit(request)
        backend = AerSimulator()
        qc = transpile(qc, backend)
        job = backend.run(qc, shots=1024)
        result = job.result()
        counts = result.get_counts()

        return {
            "counts": counts,
            "qasm": dumps(qc),
        }
    else:
        qc = build_circuit(request)
        backend = AerSimulator()
        qc = transpile(qc, backend)
            # Preferred path: qiskit-experiments
        try:
            if HAS_QISKIT_EXPERIMENTS:
                x, y, z, rho = reconstruct_single_qubit_rho_experiments(qc, target=request.targetQubit,
                                                                        shots=8192, backend=backend)
            else:
                # fallback to manual tomography if qiskit-experiments not present
                x, y, z, rho = reconstruct_single_qubit_rho_manual(qc, target=request.targetQubit,
                                                                shots=8192, backend=backend)
        except Exception as e:
            # As a robust fallback: try manual approach if experiments path failed
            try:
                x, y, z, rho = reconstruct_single_qubit_rho_manual(qc, target=request.targetQubit,
                                                                shots=8192, backend=backend)
            except Exception as e2:
                raise HTTPException(status_code=500, detail=f"Tomography failed: {str(e)} | fallback failed: {str(e2)}")
        clean_rho, bloch = clean_rho_and_bloch(rho)
        return{
        "blochs": bloch,
        "rho" : serialize_rho(clean_rho),
        }
