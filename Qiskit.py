from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.qasm2 import dumps
from typing import Tuple, Dict
from qiskit.quantum_info import DensityMatrix, partial_trace
import numpy as np
app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class Gate(BaseModel):
    type: str
    params: list[int] = []
    angle: float | None = None

class CircuitRequest(BaseModel):
    numQubits: int
    gates: list[Gate]
    initialStates: list[int] | None = None
    targetQubit : int | None = None
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

def reconstruct_single_qubit_rho(base_qc: QuantumCircuit,
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


def reconstruct_all_qubits(base_qc: QuantumCircuit,
                           shots: int = 2048,
                           backend=None):
    """Reconstruct reduced density matrices for all qubits."""
    n = base_qc.num_qubits
    backend = _ensure_backend(backend)

    circuits, meta = [], []
    for target in range(n):
        for basis in ("Z", "X", "Y"):
            circuits.append(_make_meas_circuit_variant(base_qc, target, basis))
            meta.append((target, basis))

    tcirc = transpile(circuits, backend=backend)
    result = backend.run(tcirc, shots=shots).result()

    raw_expectations = {q: {} for q in range(n)}
    for idx, (target, basis) in enumerate(meta):
        counts = result.get_counts(idx)
        p0, p1 = _counts_bit_for_qubit(counts, n, target)
        raw_expectations[target][basis] = _expectation_from_probs(p0, p1)

    bloch_list, reduced_rhos = [], []
    for q in range(n):
        x, y, z = raw_expectations[q]["X"], raw_expectations[q]["Y"], raw_expectations[q]["Z"]
        rho = _reconstruct_rho_from_xyz(x, y, z)
        bloch_list.append({"x": float(x), "y": float(y), "z": float(z)})
        reduced_rhos.append(rho)

    return bloch_list, reduced_rhos, raw_expectations

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
        x,y,z, rho = reconstruct_single_qubit_rho(qc, target = request.targetQubit,shots=8192)
        clean_rho, bloch = clean_rho_and_bloch(rho)
        return{
        "blochs": bloch,
        "rho" : serialize_rho(clean_rho),
        }
