from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.qasm2 import dumps

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

def build_circuit(data: CircuitRequest):
    n = data.numQubits
    qc = QuantumCircuit(n, n)
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
    qc.measure(range(n), range(n))

    return qc
@app.get("/")
def home():
    return {"message": "Qiskit backend is running!"}

@app.post("/run")
def run_circuit(request: CircuitRequest):
    qc = build_circuit(request)

    backend = AerSimulator()
    job = backend.run(qc, shots=1024)
    result = job.result()
    counts = result.get_counts()

    return {
        "counts": counts,
        "qasm": dumps(qc),
    }
