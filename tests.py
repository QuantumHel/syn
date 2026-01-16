import qiskit
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Operator

x = QuantumCircuit(1)
x.sdg(0)
x.x(0)
x.s(0)

op = Operator.from_circuit(x)
print(op)

x = QuantumCircuit(1)
x.y(0)

op = Operator.from_circuit(x)
print(op)
