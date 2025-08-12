from synpy.synpy_rust import PyCommand


def tuple_to_pycommand(command_tuple: tuple) -> PyCommand:
    gate = command_tuple[0]
    if gate == "H":
        return PyCommand.H(command_tuple[1])
    elif gate == "X":
        return PyCommand.X(command_tuple[1])
    elif gate == "Y":
        return PyCommand.Y(command_tuple[1])
    elif gate == "Z":
        return PyCommand.Z(command_tuple[1])
    elif gate == "S":
        return PyCommand.S(command_tuple[1])
    elif gate == "SDgr":
        return PyCommand.SDgr(command_tuple[1])
    elif gate == "V":
        return PyCommand.V(command_tuple[1])
    elif gate == "VDgr":
        return PyCommand.VDgr(command_tuple[1])
    elif gate == "Rx":
        return PyCommand.Rx(command_tuple[1], command_tuple[2])
    elif gate == "Ry":
        return PyCommand.Ry(command_tuple[1], command_tuple[2])
    elif gate == "Rz":
        return PyCommand.Rz(command_tuple[1], command_tuple[2])
    elif gate == "CX":
        return PyCommand.CX(command_tuple[1], command_tuple[2])
    elif gate == "CZ":
        return PyCommand.CZ(command_tuple[1], command_tuple[2])
    else:
        raise ValueError("Unknown gate type: " + gate)


def pycommand_to_tuple(cmd: PyCommand) -> tuple:
    if isinstance(cmd, PyCommand.H):
        return "H", cmd[0]
    elif isinstance(cmd, PyCommand.X):
        return "X", cmd[0]
    elif isinstance(cmd, PyCommand.Y):
        return "Y", cmd[0]
    elif isinstance(cmd, PyCommand.Z):
        return "Z", cmd[0]
    elif isinstance(cmd, PyCommand.S):
        return "S", cmd[0]
    elif isinstance(cmd, PyCommand.SDgr):
        return "SDgr", cmd[0]
    elif isinstance(cmd, PyCommand.V):
        return "V", cmd[0]
    elif isinstance(cmd, PyCommand.VDgr):
        return "VDgr", cmd[0]
    # Rotation gates retrieve the angle by assignment.
    elif isinstance(cmd, PyCommand.Rx):
        return "Rx", cmd[0], cmd[1]
    elif isinstance(cmd, PyCommand.Ry):
        return "Ry", cmd[0], cmd[1]
    elif isinstance(cmd, PyCommand.Rz):
        return "Rz", cmd[0], cmd[1]
    # Two-qubit gates:
    elif isinstance(cmd, PyCommand.CX):
        return "CX", cmd[0], cmd[1]
    elif isinstance(cmd, PyCommand.CZ):
        return "CZ", cmd[0], cmd[1]
    else:
        raise ValueError("Unhandled PyCommand variant: " + str(type(cmd)))


def pycommand_to_qasm(n_qubits: int, commands: list[PyCommand]) -> str:
    out = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
        f"qreg q[{n_qubits}];",
    ]

    for command in commands:
        op, *args = pycommand_to_tuple(command)
        op = op.lower()
        args = [f"q[{arg}]" for arg in args]
        s = f'{op} {", ".join(args)};'
        out.append(s)

    return "\n".join(out)
