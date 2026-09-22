"""Conversions between :class:`PyCommand` and other gate representations.

A ``PyCommand`` is what synthesis hands back. The tuple form used here is
``(name, *qubits)`` for Clifford gates and ``(name, target, angle)`` for rotations, with
``name`` one of ``H X Y Z S SDgr V VDgr Rx Ry Rz CX CZ``.
"""

from synpy.synpy_rust import PyCommand


def tuple_to_pycommand(command_tuple: tuple) -> PyCommand:
    """Build a :class:`PyCommand` from a ``(name, *args)`` tuple.

    Parameters
    ----------
    command_tuple
        ``(name, target)`` for a single-qubit Clifford, ``(name, target, angle)`` for a
        rotation, or ``(name, control, target)`` for a two-qubit gate.

    Returns
    -------
    PyCommand

    Raises
    ------
    ValueError
        If ``command_tuple[0]`` is not a recognized gate name.
    """
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
    """Flatten a :class:`PyCommand` into a ``(name, *args)`` tuple.

    The inverse of :func:`tuple_to_pycommand`. Rotation angles come last, so an ``Rz`` becomes
    ``("Rz", target, angle)``.

    Raises
    ------
    ValueError
        If ``cmd`` is a variant this function does not handle.
    """
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
    """Render a command list as an OpenQASM 2.0 program over one register ``q``.

    Gate names are lowercased and every argument is emitted as a qubit reference, so this is
    only correct for the Clifford commands -- a rotation would render its angle as ``q[...]``.
    Use it for tableau synthesis output, which contains no rotations.

    Parameters
    ----------
    n_qubits
        Size of the declared quantum register.
    commands
        Gates in circuit order.

    Returns
    -------
    str
        A QASM program, parseable with ``QuantumCircuit.from_qasm_str``.
    """
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
