# encoding: utf-8
# module synpy_core.synpy_rust

from typing import Any

def synthesize_pauli_exponential(*args: Any, **kwargs: Any) -> Any: ...

class PyCommand(object):
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    def __new__(*args: Any, **kwargs: Any) -> Any: ...
    CX: Any
    CZ: Any
    H: Any
    Rx: Any
    Ry: Any
    Rz: Any
    S: Any
    SDgr: Any
    V: Any
    VDgr: Any
    X: Any
    Y: Any
    Z: Any

class PyPauliString(object):
    def as_tuple(self, *args: Any, **kwargs: Any) -> Any: ...
    def get_qubits(self, *args: Any, **kwargs: Any) -> Any: ...
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    def __new__(*args: Any, **kwargs: Any) -> Any: ...

__all__ = [
    "PyPauliString",
    "PyCommand",
    "synthesize_pauli_exponential",
]

__loader__: Any = ...
__spec__: Any = ...
