# `synpy`: Python bindings for syn

Build and distribute using maturin:

```bash
pip install maturin
maturin build --release
```

The resulting wheel file can be installed with pip:

```bash
pip install target/wheels/synpy-*.whl
```

```python
from synpy import CliffordTableau

tableau1 = CliffordTableau(3)

tableau2 = CliffordTableau.with_pauli_columns(
    3, ["XIIZII", "IXIIZI", "IIXIIZ"], [False] * 6
)

assert tableau1 == tableau2
```