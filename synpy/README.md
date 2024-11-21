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

tableau2 = CliffordTableau.from_parts(
    ["XIIZII", "IXIIZI", "IIXIIZ"], [False] * 6, 3
)

assert tableau1 == tableau2
assert tableau1.synthesize() == []

tableau3 = CliffordTableau.from_parts(
    ["XIIZZI", "IXIIZI", "IIXIIZ"], [False] * 6, 3
)

assert tableau3.synthesize() == ['CX 0 1']
```