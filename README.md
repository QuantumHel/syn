# SynIR

## Contents
- [1) What is SynIR?](#what-is-synir?)
- [2) Getting Started](#getting-started)
    - [Installing SynPy](#installing-synpy)
    - [Installing SynIR](#installing-synir)
    - [Setting Up Jupyter](#setting-up-jupyter)
- [3) Data Structures](#data-structures)
    - [Pauli Letter](#pauli-letter)
    - [Pauli String](#pauli-string)
    - [Pauli Polynomial](#pauli-polynomial)
    - [Pauli Exponential](#pauli-exponential)
    - [Clifford Tableau](#clifford-tableau)
- [4) Algorithms](#algorithms)
    - [WIP]()

---

## What is SynIR?
WIP

## Getting Started

### Installing SynPy

1. Install `pyenv` and `pyenv-virtualenv`
```bash
brew install pyenv pyenv-virtualenv
```

2. Run the `python_setup` script. This will create a virtualenv environment
```bash
make python_setup
```

3. Run the `python_project_setup` script. This will install required dependencies from `requirements.txt`
```bash
make python_project_setup
```

4. Add the following line to your `.zshrc` to enable `pyenv-virtualenv` to automatically activate the environment when in the project directory
```zsh
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
```

### Installing SynIR

1. Install [rust (recommended)](https://www.rust-lang.org/tools/install)
2. Run set-up scripts

```bash
make rust_project_setup
```

3. Build the project
```bash
make build
```

4. Run tests
```bash
make test
```

5. Format check
```bash
make format-check
```
6. Clean
```bash
make clean
```

### Setting up Jupyter

1. Install Jupyter
```bash
pip install jupyter
```

2. Install the rust kernel for Jupyter provided by [excxr](https://github.com/evcxr/evcxr)
```bash
cargo install evcxr_jupyter
```

3. Register the kernel with Jupyter
```bash
evcxr_jupyter --install
```

4. Start a Jupyter Notebook. When creating a new notebook, choose rust as the kernel
```bash
jupyter notebook
```

## Data Structures

### Pauli Letter
The `PauliLetter` data structure is the letter by which we refer to any of the four Pauli matrices: $X$, $Y$, $Z$, and $I$.

### Pauli String

A PauliString is a sequence of Pauli Letters representing the tensor product of the respective Pauli matrices. The `PaulString` data structure is defined as follows:
```rust
pub struct PauliString {
    pub(super) x: RwLock<BitVec>,
    pub(super) z: RwLock<BitVec>,
}
```

Here, the Pauli String `"XYZ"` is represented by the two thread-safe bitvectors, `x = [1, 1, 0]` and `z = [0, 1, 1]`. This encoding follows the binary symplectic representation. In this scheme, each position in the string is determined by a pair of bits `(xᵢ, zᵢ)`:
- `(0, 0)` → Identity
- `(1, 0)` → Pauli-X
- `(0, 1)` → Pauli-Z
- `(1, 1)` → Pauli-Y

This representation allows for fast bit-level operations and efficient Clifford updates.

### Pauli Polynomial

The Pauli Polynomial is the subsequence of Pauli Strings in the Pauli Exponential that we assume to be mutually commuting. An example of this is the Phase Polynomial, which is a Pauli Polynomial consisting only of $I$ and $Z$ Pauli Letters.

The naming scheme is analogous to the Phase Polynomial. Similarly, a polynomial contains a sum of terms and since the sum of matrices is commutative, the terms in the Pauli Polynomial should also commute.

The `PauliPolynomial` data structure is defined as follows:
```rust
pub struct PauliPolynomial {
    chains: Vec<PauliString>,
    angles: RwLock<Vec<Angle>>,
    size: usize,
}
```

### Pauli Exponential

The Pauli Exponential is the universal representation that we use in the library. It contains a Clifford Tableau and a sequence of Pauli Polynomials. An exponential is a product and since the product of matrices is not commutative, the sequences in the Pauli Exponential do not commute.

The `PauliExponential` data structure is defined as follows:
```rust
pub struct PauliExponential {
    pauli_polynomials: VecDeque<PauliPolynomial>,
    clifford_tableau: CliffordTableau,
}
```

### Clifford Tableau

The `CliffordTableau` data structure is defined as follows
```rust
pub struct CliffordTableau {
    pauli_columns: Vec<PauliString>,
    signs: BitVec,
}
```

## Algorithms
WIP