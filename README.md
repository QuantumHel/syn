# synIR



## Installation

### Python

We currently use python for pre-commit hooks. We herefore have a pyenv setup, after having pyenv and pyenv-virtualenv installed, simply run:

```
make python_setup
```

and

```
make python_project_setup
```

### Rust

This project requires rust, please follow the guidelines for installation.

To run the basic rust setup, run:

```
make rust_project_setup
```

Build the project:

```
make build
```

Run tests:

```
make test
```

Format-Check:

```
make format-check
```

(You can run format and check individually)

Clean:

```
make clean
```

## Usage
You can find a basic description of the rust package and how to use it [here](examples/Rust%20library%20overview.ipynb). This is mainly meant for contributors and integrators.

If you want to use the algorithms provided in the library to compile your circuits, you can find a simple user guide for the Python bindings [here](examples/Basic%20python%20usage.ipynb).