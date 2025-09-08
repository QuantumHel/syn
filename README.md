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