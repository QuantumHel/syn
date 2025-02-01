# Variables
CARGO := cargo

python_setup:
	pyenv install -s $(shell cat .python-base-version)
	pyenv virtualenv-delete -f $(shell cat .python-version) || true
	pyenv virtualenv $(shell cat .python-base-version) $(shell cat .python-version)
	pyenv version

python_project_setup:
	pip install -r requirements.txt
	pre-commit install --hook-type pre-commit --hook-type prepare-commit-msg --hook-type commit-msg

rust_project_setup:
	rustup show

build:
	@$(CARGO) build --verbose

test:
	@$(CARGO) test --verbose

format-check:
	@$(CARGO) fmt -- --check

format:
	@$(CARGO) fmt

check:
	@$(CARGO) fmt -- --check
	@$(CARGO) clippy

clean:
	@$(CARGO) clean
