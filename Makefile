# Variables
CARGO := cargo

python_setup:
	pyenv install -s $(shell cat .python-base-version)
	pyenv virtualenv-delete -f $(shell cat .python-version) || true
	pyenv virtualenv $(shell cat .python-base-version) $(shell cat .python-version)
	pyenv version

project_setup:
	pip-compile --allow-unsafe --no-header --no-annotate --output-file=./requirements.txt requirements.in
	pip install -r requirements.txt
	pre-commit install --hook-type pre-commit --hook-type prepare-commit-msg --hook-type commit-msg

build:
	@$(CARGO) build --verbose

test:
	@$(CARGO) test --verbose

format:
	@$(CARGO) fmt

check:
	@$(CARGO) fmt -- --check
	@$(CARGO) clippy -- -D warnings

clean:
	@$(CARGO) clean
