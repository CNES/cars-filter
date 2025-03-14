# Autodocumented Makefile
# see: https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
# Dependencies : python3 venv
# Some Makefile global variables can be set in make command line
# Recall: .PHONY  defines special targets not associated with files

############### GLOBAL VARIABLES ######################
.DEFAULT_GOAL := help
# Set shell to BASH
SHELL := /bin/bash

# Set Virtualenv directory name
# Exemple: VENV="other-venv/" make install
ifndef VENV
	VENV = "venv"
endif

# Software version from setup.py and setuptools_scm
VERSION = $(shell python3 -c 'from cars_filter import __version__; print(__version__)')
VERSION_MIN = $(shell echo ${VERSION} | cut -d . -f 1,2,3)

# Browser definition
define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"

# Python global variables definition
PYTHON_VERSION_MIN = 3.9

PYTHON=$(shell command -v python3)

PYTHON_VERSION_CUR=$(shell $(PYTHON) -c 'import sys; print("%d.%d"% sys.version_info[0:2])')
PYTHON_VERSION_OK=$(shell $(PYTHON) -c 'import sys; cur_ver = sys.version_info[0:2]; min_ver = tuple(map(int, "$(PYTHON_VERSION_MIN)".split("."))); print(int(cur_ver >= min_ver))')

############### Check python version supported ############

ifeq (, $(PYTHON))
    $(error "PYTHON=$(PYTHON) not found in $(PATH)")
endif

ifeq ($(PYTHON_VERSION_OK), 0)
    $(error "Requires python version >= $(PYTHON_VERSION_MIN). Current version is $(PYTHON_VERSION_CUR)")
endif


################ MAKE targets by sections ######################

.PHONY: help
help: ## this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' | sort

## Install section

.PHONY: git
git: ## init local git repository if not present
	@test -d .git/ || git init .

.PHONY: venv
venv: ## create virtualenv in "venv" dir if not exists
	@test -d ${VENV} || python3 -m venv ${VENV}
	@${VENV}/bin/python -m pip install --upgrade pip setuptools # no check to upgrade each time
	@touch ${VENV}/bin/activate

.PHONY: install
install: venv git  ## install the package in dev mode in virtualenv
	${VENV}/bin/python -m pip install -e .[dev]
	@test -f .git/hooks/pre-commit || echo "Install pre-commit"
	@test -f .git/hooks/pre-commit || ${VENV}/bin/pre-commit install -t pre-commit
	@chmod +x ${VENV}/bin/register-python-argcomplete
	@echo "cars-filter ${VERSION} installed in dev mode in virtualenv ${VENV} with documentation"
	@echo " cars-filter venv usage : source ${VENV}/bin/activate"

## Test section

.PHONY: test
test: install ## run tests and coverage quickly with the default Python (source venv before)
	@${VENV}/bin/pytest -o log_cli=true --cov-config=.coveragerc --cov --cov-report=term-missing

.PHONY: test-all
test-all: install ## run tests on every Python version with tox (source venv before)
	@${VENV}/bin/tox -r -p auto  ## recreate venv (-r) and parallel mode (-p auto)

.PHONY: coverage
coverage: install ## check code coverage quickly with the default Python
	@${VENV}/bin/coverage run --source cars_filter -m pytest
	@${VENV}/bin/coverage report -m
	@${VENV}/bin/coverage html
	$(BROWSER) htmlcov/index.html

## Code quality, linting section

### Format with isort and black

.PHONY: format
format: install format/isort format/black  ## run black and isort formatting (depends install)

.PHONY: format/isort
format/isort: install  ## run isort formatting (depends install)
	@echo "+ $@"
	@${VENV}/bin/isort cars_filter tests

.PHONY: format/black
format/black: install  ## run black formatting (depends install)
	@echo "+ $@"
	@${VENV}/bin/black cars_filter tests

### Check code quality and linting : isort, black, flake8, pylint

.PHONY: lint
lint: install lint/isort lint/black lint/flake8 lint/pylint lint/mypy ## check code quality and linting (source venv before)

.PHONY: lint/isort
lint/isort: ## check imports style with isort
	@echo "+ $@"
	@${VENV}/bin/isort --check cars_filter tests

.PHONY: lint/black
lint/black: ## check global style with black
	@echo "+ $@"
	@${VENV}/bin/black --check cars_filter tests

.PHONY: lint/flake8
lint/flake8: ## check linting with flake8
	@echo "+ $@"
	@${VENV}/bin/flake8 cars_filter tests

.PHONY: lint/pylint
lint/pylint: ## check linting with pylint
	@echo "+ $@"
	@set -o pipefail; ${VENV}/bin/pylint cars_filter tests --rcfile=.pylintrc --output-format=parseable | tee pylint-report.txt # pipefail to propagate pylint exit code in bash

.PHONY: lint/mypy
lint/mypy: ## check linting type hints with mypy
	@echo "+ $@"
	@${VENV}/bin/mypy cars_filter tests

## Release section

.PHONY: dist
dist: clean install ## clean, install, builds source and wheel package
	@${VENV}/bin/python -m pip install --upgrade build
	@${VENV}/bin/python -m build
	ls -l dist

.PHONY: release
release: dist ## package and upload a release
	@${VENV}/bin/twine check dist/*
	@${VENV}/bin/twine upload dist/* --verbose ##  update your .pypirc accordingly

## Clean section

.PHONY: clean
clean: clean-venv clean-build clean-precommit clean-pyc clean-test clean-lint ## clean all targets

.PHONY: clean-venv
clean-venv: ## clean venv
	@echo "+ $@"
	@rm -rf ${VENV}

.PHONY: clean-build
clean-build: ## remove build artifacts
	@echo "+ $@"
	@rm -f *.so
	@rm -fr build/
	@rm -fr dist/
	@rm -fr .eggs/
	@find . -name '*.egg-info' -exec rm -fr {} +
	@find . -name '*.egg' -exec rm -f {} +

.PHONY: clean-precommit
clean-precommit: ## clean precommit hooks in .git/hooks
	@rm -f .git/hooks/pre-commit
	@rm -f .git/hooks/pre-push

.PHONY: clean-pyc
clean-pyc: ## remove Python file artifacts
	@echo "+ $@"
	@find . -type f -name "*.py[co]" -exec rm -fr {} +
	@find . -type d -name "__pycache__" -exec rm -fr {} +
	@find . -name '*~' -exec rm -fr {} +

.PHONY: clean-test
clean-test: ## remove test, logging and coverage artifacts
	@echo "+ $@"
	@rm -fr .tox/
	@rm -f .coverage
	@rm -rf .coverage.*
	@rm -rf coverage.xml
	@rm -fr htmlcov/
	@rm -fr .pytest_cache
	@rm -f pytest-report.xml
	@rm -f debug.log

.PHONY: clean-lint
clean-lint: ## remove linting artifacts
	@echo "+ $@"
	@rm -f pylint-report.txt
	@rm -f pylint-report.xml
	@rm -rf .mypy_cache/
