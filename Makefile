PYTHON ?= python
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
.PHONY: check-baseline check
check-baseline:
	$(PYTHON) scripts/check_baseline.py
check: check-baseline
	@tmp=$$(mktemp -d); trap 'rm -rf "$$tmp"' EXIT; $(PYTHON) scripts/compare_crystod.py --out "$$tmp" $(if $(CRYSTOD_ZIP),--archive "$(CRYSTOD_ZIP)") $(if $(CRYSTOD_ADAPTER),--adapter "$(CRYSTOD_ADAPTER)")
