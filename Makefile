PYTHON ?= python
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
.PHONY: check-baseline test-engine check check-release
check-baseline:
	$(PYTHON) scripts/check_baseline.py
test-engine:
	$(PYTHON) engine/tests/test_engine.py
check: check-baseline test-engine
	@echo 'CrystOD comparison deferred by owner decision B1 (DECISIONS.md); run make check-release for the release gate.'
check-release: check
	@tmp=$$(mktemp -d); trap 'rm -rf "$$tmp"' EXIT; $(PYTHON) scripts/compare_crystod.py --out "$$tmp" $(if $(CRYSTOD_ZIP),--archive "$(CRYSTOD_ZIP)") $(if $(CRYSTOD_ADAPTER),--adapter "$(CRYSTOD_ADAPTER)")
