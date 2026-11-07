PYTHON ?= python
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
# [K20] Pin the OpenBLAS kernel so float output is byte-identical on any x86-64 host (AVX2 path).
export OPENBLAS_CORETYPE=Haswell
.PHONY: check-baseline test-engine test-web check check-release
check-baseline:
	$(PYTHON) scripts/check_baseline.py
test-engine:
	$(PYTHON) engine/tests/test_engine.py
test-web:
	node web/tests/web_test.js
check: check-baseline test-engine
	@echo 'CrystOD comparison deferred by owner decision B1 (DECISIONS.md); run make check-release for the release gate.'
check-release: check
	@tmp=$$(mktemp -d); trap 'rm -rf "$$tmp"' EXIT; $(PYTHON) scripts/compare_crystod.py --out "$$tmp" $(if $(CRYSTOD_ZIP),--archive "$(CRYSTOD_ZIP)") $(if $(CRYSTOD_ADAPTER),--adapter "$(CRYSTOD_ADAPTER)")
