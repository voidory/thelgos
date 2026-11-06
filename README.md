# LGO conventions v1.2 review candidate

This is an actual generated audit package, not a claim that a file has been produced elsewhere. `CONVENTIONS.md` is the primary candidate. **Not frozen. The owner has not approved v1.2.** The input TODO is preserved in `evidence/CONVENTIONS-v1.2-TODO.md`.

## What changed

- Geometry-first Bilbao mapping replaces character/label-matched comparisons. Raw MolSym D2h fails with character error 2; the documented B2/B3 g/u label repair passes. Cartesian assertions cover D2h/C2v and the +x C2-prime axes in D3h/D6h. Ten group files include exact-value snapping and residuals.
- Actual on-site/pair/preset files drive both Hamiltonian and analytic calculations. A mutated O–H pair coefficient in a temporary copy fails the committed-spectrum guard. Twelve fixtures execute, including neutral H3 and warning assertions for partially filled degenerate shells.
- `make check-baseline` regenerates into a temporary directory, diffs groups and checks, checks all hashes, and asserts no committed data changed.
- Pair lookup/provenance enums, scale-carrier semantics, error-on-missing behaviour, preset JSON schema and the full unscaled coupling table are specified and checked.
- The complete checker environment is hash-locked; a clean hash-enforced offline wheel installation and `pip check` completed. Host is CPython 3.12.8, not the production Python 3.12.7 browser environment.
- A computed CrystOD gate cannot be satisfied by hand-editing the manifest: it recomputes the status. The original archive and a verified native adapter are unavailable; the gate is BLOCKED. Mock protocol tests are not external-checker results.

## Run the checks

Use Linux x86_64 / CPython 3.12. The tested host is 3.12.8.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install --require-hashes -r requirements-checker.lock
make check-baseline
make check
```

Expected: baseline passes; full `make check` fails with the explicit missing-CrystOD reason. The Make recipe returns nonzero. It also reports missing diatomic coverage. Actual browser/UI release tests are separate and have not been run.

To identify and later compare a real supplied archive:

```sh
make check CRYSTOD_ZIP=/absolute/path/original.zip \
  CRYSTOD_ADAPTER=/absolute/path/verified_native_adapter.py
```

The adapter interface is documented in `scripts/compare_crystod.py` and `CONVENTIONS.md` §11. It is not a native CrystOD implementation; that implementation needs the actual archive/API. Do not create a fake adapter from the expected spectra. Outputs must disclose basis differences and axis/irrep maps, and assign within-irrep branch identities independently of their energy order. More diatomic fixtures must be added before all required comparisons can pass.

For fresh source retrieval rather than the offline baseline:

```sh
python scripts/fetch_bilbao.py --out /tmp/new-bilbao
python scripts/generate_groups.py --evidence /tmp/new-bilbao --out /tmp/new-groups
```

The baseline uses the packaged freshly retrieved snapshots to avoid website availability dependence.

## Tier-A disposition

| Item | Disposition | Evidence |
|---|---|---|
| A1 | Corrected; old D2h PASS retracted | `checks/bilbao_D2h.json`, group labels, axis assertions |
| A2 | Implemented and executed | `parameter_mutation.json`, `release1_validation.json` |
| A3 | Implemented and executed | `scripts/check_baseline.py`, `SHA256SUMS` |
| A4 | Missing-archive gate and protocol implemented; real native adapter/execution explicitly unavailable | `crystod_manifest.json`, `crystod_protocol_tests.json`; never claim real CrystOD PASS |
| A5 | Column specification restored; unavailable legacy A.7 explicitly removed/replaced with reason | §12 and replacement Appendix A.7; exact old content not recovered |
| A6 | Typed contract and lookup implemented | `parameters.py`, `pairs.csv`, §9.2 |
| A7 | Schema present and validated; missing-field negative test | `presets.schema.json`, `preset_schema_validation.json` |
| A8 | Exact-value snapping, residual assertion and report | `group_validation.json` |
| A9 | Neutral-H3, four-H and six-H warnings asserted | `release1_validation.json`, energies CSV |
| A10 | Full hashed dependency graph; clean install executed | lock, `evidence/hash_install.txt` |
| A11 | Completed from actual pair data | `unscaled_reference_couplings.csv` |
| A12 | Default orientation specified and asserted, awaiting the owner's B3 decision | §2/§5, D6h group report |
| A13 | Header/review/§11 list CrystOD, diatomics and browser/engine gates | primary document |

## Owner decisions and later gates

B2=0.6, B3=benzene atoms on ±x, B5=skip zero-coupling stages are **defaults pending approval**. B1 still needs the original zip. B4 is unresolved: N2 stays excluded; N parameter values are not added. B6 sign-off has not happened. C2 diatomic source geometries/explicit blocks/regressions, C1 real checker and C4 browser tests are still open. Tier-D spec edits were not applied.

No lecture-slide alias or slide-calibrated numeric value is claimed without slides. Later d-channel parameters remain out of scope. This package corrects evidence mechanics; it is not a newly completed application engine or a release authorisation.
