# lgo_engine (provisional)

The engine the app runs on (spec sections 8 and 9). It is pure NumPy, with no file or network I/O inside the package, so the same code can run under CPython and in Pyodide. It is built against the CONVENTIONS v1.2 review candidate, which the owner has not yet signed off.

## Layout

| Module | Spec | What it does |
|---|---|---|
| `model.py` | F1.3, F1.8 | `Atom`, `BasisFunction` (s, or p along any unit direction), `System`, with validation |
| `params.py` | F1.5, F2.12 | Parameter lookup from rows shaped like `params/`. Order: override, then pair row, then channel default, then error |
| `hamiltonian.py` | F2.1, F2.2 | s/p Slater–Koster assembly with the CONVENTIONS §1 signs, one formula for any p direction |
| `solve.py` | F2.3 | Eigensolver with residual and orthonormality checks; degeneracy clustering by total span |
| `symmetry.py` | F2.4, F2.5 | Basis representation from group data, projector labels, linear-molecule decomposition |
| `analysis.py` | F2.6–F2.8, F2.10, F2.11 | Composition, bonding character, filling and warnings, HOMO/LUMO, both bond-order definitions |
| `ledger.py` | L1 | Contact ledger whose total equals the eigenvalue |
| `result.py` | F2.9 | `run(system, params, groups)` returns a JSON-serialisable Result |
| `schedule.py` | P1, P2, P5, F2.13 | Construction-slider stages from symmetry pairing ("fix one, flip its partner", then combine), B5 skipping of zero-coupling stages, H(s) and checkpoints |
| `diagram.py` | F3.1–F3.4 | DiagramSpec: symmetry-adapted fragment columns, MO column, the two line layers (dashed allowed, solid mixing ≥ 5 %), slider levels |
| `app_api.py` | — | The two entry points the browser worker calls: `case_list` and `view` |
| `cases.py` | F1.1 | Release 1 geometries from CONVENTIONS §2, with lengths and angles taken from `params/presets.json` |

`load.py`, outside the package, is the only code that reads files: `params/`, `groups/`, and the CONVENTIONS hash.

## Tests

Run `make test-engine`, which `make check` also runs. The suite covers:
- all 21 Release 1 fixtures against `evidence/reference_spectra.json`;
- closed-form Appendix A blocks written independently of `hamiltonian.py`;
- bond-order oracles (A.5) and frontier orbitals;
- the partly-filled-shell warnings;
- ledger totals;
- the spec §12 property tests: zero coupling, MO count, invariance under rotation within a degenerate set and under sign flips, and JSON round trip;
- config rejection;
- the construction schedule and DiagramSpec;
- a mutation check.

## Not yet built

- the C1 equivalent-set detection beyond plain "—" labels (F2.14);
- d and f orbitals;
- the 3D viewer (F3.7) and figure export (section 11);
- Practice and Exam modes, node markers and the rule badges (Release 2 and later).

The browser front end is in `web/` (see `web/README.md`).
