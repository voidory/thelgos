# CONVENTIONS.md: what's left before Release 1 can start

Reviewed: CONVENTIONS.md v1.1, the "corrected review candidate" dated 2026-11-07, together with its bundle (scripts, params, groups, checks).
Reviewer: Claude, for the owner. Change authority stays with the owner. The agent reports any disagreement here; it does not edit a frozen version.

## How v1.1 stands

I recomputed all of the physics by hand: the sign table, the angular vectors, every SALC coupling in Appendix A, the SF6 cis/trans and two-axis split, the six-H and four-H spectra, the bond orders, the electron counts, and all of Appendix B. All of it is correct.

The remaining problems are in three areas: how the evidence is produced, the D2h group file, and Release 1 coverage, which is larger than the header says. The work below is in four tiers:

- **A.** Agent fixes needed before the owner signs off v1.1.
- **B.** Decisions only the owner can make.
- **C.** Release 1 gates that remain after sign-off.
- **D.** Spec-doc updates I make once v1.1 is approved.

---

## A. Agent fixes needed for v1.1 sign-off (v1.2 candidate)

### A1. D2h group data is wrong *(release-critical)*

**Problem.**
- `bilbao_D2h.json` maps Bilbao's `2y` column to MolSym's `C_2'(0)`. That operation's matrix is diag(1,−1,−1), which is a rotation about **x**.
- Bilbao's `mx` is mapped to `sigma_v(0)`, which actually reflects **y**.
- Because of these two swaps, `group_validation.json` reports "PASS" for D2h when it should fail. The table in §5, which says "PASS for all five", is therefore wrong.
- The swap also flips B2u and B3u, so the x and y labels are exchanged.

**Why it matters.** D2h is the subgroup used for **every diatomic** (§2, Appendix B, H2). With the swap, the πx/πy labels come out wrong in Release 1.

**Fix.**
- Map MolSym classes to Bilbao columns by **operation geometry**: rotation axis, or mirror normal, read from the matrix. Do not map by representative label.
- Regenerate `groups/D2h.json` and assert x = B3u, y = B2u, z = B1u.
- Assert the axis labels for C2v (x = b1, y = b2, z = a1) and D3h (+x on a C2′ axis) the same way.
- Correct the claims in the review record (line 13) and in §5.

### A2. Make the checks read the data files, not hardcoded constants

- `check_release1.py` must load `params/onsite.csv`, `params/pairs.csv` and `params/presets.json`. At the moment it hardcodes their values, so the evidence shows that the *script's own* constants are consistent, not that the *bundle's data* is.
- Add a test that corrupts one value in a copy of `pairs.csv` and confirms the check then fails.

### A3. Make `make check-baseline` non-mutating and reproducible

- Generate into a temporary directory, `diff -r` the output against the committed `groups/` and `checks/`, then run `sha256sum -c SHA256SUMS`.
- The check must fail on any difference and must never overwrite committed files.

### A4. Write the CrystOD manifest with a script

- Add `scripts/compare_crystod.py`. It runs the allowed comparisons from §11: point group, degeneracy pattern, labels after axis mapping, and order within each irrep.
- For each case it records the geometry, the basis differences and the axis map.
- It writes `checks/crystod_manifest.json` with a computed status. A hand-set PASS must be impossible.
- Until the zip arrives, the script exits with BLOCKED, so `make check` fails for the right reason.

### A5. Restore or explicitly drop content lost in the rebuild

v1.1 was rebuilt rather than edited. Two things from v1.0 are now gone without a change-log entry:
- the `pairs.csv` column specification;
- Appendix A.7.

Either restore them, or record "removed, because …" in §12.

### A6. Specify `pairs.csv` explicitly

§9.2 says each row "retains source, units, distance law, scale, provenance kind and override reason", but the file doesn't fully support this. Add the following.

**Columns:**
- an explicit `s_carrier` column, saying which element triggers the 0.6 scale;
- `distance_law` as an enumerated value (`harrison_inverse_square`, `fixed`, `off`), not free text;
- `provenance_kind` as an enumerated value (`harrison`, `harrison_extrapolated_H`, `pedagogical_override`, `lecture_decoupling`).

**Lookup order:** state the precedence as preset override → element-pair row → channel default, and say what happens when no row matches. The answer must be "error", never a silent zero (§9.2: "absent channels are not fitted zeros").

### A7. Add a schema for `presets.json`

- Add `params/presets.schema.json` and validate against it in the checks.
- At minimum the schema must cover:
  - which edges are on or off;
  - which couplings are forced to zero (the water O2s–H case);
  - the H scale;
  - geometry overrides;
  - provenance tag.

### A8. Round exact group values

- Snap matrix entries and characters that are exact (0, ±1, ±½, ±√3/2, …) to their exact values before writing the JSON.
- Assert that the residual from rounding is below 1e−12.

### A9. Test that partly filled degenerate shells raise a warning

**Problem.** §7 requires a warning whenever a degenerate set is partly filled. Two fixtures currently fill such sets silently:
- six-H puts 4 electrons in its T1u;
- four-H puts 2 electrons in its T2.

**Fix.**
- Add **neutral H3** (3 electrons, so e′¹) as a fixture that asserts the Jahn–Teller warning and the "occupation not unique" warning both fire.
- Flag the six-H and four-H diagnostic fillings the same way.
- Use O2 and B2 later as real-molecule cases (see C2).

### A10. Commit a hash-pinned checker lock file

- `checker_environment.txt` lists 45 packages, but `requirements-checker.txt` pins only 6.
- Commit `requirements-checker.lock`, a complete list with `--hash=` entries, and install with `--require-hashes`.

### A11. Complete `unscaled_reference_couplings.csv`

- Add the value that justifies the 0.6 scale: O–H spσ at 0.958 Å = 15.277130 eV.
- Add every other H-containing coupling actually used:
  - O–H ssσ;
  - C–H ssσ and spσ;
  - B–H ssσ and spσ;
  - the six-H trans coupling at 2 Å;
  - the four-H edge.
- Rename the `R_A` column to `d_angstrom`.
- Write each channel in full (`ss_sigma`, `pp_pi`, and so on). The allyl row currently says `pp`, which is ambiguous.

### A12. Fix the D6h convention before D6h is used

- Write benzene's orientation into §2: which axis the atoms lie on.
- State which class is C2′ (through atoms) and which is C2″ (through bond midpoints).
- Assert the mapping as in A1.
- This depends on the owner's decision B3.

### A13. Make the header and review record match reality

The header names only CrystOD as blocking Release 1. Line 77 also concedes that the diatomic fixtures are missing; see C2. List both blockers in the header and in §11.

---

## B. The owner's decisions

| # | Decision | Default if you don't choose | Notes |
|---|---|---|---|
| B1 | **Supply the CrystOD zip** | — | This is the hard blocker. It needs the original bytes so the archive hash can be recorded. |
| B2 | **Approve the hydrogen scale 0.6** | 0.6 | Any value in the admissible range (0.4915, 0.8939) keeps both strict ordering constraints. 0.6 gives BH3 a 1.41 eV gap between the pz LUMO and the a1′ antibonding level. |
| B3 | **Benzene orientation** | Atoms on the ±x axes; C2′ through atoms | Needed only from Release 2 (π ring). Unblocks A12. |
| B4 | **Is N2 in scope?** | Ask: is N2 in your tutorial or lecture set? | v1.1 excludes nitrogen. N2 is the standard example of s–p mixing in diatomics (σg above πu). If the course uses it, add N on-sites from the same Harrison transcription (as I recall, N 2s ≈ −23.04 eV and 2p ≈ −11.47 eV; *verify against the source before use*). |
| B5 | **Construction stages when H–H coupling is off** | Skip stages whose coupling is zero, and show a "no ligand–ligand coupling in this model" note | See below the table. |
| B6 | **Sign off v1.1** (as v1.2) | — | Once tier A is closed. |

**Why B5 matters.** CH4, H2O and BH3 have their H–H edges turned off (§9.2), so any "combine the ligands" stage in the construction slider has zero coupling. The CH4 stage count from the pairing generator is K = 3. If any of those stages is a ligand–ligand stage, the slider will move without anything changing. SF6 does not have this problem, because its F–F edges are on.

The alternatives are:
- turn H–H coupling on for the ligand stage only, which is inconsistent with the energies the app reports;
- leave the dead stage in, which teaches that the ligands in these molecules don't interact.

This choice decides which stages the app shows, so it should be settled before Release 1 is built.

---

## C. Release 1 gates still open after v1.1 sign-off

### C1. Run CrystOD
Install it from the zip, record its version and SHA-256 in §11, run `compare_crystod.py` on all of the Release 1 cases, and log every disagreement for the owner. Where CrystOD and the analytic oracles disagree, the oracles win.

### C2. Diatomic fixtures *(currently missing; §2 line 77)*

The spec puts diatomics in Release 1, but none are built.

**For each of B2, C2, O2, F2 and CO** (and N2 if B4 says yes):
1. Add a geometry with a cited source (NIST CCCBDB r_e).
2. Add an Appendix A oracle as an explicit block. For homonuclear molecules use D2h subgroup blocks (σg/σu from 2s and 2pz, πu/πg from 2px and 2py); CO is C2v.
3. Add an Appendix B regression row.

**Assert these orderings, or flag them if the model disagrees:**
- B2 and C2: πu below σg(2p). These are s–p mixed.
- O2 and F2: σg(2p) below πu.

If the Harrison law gets one of these wrong, the fidelity bar's *Orbital ordering* category must show it. Do not tune parameters to make it pass.

**Expected warnings:**
- O2 has πg², which must raise the spin-state and partial-occupation warnings. This is the real-molecule test for A9.
- B2 has πu², which must raise the same warnings.

**CO checks:** the counting bond order is 3, and the HOMO is σ, weighted towards carbon.

### C3. Group data for Release 1, after A1
Run a fresh Bilbao comparison for D2h, then re-run all Release 1 fixtures. Their labels depend on the group files.

### C4. Engine and UI tests
Appendix B checks Hamiltonian assembly only. Release 1 still needs:
- an executed slider check, comparing the H(s) end-points with the oracle energies;
- an executed ledger check (the counting rule);
- an executed Watch-mode check.

Run these in Pyodide 0.27.7 in the browser, not only under CPython.

### C5. Lecture slide aliases
These don't block the release, but the alias data needs the actual slides. Add it as source-backed `lecture_alias` entries (§6).

---

## D. Spec-doc updates (I make these once v1.1 is approved)

1. Remove the "SF6 lecture simplification" wording from the coverage table, from D2 and from F1.6. In v1.1, SF6 lecture is identical to the model, with cis and trans F–F couplings both nonzero.
2. In F2.12, replace "H couplings under review" with the approved hydrogen scale and its admissible range.
3. Add the diatomic ordering assertions (C2) to the verification section.
4. Add the outcome of B5 to the slider spec (P-series and the stage generator).
5. If B4 is yes, add N2 to the coverage table.

---

## Shortest path to starting Release 1

1. The agent closes A1–A13 and sends a v1.2 candidate, while the owner answers B2–B5.
2. The owner supplies the zip (B1) and signs off (B6).
3. The agent completes C2 (diatomics) and C1 (CrystOD run).
4. I update the spec (D).
5. Release 1 build begins.

Most of tier A is mechanical. A1 and C2 are the only items that can change results: A1 changes labels, and C2 may change what the fidelity bar shows for orbital ordering.
