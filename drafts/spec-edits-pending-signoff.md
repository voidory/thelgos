> **Applied 2026-11-06** (D1–D13) to the spec under the owner's standing permission. Kept as the review record.

# Spec edits held until owner sign-off (B6)

These edits to `spec/LGO_Method_Explorer_Functional_Specification.md` are ready to apply. They come from the v1.2 TODO tier D and from this session's review. They are held, not applied, because tier D is reserved until the owner signs off v1.2. Each item says where it comes from and what it depends on.

| # | Where in the spec | Edit | Source | Depends on |
|---|---|---|---|---|
| D1 | §6 coverage table (SF6 row), §5 D2, F1.6 | Remove "LGO splitting from trans F–F contacts about 3 Å apart" as a lecture simplification. In the model, cis contacts (2.21 Å) and trans contacts (3.12 Å) both act, and cis contacts are what split eg from a1g; the SF6 lecture preset equals the model. D2 loses its SF6 switch, or it is redefined. | TODO D1; review point 1 | — |
| D2 | F2.12 | Replace "Couplings involving H … under review" with the approved hydrogen scale, its admissible range (0.4915, 0.8939), and a note that it is fitted to orderings, so those orderings cannot count as green fidelity checks. | TODO D2; review | B2 |
| D3 | §12 verification | Add the diatomic blocks and reference-ordering assertions (see `drafts/CONVENTIONS-diatomics-draft.md`). | TODO D3 | C2 accepted |
| D4 | P-series, F2.13 | Add the outcome of B5: stages whose newly enabled coupling block is zero are skipped, with the note "No ligand–ligand coupling in this model". | TODO D4 | B5 |
| D5 | §6 coverage table, T0 template row | Add N2 as a Release 1 diatomic, and drop "N2 once nitrogen parameters exist". | TODO D5 | B4 = yes (owner) |
| D6 | §7 fidelity table, Orbital-ordering row; §6 diatomics row | Replace "B2 red without mixing, green with the R2 switch" with the actual model result. Under the current scheme, B2, C2 and N2 stay red for orbital ordering with s–p mixing on, and B2 and C2 also for spin state. The "diboron's magnetism via the R2 switch" demonstration (fidelity item B2) cannot work as written. | K14, K16 findings | M1 (decided under delegation) |
| D7 | §6 CH4 row | "Photoelectron-spectrum 6:2 intensity check" is not valid as stated. The 2a1 band (~23 eV) is beyond He I (21.2 eV), and band intensities follow cross-sections, not occupancy. Reword to "two bands, degeneracies 3:1", badged approximate. | Review point 3 | — |
| D8 | §2 rule map, R2 | "Mixing scales roughly as coupling² / energy gap" describes the energy shift; the admixed weight scales as (coupling/gap)². State which one is meant. | Review point 4 | — |
| D9 | §5 D10 | Measure the 100 ms budget warm. Add a separate first-load budget for Pyodide plus the NumPy download. | Review point 6 | — |
| D10 | §14 linear molecules | Move the D∞h/C∞v label lookup into Release 1. It now exists in the checker (K14) and is needed for the Release 1 diatomics. | Review point 2; K14 | — |
| D11 | §12 CrystOD cross-check | Replace the energy-order criterion with the descriptor-order criterion, the contract-basis weight cut-off and the blocking/review split (K7, K8, K11), if the owner accepts them into CONVENTIONS §11. | K7, K8, K11 | Owner accepts K-patches |
| D12 | §7 "Green" definition | Add that a check whose outcome a parameter was fitted to (for example the 0.6 hydrogen scale) cannot be green. | Review | B2 |
| D13 | §6 release plan, Release 1 gate; §12 CrystOD cross-check | CrystOD comparison deferred until the programme is finished (owner decision B1). Release 1 gates on the baseline and engine tests; the CrystOD gate runs as `make check-release`. | B1 | — |
