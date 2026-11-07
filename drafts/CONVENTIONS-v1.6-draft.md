# CONVENTIONS v1.6 — draft for owner sign-off

**Status: draft, not applied.** `CONVENTIONS.md` v1.5 stays frozen until the owner accepts this. Prepared 2026-11-07 by Claude during Release 3 Milestone A.

## One-page summary

The spec's F2.12 says "any element in the scheme's table runs without hand-tuning". Today it does not.

- **On-sites:** the table (`params/onsite.csv`) holds 11 elements, but Miessler's table, already the on-site source (M7), lists H to Ca.
- **Couplings:** every coupling needs its own row in `params/pairs.csv`, and only 33 element pairs have one. All rows use the same four Harrison η and the same hydrogen factor, so a missing row is never a different number, only a missing copy of the same law. The CONVENTIONS default list is empty (§9.2), so a molecule such as CH3Cl (no C–Cl row), HOCl (no O–Cl, H–Cl) or AlH3 (no Al–H) cannot be built.

| # | Change | Effect |
|---|---|---|
| W1 | **Complete the on-site table** from the same Miessler snapshot (`evidence/libretexts_diatomic_MOs.html`): add He, Li, Be, Ne, Na, Mg, Ar, K, Ca. Li, Be, Na, Mg, K and Ca have an s value only, so a molecule needing their p is still "not built yet: element without parameters" (owner decision R3-1). | No existing value changes. |
| W2 | **Universal-law defaults.** `params/channel_defaults.json` gains the Harrison law for the four channels: η ssσ −1.40, spσ +1.84, ppσ +3.24, ppπ −0.81; K 7.62 eV Å²; scale 0.6 when the pair contains H, 1 otherwise, the same rule as every H row. | Every existing pair keeps its row, which takes precedence (§9.2), so no fixture value changes. |
| W3 | **When the defaults apply.** Only to an element pair with no row at all, and only when both elements are in the on-site table. A pair that has rows keeps the current rule that a missing channel is an error: O–H stays without ppπ (K3). An element absent from the table still errors (K16 uses Ne, which W1 adds, so the test moves to an element outside the table, e.g. Kr). | The no-silent-zero guarantee is unchanged. |

Code changes this needs (checker patch K22):
- the default lookup in `lgo_engine/params.py` and `scripts/parameters.py` chooses the H or non-H default by the pair, and applies W3;
- the K16 test element changes;
- a negative test confirms that O–H ppπ still errors.

## What it unlocks

- **T3:** CH3Cl, with C3v and inequivalent ligands. It is the spec's method-fit example (B5).
- **T2:** AlH3.
- **T4:** SiF6²⁻ and PF6⁻, once sourced geometries exist (owner decision R3-3); CCCBDB has few ion geometries, so these may stay "not built yet: no sourced geometry".
- **Packs with Cl partners:** HOCl and others.

Without W1–W3, Release 3's matcher would report "element pair without parameters" for these, although every number they need is already in the scheme.

## Sections touched

- §9.1 table: W1 rows.
- §9.2: replace "the supplied default list is empty" with the W2 defaults and the W3 rule.
- §12: change-log row.

No Appendix B row changes.
