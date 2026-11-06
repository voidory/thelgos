# Decision log

The owner decided B1 and B4 and fixed the ring orientation part of B3. The owner delegated B2, B3 (in-plane detail), B5 and the follow-on modelling questions to Claude, asking that each decision be documented for later review (2026-11-06). Every entry below says what was decided, what it means in plain terms, why, what else could have been done, and how to reverse it. Nothing here edits `CONVENTIONS.md`, which stays the owner's to change. Proposed CONVENTIONS text is in `drafts/`.

Status key: **owner** = decided by the owner; **delegated** = decided by Claude under the owner's delegation, open to review; **owner only** = still needs the owner.

| # | Decision | Status |
|---|---|---|
| B1 | CrystOD comparison deferred until the programme is finished | owner |
| B2 | Hydrogen coupling scale stays 0.6 | delegated |
| B3 | Rings flat, above and below the metal; atoms on ±x in the ring plane | owner (orientation), delegated (in-plane axis) |
| B4 | N2 is in Release 1 | owner |
| B5 | Construction stages with zero coupling are skipped, with a note | delegated |
| M1 | B2, C2 and N2 keep their red orbital-ordering flags; parameters not changed | delegated |
| M2 | Harrison's 1980 coupling coefficients stay | delegated |
| M3 | SF6 F–F couplings stay at the full nearest-neighbour law | delegated |
| M4 | Preset data unchanged; misleading "lecture" tags flagged for the spec, not renamed | delegated |
| B6 | Sign-off of CONVENTIONS v1.2 | owner only |

---

## B1. CrystOD comparison deferred (owner)

**Decision.** The cross-check against the external CrystOD program waits until the programme is finished.

**Effect.**
- `make check` now runs the baseline and the engine tests and passes. It prints that the CrystOD step is deferred.
- The CrystOD gate moved to `make check-release`, which still stops with BLOCKED until the original zip and a real adapter are supplied. That gate's anti-tampering logic is unchanged.

**Why the gate is kept, not removed.** The spec makes CrystOD a Release 1 gate, and the work to compare against it already exists. Moving it to its own target lets day-to-day checking pass without pretending the cross-check happened.

**Reverse.** Point `check` back at the CrystOD step in `Makefile`.

## B2. Hydrogen coupling scale stays 0.6 (delegated)

**Plain terms.** Couplings involving hydrogen use Harrison's formula, which was never fitted to hydrogen. They are multiplied by 0.6 so that the textbook level orderings come out right.

**Decision.** Keep 0.6.

**Why.**
- Any value from 0.4915 to 0.8939 gives the required orderings, so the number is not fine-tuned within that range.
- 0.6 is what every committed regression value, oracle and reviewed number already uses, so choosing it changes nothing that has been checked.
- It leaves BH3's empty Lewis-acid orbital (the LUMO) 1.41 eV clear of the next level, comfortably separated for teaching.

**Condition attached.** Because 0.6 was chosen to make orderings come out right, those orderings cannot later count as independent "green" evidence in the fidelity bar. This is held spec edit D12.

**Alternatives.**
- The midpoint of the range (about 0.69) is no better justified.
- 1.0, meaning no scale, breaks the water ordering.

**Reverse.** Change `scale` on the H rows in `params/pairs.csv` and `hydrogen_scale` in the presets, then regenerate the reference spectra.

## B3. Ring orientation (owner: flat above and below; delegated: in-plane axis)

**Owner's decision.** In the sandwich (uranium) case, each ring lies flat, one above and one below the metal. So the ring planes are perpendicular to the metal–ring axis, the z axis, with the metal at the origin.

**Delegated detail.** Within each ring's plane, one ring atom lies on the +x axis.

**Why.**
- This is the orientation the committed D6h group file was built and checked against Bilbao in. Its C2′ axes pass through atoms, and its C2″ axes through bond midpoints.
- The committed D8h file, used for uranocene, also has a C2′ rotation about x, so the same rule works there.
- Choosing it changes no committed data.

**Follow-on, not yet decided.** Whether the two rings are eclipsed or staggered relative to each other is fixed by the molecule, not by convention. Uranocene is D8h (eclipsed) and bis(benzene)chromium is D6h (eclipsed). This is from my general knowledge; it will be checked against a cited structure when the sandwich templates are built in Release 4.

**Reverse.** The in-plane axis is one line in CONVENTIONS §2 plus the D6h axis assertions in `scripts/generate_groups.py`.

## B4. N2 in Release 1 (owner)

**Effect.**
- Nitrogen on-site energies added: 2s −23.04 eV, 2p −11.47 eV, from the same Harrison transcription as B, C, O and S (`evidence/harrison.txt`).
- Nitrogen pair rows added for every element present, with the same coefficients and rules as the existing rows. N–H uses the 0.6 hydrogen scale.
- N2 fixture added at r_e = 1.0977 Å (NIST CCCBDB, Huber & Herzberg 1979; `evidence/cccbdb_N2.html`).

**Result.**
- Counting bond order 3.
- Closed-shell singlet, matching the reference ground state.
- σg(2p) lies below πu. In real N2, πu lies below σg(2p) because of s–p mixing, so N2 joins B2 and C2 under M1.

**Side effect.** The v1.2 test that a missing element raises an error used N–N as its example. It now uses Ne–Ne (CHECKER-CHANGES K16).

## B5. Construction stages with no coupling (delegated)

**Plain terms.** In CH4, H2O and BH3 the model has no hydrogen–hydrogen coupling. So the "combine the ligands" step of the construction slider would move without anything changing.

**Decision.** Skip any stage whose newly added couplings are all exactly zero, and show the note "No ligand–ligand coupling in this model." SF6 keeps its combine stage, because its F–F couplings are on.

**Why.** The two alternatives each mislead:
- switching H–H coupling on just for that stage shows energies the model never uses;
- keeping a dead stage suggests the ligands genuinely don't interact.

Skipping with a note says exactly what the model assumes.

**Reverse.** A flag in the stage generator (F2.13, not yet built).

## M1. B2, C2 and N2 orbital ordering (delegated)

**Plain terms.** For B2, C2 and N2 the model puts two orbitals in the wrong order compared with experiment, and for B2 and C2 it gets the spin wrong too. O2, F2 and CO come out right.

**Decision.**
- Keep the parameters.
- Show B2, C2 and N2 with red Orbital-ordering chips, plus red Spin-state chips for B2 and C2.
- Reword the spec's plan to demonstrate "B2 red without mixing, green with mixing". This is held spec edit D6.

**Why.**
- The failure is a known limit of this simple model, and showing it is what the fidelity bar is for.
- Changing parameters to fix these molecules would be fitting to the answer, which CONVENTIONS and the v1.2 review both rule out.
- It would also move every other molecule's numbers.

**Alternative left open.** A deliberate change of parameter scheme, as a choice about the model rather than about these molecules (M2), would need a full regression rerun.

## M2. Harrison coupling coefficients (delegated)

**Decision.** Keep Harrison's 1980 book values (ssσ −1.40, spσ +1.84, ppσ +3.24, ppπ −0.81).

**Why.**
- They are the values cited in the spec (F2.12) and used by every reviewed number.
- `evidence/harrison.txt` also lists a second set, which it attributes to Shi & Papaconstantopoulos (2004) as Harrison's values (−1.32, +1.42, +2.22, −0.63). That set has no stated advantage for these molecules, and switching would move every result.

**Documentation owed.** CONVENTIONS §9.2 should name which set is used and why. Draft text is in `drafts/CONVENTIONS-diatomics-draft.md`.

## M3. SF6 F–F couplings (delegated)

**Decision.** Keep F–F couplings at the full Harrison law, including the non-bonded cis pairs 2.21 Å apart.

**Why.**
- v1.1 deliberately set SF6 lecture equal to the model with these couplings, and it was reviewed in that form.
- Weakening non-bonded contacts would be a new parameter with no source.

**Concern recorded.** These couplings spread the fluorine combinations over about 19 eV, larger than the lecture picture suggests. If a source for non-bonded damping turns up, revisit.

## M4. "Lecture" tags on presets identical to the model (delegated)

**Decision.** Leave `params/presets.json` unchanged.

**Why.** `SF6_lecture` and `sixH_lecture` are tagged "lecture" but no slide has been checked against them. Renaming the ids would ripple through checks and evidence for a labelling issue. Instead, held spec edit D1 removes the claim that SF6 has a separate lecture simplification, and no exported label will claim slide fidelity until the slides are checked.

## B6. Sign-off (owner only)

Not delegated, and not something Claude can do. When you are ready, the drafts in `drafts/` show exactly what would change in CONVENTIONS.md and the spec.
