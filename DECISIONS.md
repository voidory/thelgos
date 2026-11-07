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
| M5 | Period-3 diatomics Al2–Cl2 added; sources documented | owner request, details delegated |
| M6 | One entry per molecule with a model/lecture toggle | owner |
| M7 | All on-sites from the Miessler teaching table | owner |
| B6 | CONVENTIONS v1.4 signed off and frozen (2026-11-07) | owner |

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

**Decision.** Keep 0.6. Re-confirmed by the owner on 2026-11-07 ("do whatever aligns with Miessler"). Dropping the scale (1.0) would flip two textbook LUMOs: water's would become 2b2* instead of 3a1*, and CH4's t2* instead of a1*. Miessler's figures, as I recall them, show a1* as the LUMO in both, which 0.6 reproduces. This should be checked against the book. The value range that keeps all four orderings is 0.4235–0.6668.

**Why.**
- Any value from 0.4235 to 1.0313 gives the required orderings, so the number is not fine-tuned within that range. (The range was 0.4915–0.8939 before the on-site switch, M7.)
- 0.6 is what every committed regression value, oracle and reviewed number already uses, so choosing it changes nothing that has been checked.
- It leaves BH3's empty Lewis-acid orbital (the LUMO) 2.30 eV clear of the next level (1.41 eV before M7), comfortably separated for teaching.

**Condition attached.** Because 0.6 was chosen to make orderings come out right, those orderings cannot later count as independent "green" evidence in the fidelity bar. This is held spec edit D12.

**Alternatives.**
- The midpoint of the range (about 0.69) is no better justified.
- 1.0, meaning no scale, broke the lecture-water ordering under the earlier on-sites. Under M7 it no longer does (the window now reaches 1.0313), so the scale is no longer strictly needed for these two orderings. See M7 for why 0.6 is kept for now.

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

## M5. Period-3 diatomics: which ones, and from which source (owner request; details delegated)

**Owner request (2026-11-06).** Add the period-3 diatomics.

**Decision.** Add Al2, Si2, P2, S2 and Cl2. Si2, P2 and S2 went in on 2026-11-06. Al2 and Cl2 followed on 2026-11-07, once the owner opened access to chem.libretexts.org.

**Why.**
- *Source for Si, P and S.* The Harrison transcription (`evidence/harrison.txt`, lines 93–95) has Si, P and S, the same source as S already. So Si2, P2 and S2 need no new source. Geometries and ground states come from the NIST CCCBDB pages, snapshotted in `evidence/cccbdb_{Si2,P2,S2}.html` (r_e from Huber & Herzberg 1979).
- *Al and Cl.* The Harrison transcription has no Al or Cl, so they come from the teaching table already used for H and F, now snapshotted in `evidence/libretexts_diatomic_MOs.html`. Values: Al 3s −11.32, 3p −5.98; Cl 3s −25.23, 3p −13.67 eV. The page credits the table to Miessler, Fischer & Tarr (Pearson, 2014).
- *Scale (superseded by M7).* The owner first chose to keep the mix, then switched to one table (M7). The two tables are not on one scale: the teaching table is deeper (Si 3s −15.89 against Harrison's −13.55). Following the existing composite rule (Harrison where available) keeps every earlier result unchanged. The cost is that Al2 and Cl2 sit a few eV deeper than they would on Harrison's scale, so absolute energies along Al2 → Cl2 are not comparable; shapes within each molecule are unaffected. The alternative was to move every element to the teaching table, which would change every existing fixture. **Owner decision (2026-11-07): keep the mix.** The app's reading note on Al2 and Cl2 says this.
- *Al2 spin disagrees.* The model gives a closed shell, with σg(3p) about 1.7 eV below πu, while experiment gives a triplet. This is recorded as a fidelity warning, not tuned (as for B2 under M1). CCCBDB lists Al2's lowest level as ³Σg with degeneracy 4, which is not self-consistent for a ³Σ state. Every candidate ground state is a triplet, though, so the spin comparison stands.
- *Not added.* Na2, Mg2 and Ar2 are not standard MO-diagram examples, and Mg2 and Ar2 are only van der Waals bound.
- *No orbital-ordering reference.* None is asserted for period 3: CONVENTIONS A.6 sources only period 2, and I have no sourced period-3 ordering in hand. Only the ground-state spin is compared. Si2, P2, S2 and Cl2 agree with CCCBDB (Si2 ³Σg, P2 ¹Σg, S2 ³Σg, Cl2 ¹Σg); Al2 does not.
- *What the model gives.* σg(3p) below πu in all five, with Si2 ending in a half-filled πu (triplet).

**Knock-on.**
- `checks/crystod_manifest.json` now lists the three new cases among those the deferred CrystOD comparison must cover (it is derived from the fixture list).
- Checker patches K18 and K19 (CHECKER-CHANGES.md).

## M6. One entry per molecule with a model/lecture toggle (owner decision)

**Owner decision (2026-11-06).** H2O and SF6 appear once each, with a toggle.

**How.**
- Only the app's case list merges them. The engine, checks, evidence, share links and exports keep the two fixture names ("H2O model", "H2O lecture"), so nothing frozen changes.
- The toggle states what the lecture preset changes, worked out from the presets: for H2O, the O 2s–H 1s coupling is switched off. For SF6 it says the two presets are identical (see M4), rather than implying a difference that does not exist.

## M7. One on-site table for every element: Miessler (owner decision)

**Owner decision (2026-11-07).** Take every on-site energy from the Miessler, Fischer & Tarr table as given on LibreTexts (snapshot `evidence/libretexts_diatomic_MOs.html`), because the owner's module uses Miessler. This reverses the "keep the mix" choice recorded under M5 earlier the same day.

**Why the tables differ.**
- Harrison's values are computed atomic term values (as I recall, from Herman–Skillman 1963 calculations), which come out too shallow.
- The teaching values are measured ionisation energies averaged over each configuration's spectroscopic states. I believe they come from Mann, Meek & Allen 2000; the LibreTexts page itself credits only Miessler.
- The teaching values run 1.5–4 eV deeper. The s–p gaps are similar.

**What changed.**
- `params/onsite.csv`: values for B, C, N, O, Si, P and S, and the source column for all rows.
- `evidence/reference_spectra.json`, `checks/release1_energies.csv`, `checks/release1_validation.json`, `checks/crystod_cases.json`, CONVENTIONS Appendix B and §9.1–9.2, and spec F2.12.
- Unchanged (they use no changed element): H2, H3+, neutral H3, six H, four H, F2, Al2, Cl2.

**What did not change, verified before switching.**
- Every checker assertion passes as before.
- Every diatomic's ordering agreement, spin agreement and counting bond order is the same.
- The CO HOMO is still 98 % carbon.
- The engine and app code needed no changes: on-sites are data.
- The SF6 levels shift (S is deeper, so a1g goes from −11.79 to −13.12 and t1u* from +1.43 to +0.64), but its ordering, labels and construction stages are unchanged. The eg level is F-only and does not move.

**Flip audit (owner request, 2026-11-07).** `engine/tools/compare_params.py` compared every case under the two tables (old params from commit 64a2202). It checks level order, electrons, HOMO/LUMO, warnings, the bonding/antibonding star, comparisons with experiment, bond orders, construction stages, the order at every construction checkpoint, and which fragment orbital dominates each MO. 5 changes across 23 cases, none in a finished diagram's level order:
- *BH3.* B 2s moves from 1.07 eV above the H 1s set to 0.44 eV below it. So at the separate-atoms step the order is B 2s < H set, and the lowest a1′ is now mostly B 2s (51 %, was 47 %); a1′* swaps the other way.
- *SF6, model and lecture.* At construction step 1, S 3s (−22.71) drops below the lower three fluorine pair combinations (−21.18). Before, S 3s sat above them at −20.80. Only that intermediate step changes.
- *N2.* 1σg and 2σg are near-even 2s/2p mixes, and their majority character tips across 50 %: 1σg goes from 50 % 2p to 52 % 2s, 2σg the reverse. The cause is Miessler's wider 2s–2p gap (12.38 against 11.57 eV), which means slightly less s–p mixing.

Each of these has a "Miessler vs Harrison energies" reading note in the app (`params/case_notes.json`).

**Couplings are unchanged** (Harrison η, M2). The model now pairs measured-ionisation on-sites with Harrison's couplings. It was already a pedagogical mix; this makes it no less one.

**Hydrogen scale.** The window widened to 0.4235–1.0313, which now includes 1.0. Keeping 0.6 leaves every H-containing result as it would be under the decision already approved (B2), and changes one thing at a time. Whether to drop the scale altogether is a separate owner question.

## Release 2, Milestone A: design choices (delegated, 2026-11-07)

The owner asked for Milestone A first and for D10 to be measured on the owner's phone. The spec fixes each item's purpose and "done when"; these are the details it leaves open, and none needs a CONVENTIONS change.

- **D2 approximation slider.** The path from "lecture" to "model" is linear in the Hamiltonian: H(t) = H_lecture + t·(H_model − H_lecture). Only the couplings that differ move; for water that is the O 2s–H 1s coupling. The ends are exactly the two presets (tested to 1e−12 eV). SF6's presets are identical (M4), so it gets no slider.
- **D4 degenerate-set dial.** It rotates between the first two components of a degenerate set (0–90°). A 3-fold set has more freedom than that, but one rotation is enough to show the invariance: the energy and the set's per-atom weights stay fixed (tested on every degenerate MO and fragment orbital in all 23 cases). For the two-column layout, the separate atoms are degenerate only under the on-site part of H, so their dial and ledger use that part.
- **D5 lone-pair dial.** Water's highest occupied a1 and b1 levels rotate by up to 45°. At 45° they are equal-energy mirror images through the molecular plane. The sum of the two energies and the total density are unchanged (tested). The sketches are drawn side-on, because seen face-on the two mirror images look identical.
- **D6 electron stepper.** The charge is worked out from neutral-atom counting. An atom with its full valence basis gives its valence electrons; an atom represented only by p functions (π-only carbon, an inward ligand p) gives one per function. This gives the expected charges for all 23 fixtures (H3+ and the allyl cation +1, all others 0).
- **F3.6 symmetry-only view.** Levels are spaced evenly by energy rank, with no energy scale. Labels and lines are unchanged (tested). The order itself is still the model's, because symmetry alone does not fix it; the axis says "order only".
- **D10 speed test.** A footer button times every case twice and reports the warm pass: engine time, and engine plus redraw, against the spec's 100 ms budget. The owner runs it on the reference phone and sends the copied result.

## B6. Sign-off (owner only)

Not delegated. On 2026-11-06 the owner gave standing permission for Claude to apply the drafted changes. They are now in `CONVENTIONS.md` v1.3 and the spec, and the drafts in `drafts/` remain as the record. Applying them is not sign-off: v1.3 stays "not frozen" until the owner reviews it.

**Signed off, 2026-11-07.** After reading `drafts/signoff-summary-v1.4.md`, the owner wrote "I sign off v1.4". `CONVENTIONS.md` v1.4 is frozen. The signed text is commit 81cf760, sha256 `ca861a1d6f47f5fe9962eb229a3fbdf4313c22ed2e83083dcfa8fdbe71b98d36`. Any later change needs a new version and the owner's approval. Three points the summary flagged as resting on weaker sources remain open for correction in a later version:
- the Miessler water and CH4 LUMO orderings behind the hydrogen factor's upper limit;
- CH4's second photoelectron band near 23 eV;
- the inconsistent CCCBDB entry for Al2.
