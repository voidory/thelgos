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

## M8. Exam and Practice modes check against the Miessler levels (owner decision, 2026-11-07)

**Decision.** Exam and Practice modes (Milestone B) check a student's answers against the app's own levels, which come from the Miessler table (M7). They do not take level orders from a question's hint.

**Why.** Only one real exam question has been seen so far. That is too few to design an alternative convention around.

**Known consequence.** That question asks about NH and hints that H 1s lies slightly above N 2p. Miessler puts it slightly below (−13.61 vs −13.18 eV). Miessler's value is averaged over the whole 2p³ configuration, whereas the hint follows N's ground state (first IE 14.53 eV). The expected answer is the same under both: one σ bond, a nonbonding π pair on N, and a triplet. Only the drawn position of H 1s relative to N 2p differs. If NH is added later (it needs CONVENTIONS v1.5), its reading note should explain the difference. Revisit this decision if more exam questions show the lecturer's ordering differing from Miessler.

## Release 2, Milestone B: Practice and Exam modes, design choices (delegated within the owner's defaults, 2026-11-07)

The owner approved these defaults: the current molecules only; two question types modelled on the owner's past paper (Q2(a) "draw the MO diagram, give the spin multiplicity" and Q2(d) "give and explain the bond order"); Practice checks step by step, Exam at the end; questions like Q2(b) and (c) left out; answers checked against the Miessler levels (M8).

- **Where checking happens.** Questions and checking run in the engine (`engine/lgo_engine/exam.py`), like everything else.
- **"Order the atomic orbitals" uses the on-site energies only.** SF6's ligand group orbitals spread over about 19 eV because of the model's F–F contacts (M3), where a textbook puts them all near F 2p. Grading LGO order would penalise the textbook answer. The paper's hint was also about atomic orbitals.
- **"Which orbitals combine."**
  - *Central-atom cases:* the symmetry-adapted central orbitals and the LGOs. Two orbitals combine when they share an irrep and are coupled in the chosen preset, so in lecture water O 2s correctly stays nonbonding.
  - *Diatomics:* each atom's s, p along the bond (σ) and perpendicular p (π). A cross s–p pairing on the same σ symmetry is accepted but not required (s–p mixing).
  - *Molecules with no central atom:* no such step.
- **MO order** is asked as a sequence of symmetry labels, so equal labels cannot be told apart and do not need to be.
- **Bonding character.** The model's label is accepted. "Nonbonding" is also accepted when an orbital's contacts sum to under 0.5 eV, which covers CO's HOMO, the textbook "carbon lone pair". Bond orders are unaffected.
- **Spin multiplicity** follows Hund's rule from the filled levels. The engine itself still flags open shells as "spin not determined".
- **Left out of the quizzes:** the two hydrogen-only diagnostic sets (six H Oh, four H Td), which are not molecules.
- **Water model vs lecture.** The model counts 3a1 as bonding and gives water a bond order of 1.5; the lecture preset gives the textbook 1. The quiz uses whichever preset is selected. This is model behaviour, not a quiz bug.

## Release 2, Milestone C: node markers, rule checks, phase painting, design choices (delegated, 2026-11-07)

- **Where they apply.** Node markers (L5), the rule check (L4) and phase painting (D3) work on "one orbital per atom" sets with equal on-site energies, the setting R6 and R9 are stated for:
  - the ligand set of a central-atom molecule (its LGO side);
  - a whole molecule with one orbital per atom (H2, H3+, neutral H3, allyl π, the hydrogen-only sets).

  Diatomics with s and p, and the MOs of central-atom molecules, mix different orbital types on one atom. Counting nodes there is not well defined, so these features are not offered.
- **Node count (L5).** The chain or ring is read from which atoms are coupled.
  - *Chain:* the number of sign changes along it, where a zero coefficient counts once (a node through that atom).
  - *Ring:* the number of nodal planes, which is sign changes divided by 2.
  - *3D sets:* markers only, with no count.

  Node markers link to the ledger: between-atom nodes to out-of-phase rows, through-atom nodes to the zero-weight rows, which the ledger now shows (it used to hide zero rows).
- **R6 (L4).** The prediction is ε + (2/N) Σ ± t over the contacts of the nonzero atoms. When every contact is the same size, this is the lecture's ε + (2t/N)(n_in − n_out). It is exact exactly when all nonzero coefficients have equal size. The stated reasons when it is not:
  - unequal coefficients, shown scaled to the smallest one, e.g. (1, 1.41, 1) for allyl ψ1 and (2, 2, −1, −1, −1, −1) for octahedral e_g;
  - unequal contacts (neighbouring vs opposite atoms).
- **R9 (L4).** Ordering by node count is compared with the exact ordering, badged "exact" for chains and "exact up to degenerate pairs" for rings. 3D sets are "heuristic", and sets with no contacts are "not applied".
- **Default component of a degenerate set.** The first component is now the projection of the basis function with the largest weight, so it "points at" an atom. This gives the textbook patterns, (2, −1, −1) for e′ and the d_z²-like e_g, instead of the eigen-solver's arbitrary (and machine-dependent) choice. The D4 dial rotates from there.
- **Painting (D3).** Clicking an atom or one of its lobes cycles + → − → 0. The painted combination's energy, ledger, nodes and rule check follow. Its decomposition onto the set's levels (LGOs, or MOs for a whole-molecule set) shows "a true LGO/MO" when one level holds all the weight.
- **Regression fixed.** Since the p functions on an atom were merged into one oriented lobe (design pass), a p lobe pointing at the viewer was always drawn red, because the merged size is positive and the sign is in the direction. So the allyl π orbitals lost their colours. The face-on lobe is now coloured by whether its positive end faces the viewer, and a browser check guards it.

## Release 2, Milestone D: morph paths and Walsh plots, design choices (delegated, 2026-11-07)

- **Paths built** (`engine/lgo_engine/morph.py`). Each passes through a worked case and reproduces it exactly there (tested to 1e−10 eV):
  - water bending 90°–180°, for both presets;
  - CH4 flattened from tetrahedral to square planar;
  - the allyl chain closed into a ring;
  - H2 + H → linear H3 → triangular H3+;
  - the six-fluorine set → SF6, with the sulfur switched on;
  - a central-orbital energy sweep, for the s or p orbital of H2O, CH4, BH3 and SF6.
- **Not built yet.** BH3 → NH3 needs C3v group data and an NH3 fixture, so it waits for CONVENTIONS v1.5. The spec's literal "six H → SF6" would also change the ligand element (H s to F p); the six-fluorine → SF6 path keeps the ligands and switches on the central atom instead.
- **Geometry paths and coupling paths.** Water and CH4 change geometry; the couplings follow from the same distance law. Allyl, H3 and SF6 switch contacts on, because in a tight-binding model bringing atoms into contact means switching their coupling on. Allyl and H3 use the full neighbour coupling at the end, so the ring and triangle endpoints are the textbook ε + 2t, ε − t, ε − t patterns.
- **Labels along a path.** Each path is labelled in the symmetry it keeps throughout: water C2v, CH4 D2d, allyl C2v, SF6 Oh. The H3 assembly keeps none, so its frames are unlabelled (C1) except at the H3+ point.
  - *Walsh plots:* the Hamiltonian is block-diagonalised by irrep, so different irreps may cross while the same irrep never does. Where there is no symmetry, levels are followed by overlap.
  - *Exact degeneracies:* where a path produces one the group does not require (linear water's πu = a1 + b1 in C2v), the level is labelled "a1+b1" instead of being rejected.
- **A morph frame is a finished diagram.** The construction slider and timeline are hidden during a morph. The Lecture/Model slider is hidden too, because the two cannot be combined.
- **Model vs textbook: planar CH4.** In this model the occupied nonbonding orbital of planar CH4 is a hydrogen-only combination, at the H 1s energy, and C 2pz is the empty LUMO. Textbooks make C 2pz the occupied lone pair, because H–H repulsion (absent here, B5) raises the hydrogen combination. The path's note says so.

## Release 2, Milestone E: doubt toggles, alternative schedules, compare and linked highlighting, design choices (delegated, 2026-11-07)

- **Alternative schedules (P5)** (`engine/lgo_engine/schedule.py`). Besides the lecture order, a case may offer:
  - *a different pairing:* another symmetry operation that pairs as many atoms. Pairings related by a symmetry operation give the same plateaus, so only one per class is kept. SF6 gets the cis pairing (2-6, 3-7, 4-5) beside the lecture's trans pairs; H3+'s three choices are all equivalent, so it gets none;
  - *central atom first*, then the ligand stages;
  - *every coupling at once*;
  - *atom by atom*, in numbering order.
  Options that give the same stages as an earlier one are dropped: with ligand–ligand couplings off (B5), CH4's "central first" is the lecture order. Every option partitions the same couplings, so the endpoint is the same Hamiltonian; the app shows the largest endpoint difference, about 1e−14 eV, and the tests require 1e−10.
- **Doubt toggles (P4)** (`engine/lgo_engine/doubt.py`). "Not convinced by" chips list the rules the case uses on the selected schedule. A doubted rule becomes a checkpoint: a step on the timeline that Watch stops at, carrying the engine's check of the rule and something to try. Doubts persist across cases, since they are about the rules, not the molecule.
  - *Where each rule sits:* R5, R6, R9 when the ligand set is finished (for a one-set molecule such as H3+ or allyl, at the end); R7 at the start of the central-atom stage; R2 at the end of the last stage that couples orbitals of different energy; R3 where the omitted couplings would have switched on; R1, R4, R8 on the finished diagram.
  - *When a rule counts as used:* R5 and R7 belong to the symmetry route, so "every coupling at once" and "atom by atom" do not use them. R2 needs a coupling across an energy gap, so H3+ and allyl do not use it. R3 needs an atom pair left uncoupled, R6 and R9 a one-orbital-per-atom set. R10 (sandwiches) is used by no case yet.
  - *What each check shows:*
    - R1: the same-colour and opposite-colour contact totals per level;
    - R2: for each coupled pair, the gap, coupling and their ratio, with the perturbative shift and admixture against the two-level exact values. Where |t| exceeds half the gap, the estimate is reported as failing (water, CH4, SF6, N2 s–p: every pair here). The full-calculation spread of each central orbital is also shown;
    - R3: the omitted pairs with their distances, and what switching them on at their distance-law values would do. Water's level order changes; CH4's and allyl's do not;
    - R4: the orbital count;
    - R5: the LGOs from projection;
    - R6: the counting rule against the exact energy for every combination. For a degenerate set it reports agreement if a ±1 member agrees, and keeps the reason for the member shown elsewhere. So H3+'s e′ is exact, as L4's done-when requires, while SF6's (2, 2, −1, −1, −1, −1) reason stays visible;
    - R7: the non-zero matching tallies and that all join equal irreps;
    - R8: counting and density-matrix bond orders side by side;
    - R9: node counts against the exact order.
  - *Alternatives:* paint phases, another schedule, the lecture/model slider, a morph (central sweep; closing the allyl ring), the symmetry-allowed lines, or the electron stepper, whichever exists for the case.
- **Speed (D10).** The schedule options depend only on geometry, symmetry, edges and which couplings are non-zero, so they are cached on those, and the pairing search keeps only the best involutions per set. A cold first view of SF6 costs about 30 ms in CPython, and the browser's worst warm recompute stayed at 47 ms.
- **Compare (D9).** A *Compare with* menu under the morph bar. It lists suggested partners first, then any case:
  - the allyl cation against the radical and the anion (the same diagram with 3 or 4 π electrons);
  - each period-2 diatomic against its period-3 partner (N2/P2, O2/S2, F2/Cl2, C2/Si2, B2/Al2), and CO against N2 (isoelectronic);
  - H2O model against lecture, H3+ against neutral H3, CH4 against its four hydrogens alone.
  The comparison is a second, finished diagram beside the first, on the same energy range. Levels with the same symmetry label (number and irrep, e.g. 1πu) highlight together in both, and a selected level's namesake stays highlighted. The spec's H2O vs H2S and CH4 vs SiH4 need H2S and SiH4, which wait for CONVENTIONS v1.5.
- **Linked highlighting (D7).** One mechanism for every pane. Each item carries identity tokens, plus tokens for its counterparts, and hovering or focusing an item lights everything that shares one. A line (a correlation line, a structure contact) lights only when one of its own endpoints is the item. The links:
  - an MO level ↔ the fragment orbitals it is built from and the mixing lines joining them, its photoelectron band, its Walsh line during a morph, and the matching-table rows and columns of its fragments;
  - a fragment orbital ↔ its MOs and its matching-table row or column;
  - a matching-table cell ↔ both its fragment orbitals;
  - an atom or lobe in the structure ↔ the ledger rows it takes part in;
  - a ledger row ↔ its contact line, its atoms, and the 3D view (as before).
  Electrons sit on their level and link with it. The lighting is re-applied after every redraw, so it survives a recompute under the pointer.
- **Engine output added for D7 and display:**
  - matching-table rows and columns name their fragment orbital (`frag`);
  - the ledger gives each basis function's atom (`basis_atoms`);
  - an atomic column whose group is one orbital on several atoms gets a display name ("C 2p" rather than "C1px/C2px/C3px"); the label is unchanged.
  No energy, label or check changes.
- **Fixed along the way:**
  - mid-construction MO rows with several levels drew stray commas between empty labels;
  - in the two-column layout the atomic column's labels sat in the gutter where the correlation lines leave the bars; they now go in the outer gutter.

## B6. Sign-off (owner only)

Not delegated. On 2026-11-06 the owner gave standing permission for Claude to apply the drafted changes. They are now in `CONVENTIONS.md` v1.3 and the spec, and the drafts in `drafts/` remain as the record. Applying them is not sign-off: v1.3 stays "not frozen" until the owner reviews it.

**Signed off, 2026-11-07.** After reading `drafts/signoff-summary-v1.4.md`, the owner wrote "I sign off v1.4". `CONVENTIONS.md` v1.4 is frozen. The signed text is commit 81cf760, sha256 `ca861a1d6f47f5fe9962eb229a3fbdf4313c22ed2e83083dcfa8fdbe71b98d36`. Any later change needs a new version and the owner's approval. Three points the summary flagged as resting on weaker sources remain open for correction in a later version:
- the Miessler water and CH4 LUMO orderings behind the hydrogen factor's upper limit;
- CH4's second photoelectron band near 23 eV;
- the inconsistent CCCBDB entry for Al2.
