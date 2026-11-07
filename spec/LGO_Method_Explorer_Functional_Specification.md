# LGO Method Explorer: Functional Specification

Nov 6, 2026 · @void · revised 2026-11-06 under owner permission (drafts/spec-edits-pending-signoff.md D1–D13)

## 1. Purpose and design principles

The app visualises the LGO technique taught in CM2021 (separate central and terminal atoms, fix one orbital red, flip its partner, combine equal-energy sets, put the central atom back, match colours, fill, count bond order) and shows how each of the instructor's rules follows from, or approximates, a simple orbital-interaction model. The user works in the instructor's vocabulary; the model runs underneath as a referee.

- **Method first, model underneath.** Every screen starts from a step of the lecture procedure. Matrices, irreps and parameters sit one click deeper.
- **Every rule is visible arithmetic.** Colour matching, contact counting and overlap checks are all shown as the same signed sum (section 2).
- **Approximations are switches, not footnotes.** Where the lecture simplifies, the app says so and lets you turn the simplification off and watch what changes.
- **Everything is a continuous family.** Static diagrams are snapshots of paths you can scrub: bending, ring closure, adding a central atom, changing electron count (section 5).
- **The model is also a model.** The app shows which rules follow from a simple coupling model and which add assumptions. It does not present the model as the final truth.

IDs: user-facing functions use L (ledger), M (matcher), P (method), D (dynamic). Engine functions keep their F numbers (sections 8 onward).

**Related tools.** The nearest existing project is CrystOD (2026, MIT licence), a command-line toolkit that finds a molecule's point group, computes SALCs and draws an extended-Hückel MO diagram as an HTML page with orbital sketches on hover. It does not show the stepwise LGO construction, map energies onto teaching rules, or run in a browser, which are this app's purpose. It is used here as an independent checker (section 12). Other tools are components rather than alternatives: symmetry back-ends (MolSym, libmsym, pgsym), Hückel π calculators (HMO), solid-state Slater–Koster libraries (chinook, pysktb, tightbinder), and diagram drawers (PyMoDia, mohelper).

## 2. Core idea: contact ledger and rule map

Every quantity the lecture reasons about is a sum of signed atom–atom contacts. For a combination with coefficients c\_i, its energy is

```latex
E = \sum_i c_i^2\,\varepsilon_i \;+\; 2\sum_{i<j} c_i c_j\, t_{ij}
```

where t\_ij is the coupling between orbitals i and j (negative for a stabilising in-phase contact). For N equivalent atoms with coefficients ±1/√N and equal neighbour coupling t, this reduces to the lecture's counting rule:

```latex
E = \varepsilon + \frac{2t}{N}\,(n_{\text{in-phase}} - n_{\text{out-of-phase}})
```

CH4 check: the all-red LGO has 6 in-phase contacts, giving ε + 3t; each two-red-two-blue LGO has 2 in-phase and 4 out-of-phase, giving ε − t. Those are exactly the model's eigenvalues.

**Rule map.** Each instructor rule, what it is in the model, and the badge it earns when used.

| Rule | In the model | Badge |
| --- | --- | --- |
| R1. Same colour → bonding | Sign of each contact's contribution | Exact, given a fixed phase convention |
| R2. Similar energy → interact | The energy shift scales roughly as coupling² / energy gap, and the admixed weight as (coupling / gap)² | Approximation: never strictly on/off |
| R3. Close enough → interact | Coupling magnitude falls with distance | Qualitative: it is a parameter |
| R4. n atomic orbitals → n MOs | Basis dimension | Exact |
| R5. Fix one, flip the partner, combine equal-energy sets | Stepwise symmetry adaptation through smaller subgroups | Exact; may differ from textbook LGOs by a rotation within a degenerate set |
| R6. Count in-phase vs out-of-phase contacts | Formula above | Exact for ±1 coefficients and equal contacts; approximate otherwise (unequal coefficients such as allyl's ½, 1/√2, ½ or octahedral 2, −1, −1; cis vs trans contacts) |
| R7. Put the central atom back at its original position | Couplings depend on direction through angular factors | Exact |
| R8. Bond order = (bonding − antibonding) / 2 | Counting definition | Convention, not an observable |
| R9. More nodes → higher energy | Sign changes in the coefficients. A node between atoms is an out-of-phase contact and costs energy; a node through an atom is a zero coefficient, so that atom's contacts drop out and cost nothing | Exact ordering for open chains with nearest-neighbour coupling (oscillation theorem for tridiagonal matrices); exact up to degenerate pairs on rings; ordering heuristic for 3D sets. Gives order, not gaps |
| R10. Sandwich matching: count nodal planes through the axis | A ring orbital with k nodal planes through the axis matches central orbitals with \|m\| = k; top-plus-bottom combinations match even central orbitals (s, d), top-minus-bottom match odd ones (p, f) | Exact for an ideal sandwich; \|m\| is only defined modulo the ring size, so a 6-ring's k = 3 orbital matches only one of the f (m = ±3) pair |

| ID | Function | Done when |
| --- | --- | --- |
| L1 | **Contact ledger.** For the current orbital or step, list every contact: the pair, in-phase or out-of-phase, weight (c\_i c\_j t\_ij), running total. | Ledger total matches the engine's energy or matrix element to 1e−10. |
| L2 | **Contacts on the structure.** Lines between atoms: solid for stabilising, dashed for destabilising, thickness ∝ size. Red and blue stay reserved for orbital phase. | Every ledger row has a matching line and vice versa. |
| L3 | **Rule badges.** Each step carries the R-number it uses and its badge; clicking opens the switch that removes that approximation (section 5). | Every step in every covered case has a badge. |
| L4 | **Rule check.** When R6 or R9 is applied, show its prediction beside the exact result and, if they differ, the reason. Where one rule is approximate and another exact, show both. | Exact agreement for CH4 and H3+; for allyl ψ1, R6 badged approximate and R9 exact, both shown; the stated reason appears for octahedral eg. |
| L5 | **Node markers.** Nodal planes drawn on the structure, with through-atom and between-atom nodes styled differently and linked to their ledger rows (zero-weight rows and out-of-phase rows respectively). Node count shown beside each level. | Allyl ψ2 shows one through-atom node on C2 linked to zero-weight rows; ψ3 shows two between-atom nodes linked to out-of-phase rows; every node count equals the sign changes in the engine's eigenvector. |

## 3. Molecule input and case matcher

The user types or builds any molecule. The app returns one of three verdicts (drawn, not built yet, or outside the method) and always shows why. The matcher is rule-based and deterministic; it never calls a language model.

**Case templates.** A case is defined by topology, symmetry, basis and electron count, not by formula.

| Template | Shape | Worked example | Other members |
| --- | --- | --- | --- |
| T0 | Two atoms | Lecture diatomics (H2, second-row) | B2, C2, N2, O2, F2, CO |
| T1 | Central atom + 2 terminals, bent | H2O | H2S, H2Se, NH2− |
| T2 | Central atom + 3 terminals, trigonal planar | BH3 | AlH3, BF3 (σ only, warned) |
| T3 | Central atom + 4 terminals, tetrahedral | CH4 | SiH4, NH4+, BH4−, CF4 (σ only, warned) |
| T4 | Central atom + 6 terminals, octahedral | SF6 | SeF6, SiF6 2−, PF6− |
| T5 | Ring of n equivalent s orbitals | H3+ | H3, H3− (D3h), square H4 once D4h ring is enabled |
| T6 | Open chain of n p\_π orbitals | Allyl cation | allyl radical and anion, butadiene (with conformation caveat) |
| T7 | Ring of n p\_π orbitals | — | cyclopropenyl cation, benzene once D6h data is added |
| T8 | Ligand set alone, no central atom | Six H in O\_h | four H in T\_d (the CH4 LGOs) |
| T9 | Central atom between two parallel rings (sandwich) | Beyond the course; reference case uranocene U(C8H8)2 | ferrocene, bis(benzene)chromium, an arene–U–arene sandwich |

Sandwiches (T9) use two tiers of LGOs: each ring's π orbitals (from T7), then top-plus-bottom and top-minus-bottom combinations (the lecture's fix-and-flip applied to whole rings), matched to central orbitals by R10. The central basis gains f orbitals for actinides; couplings are set per ring orbital and per |m| type rather than per atom pair, which avoids f-orbital Slater–Koster tables. Every T9 case runs with the fidelity bar (section 7) showing its heavy-element, spin-state and electron-counting chips amber or red.

| ID | Function | Done when |
| --- | --- | --- |
| M1 | **Parse input:** formula with charge ("C3H5+"), common name, or SMILES. | Each fixture input parses to the right atoms and charge. |
| M2 | **Get a structure.** Curated library first; otherwise VSEPR from the formula for central-atom molecules; π systems need connectivity from SMILES or a library entry. | Ambiguous formulas (C3H5+: allyl or cyclopropyl cation) prompt a choice rather than a guess. |
| M3 | **Classify** against the templates by central atom, equivalence of terminal atoms, point group and electron count. | The fixture list returns its expected verdicts. |
| M4 | **Three-way verdict**, always with the evidence (shape, equivalence, group, electron count): **drawn** (worked case, template member, or flagged case below); **not built yet** (the method applies; names what is missing); **outside the method** (names why). | Every fixture lands in the expected class with its evidence listed. |
| M5 | **Nearest covered case.** For anything not drawn, name the closest drawn case and what differs, with a morph path if one exists. | NH3 points to BH3 with the pyramidalisation morph. |
| M6 | **Template builder.** Pick a template, swap atoms (e.g. O → S, H → F), adjust bond lengths and angles within the template's symmetry. The result goes straight to the matcher. | Every builder output is drawn; the builder cannot produce a dead end. |
| M7 | **Free editor (later).** An embedded molecule editor (JSME or Ketcher) outputs SMILES to the matcher, so anything can be entered and sorted into the three verdicts. | Editor output and typed SMILES give identical verdicts. |

**Reason codes**, by verdict class:

*Not built yet* (the method applies; data or a template is missing):

- **Point group not in library**, e.g. NH3 (C3v). Fix: a group data file.
- **π from terminal atoms**, e.g. CO2, NO3−, the π part of BF3. A σ-only drawing is offered meanwhile.
- **Multicentre template**, e.g. B2H6.
- **Groups with complex characters**, e.g. C3, S4, T (section 14).
- **Element without parameters**: needs an entry in the parameter scheme (F2.12).
- **Too large for step-by-step** (more than about 20 basis functions): drawn as a result only, without method mode.

*Outside the method* (LGO construction does not meaningfully apply):

- **Extended solids**: SiO2 networks, salt lattices, metals. These need band theory, not a finite molecule.
- **No determinable structure**: the user is asked to pick an isomer or give SMILES; if none exists, the input is outside.
- **Not an orbital question**: bond-enthalpy comparisons, oxidation states, inert-pair trends asked as such (tutorial E–G).

**Flagged but drawn.** Two kinds of case are drawn with warnings rather than turned away:

- **Inequivalent ligands** (CH3Cl, SF4, PCl5). The engine solves them like anything else; symmetry labels come from the lower point group, or the C1 fallback (F2.14) if that group is not in the library. What changes is the method: fix-and-flip pairing, "combine equal-energy sets" and the counting rule R6 apply only within each set of equivalent atoms, and bond orders are reported per bond from the density matrix because equal sharing (R8) no longer holds. The method-fit indicator (section 7) marks the steps that go beyond what has been taught.
- **A degenerate level partly filled** (neutral cyclic H3), shown with a Jahn–Teller warning.

Automatic point-group detection for arbitrary geometries needs tolerance-based detection; MolSym or libmsym are candidates here.

## 4. Method mode: the instructor's procedure as steps

The lecture procedure is a schedule for switching on couplings in stages: contacts within each pair (fix and flip), then contacts between pairs (combine equal-energy sets, in one or more rounds), then the central-atom coupling (put the central atom back). The app plays it as one continuous construction slider, not eight screens. It runs in the lecture's vocabulary (red and blue lobes, LGO1′ numbering, bonding / nonbonding / antibonding), with irrep labels as an optional overlay.

```latex
H(s) = H_0 + \sum_{k=1}^{K} \lambda_k(s)\, H_k, \qquad s \in [0, 1]
```

H\_k is the set of couplings switched on at stage k, as generated by F2.13. The generator proposes, for CH4, K = 3 (within pairs, between pairs, central atom) and, for SF6 as lectured, K = 4 (within each axis pair, combining two axes, adding the third, central atom). A stage whose newly enabled couplings are identically zero is skipped and the slider shows "No ligand–ligand coupling in this model" (DECISIONS.md B5): CH4, H2O and BH3 have H–H couplings off, so CH4 runs with one effective stage, while SF6 keeps all four because its F–F couplings are on. Each λ\_k ramps from 0 to 1 in turn as s advances. At s = 0 all atoms are separate; at s = 1 the diagram is final. The endpoint does not depend on the stages or their order, which is the rigorous form of the lecture's "it's up to you".

| ID | Function | Done when |
| --- | --- | --- |
| P1 | **Construction slider.** Starts at s = 1 (finished diagram); scrubbing replays the construction with structure, ledger, node markers and diagram moving together. | Endpoint equals the engine result; intermediate plateaus reproduce the lecture's LGO′ energies (CH4: LGO1′ level with LGO3′, LGO2′ with LGO4′). |
| P2 | **Checkpoints.** The lecture steps are labelled ticks on the slider (separate, fix and flip, combine, place the central atom, match, fill and count); snapping to one shows its narration. | Every tick lands on a plateau of the relevant ramp. |
| P3 | **Why-jump.** Tap any level, degeneracy, label or line and the slider jumps to the moment it appeared, with that moment's ledger and node markers. | Every displayed feature has a recorded origin on the slider. |
| P4 | **Doubt toggles.** Each rule (R1–R10) is marked convinced or not. Convinced rules play as animation; a rule marked not convinced expands its segment into a discrete step with the ledger and alternatives to try (doubting R2 replaces the on/off gap rule with actual mixing). | Expanded checkpoints equal the distrusted rules the case uses. |
| P5 | **Alternative schedules.** Pick a different pairing or switching order; the path differs, the endpoint does not. | Endpoint energies identical to 1e−10 for every schedule. |
| P6 | **Matching view.** As the central coupling ramps, per-lobe overlap tallies show which central orbitals pair with which LGOs; zero tally means no interaction (e.g. O 2p\_x with water's LGO1). | Matches agree with the engine's symmetry blocks. |
| P7 | **Fill and count.** Electrons from both sides (neutral-atom counting as in lecture, or a stated alternative); HOMO, LUMO, bond order per bond by the counting definition (R8), density-matrix value as an advanced toggle. | Water 1 per O–H, SF6 ⅔ per S–F, allyl cation 0.5 per C–C (π), H3+ ⅓ per H–H, BH3 1 per B–H. |

Three ways to run it:

- **Watch:** the app performs every step with narration.
- **Practice:** the user performs each step; a wrong colour or pairing gets an explanation from the ledger, not just a red cross.
- **Exam:** a blank diagram to fill in, checked at the end.

## 5. Making it dynamic

The cases are few and fixed, so the interaction carries the insight: each worked case is one point on paths the user can scrub, paint and compare. All transitions animate. Levels are tracked by symmetry label wherever the path keeps a symmetry, so crossings and avoided crossings read correctly; where it does not (free dragging, a C1 stretch), each level is followed by maximum overlap between successive eigenvectors, so no level swaps identity mid-animation.

**Morph paths between cases (D1).** Scrubbable sliders that turn one case into another:

| Path | What you see |
| --- | --- |
| Allyl chain → cyclic ring | Closing the chain adds one contact: energies go from ε ± √2t, ε to ε + 2t, ε − t, ε − t. Tutorial A and B become one problem; node markers show the one-node and two-node chain orbitals becoming the ring's degenerate pair. |
| H2 → linear H3 → cyclic H3+ | The third atom's contacts switch on one by one; the e′ pair becomes degenerate only at the triangle. |
| H2O bending 90° → 180° | Walsh diagram: why water is bent in MO terms, as the lecture asks. |
| BH3 → NH3 pyramidalisation | The empty p\_z (Lewis-acid orbital) mixes in as the molecule bends (needs C3v data). |
| CH4 tetrahedral → square planar | A bonding orbital turns into a nonbonding C 2p\_z as the geometry flattens. |
| Six H in O\_h → SF6 | Insert the central atom and raise its coupling from zero: LGOs split into bonding/antibonding pairs; the two with no partner stay put. |
| Central-atom energy sweep | Move a central orbital's energy through the LGO levels: mixing peaks near resonance and character swaps across the avoided crossing. R2 as a continuum. |

| ID | Function | Done when |
| --- | --- | --- |
| D1 | **Morph paths**, as above, each with the diagram, structure and ledger updating live. | Endpoints reproduce their worked cases exactly. |
| D2 | **Approximation switches.** Each badge-able simplification (for example water's O 2s as nonbonding) becomes a slider from "lecture" to "model". | At "lecture" the diagram matches the slide; at "model" it matches the full engine. |
| D3 | **Phase painting.** Click lobes to set colours yourself; the app gives that combination's energy and node count, then shows whether it is a true LGO or a mixture of LGOs, with the mixture's weights. | Painting a valid LGO shows weight 1 on one level; any other pattern shows its decomposition. |
| D4 | **Degenerate-set dial.** Rotate continuously between equivalent choices of degenerate LGOs (e.g. the lecture's top/bottom split for CH4 and another split). | Energies and degenerate-summed weights stay fixed throughout. |
| D5 | **Lone-pair dial.** Mix water's two nonbonding orbitals from the canonical pair to two equivalent VSEPR-style lone pairs. | Total electron density and energy shown unchanged at every angle. |
| D6 | **Electron stepper.** Add or remove electrons (allyl cation → radical → anion; H3+ → H3 → H3−) with HOMO, LUMO and bond orders updating. | Jahn–Teller warning appears at neutral H3. |
| D7 | **Linked highlighting.** Hover any lobe, contact, LGO, level or electron to highlight its counterparts in every pane. | Each item has counterparts in all panes where it exists. |
| D8 | **Case map.** A graph of all cases joined by the morph paths; the current molecule glows, and the matcher (section 3) places the user's input on it, or just outside with its reason code. | Every D1 path is an edge; every T-template has a node. |
| D9 | **Compare.** Two molecules side by side with aligned diagrams (H2O vs H2S, CH4 vs SiH4, allyl cation vs anion). | Both diagrams share the energy axis and highlighting. |
| D10 | **Responsiveness.** Recompute and redraw in under 100 ms for bases up to about 20 functions on a mid-range phone, measured warm. First load (Pyodide runtime plus the NumPy download) has its own budget, set when the reference device is chosen. | Measured on the reference device, warm and cold separately. |

**How the dynamic features are derived.** Every D-function above came from the same recipe. Run it on any new case to get its dynamic features:

1. **Fix the worked cases as anchor points.** The lecture and tutorial diagrams are the endpoints every path must reproduce exactly.
2. **List each case's hidden parameters.** Geometry (bond angle, ring closure, pyramidalisation), couplings (central atom on or off, ligand–ligand), on-site energies (the gap), electron count, and phase choices.
3. **Connect cases through those parameters.** Where two cases share a template, or one is a limit of the other (chain closing into a ring, a ligand set gaining a central atom, water straightening to linear), make the path a scrubbable morph and an edge on the case map.
4. **Expose what each rule silently fixes.** For each R-rule, find the parameter it holds constant or switches off (R2 the gap, R3 the distance, each lecture simplification a coupling set to zero) and turn it into a slider from "lecture" to "model".
5. **Turn every arbitrary choice into a dial.** Wherever the lecture says "it's up to you" (pairing, which degenerate LGOs, canonical vs localised lone pairs), make the choice continuous and display what stays invariant.
6. **Let the user be wrong.** Wrong colours, a misplaced central atom, combining unequal energies: the engine decomposes the result and explains it instead of rejecting it.
7. **Keep every path honest.** Track levels by symmetry (or by eigenvector overlap where the path has none) so crossings and avoided crossings read correctly; show node counts and ledger totals throughout; animate the interpolation.
8. **Link all views.** One change propagates to structure, ledger, node markers, diagram and method step at once.
9. **Test each path.** Endpoints match the worked cases; invariants stay fixed along every dial; recompute stays within D10's budget.

## 6. Coverage: CM2021 lectures and Week 7 tutorial

**Release plan.** Each release ships only when its tests pass; later releases add functions, they never rework earlier ones.

| Release | Contents | Gate |
| --- | --- | --- |
| Release 1 · Core | Engine (sections 8–9, with F2.12 limited to the Release 1 elements and F2.13 for slider schedules); CONVENTIONS.md; ledger L1–L3; slider P1–P3 and P6–P7; Watch mode; diagram F3.1–F3.4, F3.7, F3.9; export (section 11); the seven cases below plus the lecture diatomics | CONVENTIONS.md frozen with its open items closed; all section 12 oracles and property tests pass (`make check`); the CrystOD cross-check is deferred until the programme is finished (owner decision B1) and runs as `make check-release`; every case below runs end to end with its lecture simplifications badged |
| Release 2 · Dynamic | L4–L5; P4–P5; Practice and Exam modes; D1–D7, D9–D10; F3.6, F3.8 | Morph endpoints reproduce their cases; invariants hold along every dial; lecture simplifications switchable (D2) |
| Release 3 · Open input | Matcher M1–M6; molecule packs; full parameter scheme F2.12; F2.14–F2.15; fidelity bar B1–B5; case map D8 | Fixture verdicts correct; every green or red chip cited |
| Release 4 · Beyond the course | Sandwiches (T9); free editor M7; advanced complex mode; further point groups such as C3v | Per-template oracles pass |

Release 1 cases:

| Case | Source | The app must show | Lecture simplifications (badged; switchable from Release 2) |
| --- | --- | --- | --- |
| Diatomics (H2, second-row) | Lecture | σ/σ\* from head-to-head and π/π\* from side-on overlap by the colour rule; filling; bond order | s–p mixing is included; under the current parameter scheme B2, C2 and N2 still order σg(2p) below πu, and B2 and C2 get the spin state wrong. Shown as red fidelity chips, not tuned away (DECISIONS.md M1) |
| H2O | Lecture | Two LGOs from the two H 1s; O 2p matching; 2 bonding, 2 nonbonding, 2 antibonding; bond order 1 per O–H | O 2s treated as nonbonding (energy gap) |
| CH4 | Lecture | Top/bottom pairing on the cube; one low LGO and three equal; 4 bonding MOs; photoelectron spectrum shows two bands with degeneracies 3:1 (badged approximate: band intensities follow photoionisation cross-sections, not occupancy, and the 2a1 band near 23 eV lies beyond He I) | none |
| SF6 | Lecture | Six inward F 2p orbitals combined by axis; S 3s and 3p only; two LGOs left nonbonding; bond order ⅔; no d-orbital contribution | none: the lecture preset equals the model. Cis F–F contacts (2.21 Å) and trans contacts (3.12 Å) both act, and the cis contacts are what separate eg from a1g |
| Allyl cation | Tutorial A | Three π levels ordered by node count 0, 1, 2 (R9, exact for chains); HOMO bonding, LUMO nonbonding with a node on the middle carbon; π bond order 0.5 per C–C | σ framework assumed single bonds |
| H3+ | Tutorial B | One low LGO-like level below an equal pair; 2 electrons; bond order ⅓ per H–H | none |
| BH3 | Tutorial C | Bonding a1′ and e′ sets filled; empty nonbonding B 2p\_z as LUMO, flagged as the Lewis-acid orbital; bond order 1 per B–H | none |
| Six H in O\_h | Tutorial D | Six LGOs with phases; order: all-in-phase lowest, then three, then two | Needs H–H coupling switched on; without it all six are degenerate |

Out of scope from the tutorial: question E (C vs Si oxides by bond energies), F (CCl4 vs PbCl4, inert-pair effect) and G (nitrogen oxide structures and oxidation states). None is an MO-diagram question.

## 7. Fidelity bar

The course rules draw a diagram for any geometry they are given; they never decide whether a molecule exists or what shape it takes. Where they are used to predict (bond order as stability, Walsh-style comparison of geometries), they fail in recurring ways. A permanent bar shows one chip per failure category for the loaded molecule.

- **Green:** checked against curated data, and the course rules get it right. A check whose outcome a parameter was chosen to reproduce (for example the orderings that fix the 0.6 hydrogen scale) cannot be green.
- **Red:** checked against curated data, and they get it wrong.
- **Amber:** a detector fired, but nothing has been checked.
- **Grey:** not applicable.

Green means the course rules got this molecule's checked property right, not that the model is right in general.

| Category | What goes wrong | Detector (amber) | Curated checks |
| --- | --- | --- | --- |
| Existence | A diagram is drawn for anything; bond order above zero is not existence, since there is no energy scale and no comparison with dissociation products | Terminal atoms with no unpaired electron (noble gases); low bond order per bond | He2 green (bond order 0, unbound); H3+ green (exists); arene–U–arene amber |
| Geometry | Structure is an input; Walsh comparisons can mislead | Degenerate level partly filled (Jahn–Teller); heavy alkaline-earth dihalides | H2O bent green; neutral H3 red as a stable triangle; BaF2 red (predicted linear, observed bent; mechanism debated) |
| Orbital ordering | The on/off energy-gap rule misses s–p mixing | Small s–p gap relative to couplings (B, C, N) | B2, C2 and N2 red: even with s–p mixing on, the current parameter scheme puts σg(2p) below πu (DECISIONS.md M1); O2 and F2 green |
| Spin state | Bottom-up filling fails for nearly degenerate levels and open d or f shells | Small frontier gap relative to coupling; open d or f shell | O2 triplet green (with Hund's rule); B2 red (model closed shell, reference ³Σg); C2 red (model partly fills πu, reference ¹Σg⁺); CH2 red (triplet ground state) |
| Heavy elements | No relativistic effects; course basis lacks d and f | Z above about 50; any f-block atom | PbCl4 vs CCl4 red (tutorial F); uranium sandwiches red |
| Aggregation | Single-molecule picture; cannot compare monomer and condensed forms | π bonds between atoms below period 2; Lewis-acid and Lewis-base sites on one molecule | S2 vs S8 green once R3 ties π coupling to distance; SO3 polymerisation amber |
| Bond-order meaning | A convention; equal sharing assumes equivalent bonds | Inequivalent bonds | Informational only (e.g. SF6's ⅔ is a counting result) |
| Electron counting | Neutral-atom counting fails for closed-shell donors and charged ring ligands | Closed-shell terminal atoms; charged rings (C8H8 2−) | Uranocene switches to ionic counting, flagged |

| ID | Function | Done when |
| --- | --- | --- |
| B1 | **Bar always visible**, one chip per category, updating as the molecule or slider changes. | Every case in section 6 shows all eight chips. |
| B2 | **Tap a chip** to see the reason and a link to the slider, switch or morph that demonstrates it (e.g. the R2 switch shows how far s–p mixing moves B2's levels, and that under the current parameters it is not far enough). | Every amber or red chip links to a demonstration or states that none exists. |
| B3 | **Curated checks** stored as data, each with its source reference. No green or red without a citation. | Every green and red entry resolves to a cited source. |
| B4 | **Detectors** are documented deterministic rules with stated thresholds. | Each detector has fixture molecules that fire it and that don't. |

**Method fit, separate from accuracy.** Inequivalent ligands do not make predictions less accurate; they take the case beyond the taught method. Mixing the two would make the traffic lights mean two things, so method fit gets its own indicator beside the bar: blue for "taught method applies as taught", purple for "some steps go beyond what has been taught", listing those steps (e.g. pairing only within equivalent sets, per-bond bond orders).

| ID | Function | Done when |
| --- | --- | --- |
| B5 | **Method-fit indicator** beside the bar, listing the steps that go beyond the taught method for this molecule. | CH3Cl and SF4 show purple with their affected steps; every section 6 case shows blue. |

---

*Sections 8–14 specify the engine that sections 1–6 run on. It is unchanged in substance from the earlier spec: it computes, checks and labels; it never drives the interface. Tech stack, pinned in CONVENTIONS.md: CPython 3.12.7 + NumPy 2.0.2 engine; Pyodide 0.27.7 in a Web Worker; vanilla ES2022 + HTML/CSS/SVG + WebGL2 frontend; static hosting; no backend.*

## 8. Engine: system definition

A course case, the matcher, the builder or a molecule pack supplies the geometry, basis, parameters and electron count; the engine validates them into one configuration object that everything downstream reads.

| ID | Function | Done when |
| --- | --- | --- |
| F1.1 | **Geometry presets.** Release 1: H2O (C2v; bond length and HOH angle as inputs), CH4 (T\_d), SF6 (O\_h), allyl (C2v, π only), cyclic H3 (D3h), BH3 (D3h), six H in O\_h with no central atom, and the lecture diatomics. Later releases: template members via molecule packs; complexes (O\_h ML6, T\_d ML4, D4h ML4) in advanced mode; sandwiches. | Atom positions match the table in CONVENTIONS.md exactly. |
| F1.2 | **Distortion paths** (Release 2), as used by the D1 morphs: allyl ring closure, H2 → linear H3 → cyclic H3+, H2O bending 90° → 180°, CH4 flattening (T\_d → D2d → D4h), six H gaining a central atom, isosceles H3 (D3h → C2v); BH3 pyramidalisation once C3v is added. Complex distortions belong to advanced mode. | Detected point group drops as expected and degeneracies split as the correlation table predicts. |
| F1.3 | **Per-atom basis.** Central atom: s and p, with d off by default for every main-group atom; f orbitals only for actinides in sandwiches (Release 4). Peripheral atoms: σ donor only; σ donor + two π functions perpendicular to the bond; full s and p; or a single p\_π for π-only systems. | Basis size and output MO count update together. |
| F1.4 | **Coupling model.** Atom-pair Slater–Koster over a neighbour list (default); AOM, central–peripheral pairs only (advanced complex mode); ring-orbital couplings per \|m\| type (sandwiches, Release 4). | On an ML6 complex with peripheral–peripheral terms off, Slater–Koster and AOM give identical results. |
| F1.5 | **Parameters.** On-site energies per element and orbital type; Slater–Koster pair parameters per element pair and bond class, from the sourced parameter table (F2.12). | Every parameter appears in the exported config with name, units and source. |
| F1.6 | **Presets.** For each course case, a *lecture* preset reproducing the slide (e.g. water's O 2s decoupled; for SF6 the lecture preset equals the model) and a *model* preset from the parameter table; complex presets (σ-only, π-donor, π-acceptor) in advanced mode. Every preset is tagged with which it is. | Loading a preset shows its tag on screen and in exports; the lecture preset reproduces the slide's diagram. |
| F1.7 | **Electron count.** Entered directly, or derived from a named convention: neutral-atom counting (the lecture's) for molecules, oxidation-state counting for complexes. The convention is stored with the result. | SF6 σ-only gives 12 electrons under both conventions; allyl cation, radical, anion give 2, 3, 4. |
| F1.8 | **Config validation.** Schema check before anything runs. | Unknown keys, missing parameters and empty bases are rejected with a readable message. |

In the advanced complex mode, the familiar AOM parameters emerge in the weak-coupling limit, e.g. e\_σ ≈ (sdσ)² / (ε\_d − ε\_Lσ), treating the ligand σ donor as an s-type function. The site should say so wherever both appear.

## 9. Engine: Hamiltonian and analysis

The engine is a pure, deterministic Python package (numpy only, no I/O) that turns a configuration into one serialisable Result. The same package runs in the notebook and in the browser via Pyodide.

The Hamiltonian in the atomic basis:

```latex
H_{a\mu,\,b\nu} = \varepsilon_{a\mu}\,\delta_{ab}\delta_{\mu\nu} + (1-\delta_{ab})\sum_{k\in\{\sigma,\pi,\delta\}} V_{k}(a,b)\,E^{k}_{\mu\nu}(\hat{\mathbf r}_{ab})
```

Here a, b are atoms, μ, ν their orbitals, V\_k the Slater–Koster pair parameter of bond type k (zero for pairs outside the neighbour list), and E the standard direction-cosine factors of the unit vector (l, m, n) from a to b. Examples: s–p\_x gives l·spσ; p\_x–p\_x gives l²·ppσ + (1 − l²)·ppπ; s–d\_z² gives (n² − (l² + m²)/2)·sdσ. With the ligand σ donor treated as an s-type function, the last reproduces the AOM factor F(d\_z²) = (3n² − 1)/2, so the complex workflow is unchanged. Ligand π functions are p orbitals perpendicular to the bond, so their frames follow from geometry.

| ID | Function | Done when |
| --- | --- | --- |
| F2.1 | **Build basis.** Atoms, per-atom orbitals (F1.3), ordered labels. | Labels and ordering match CONVENTIONS.md. |
| F2.2 | **Assemble H** from on-sites and Slater–Koster couplings over the neighbour list. | H symmetric to 1e−12; angular factors match an established Slater–Koster implementation (chinook or pysktb, test-only dependency); with peripheral–peripheral terms off, σ d-factors squared sum to 1 per ligand (AOM check). |
| F2.3 | **Diagonalise** with a symmetric eigensolver; sort ascending; group degenerate levels within a stated tolerance. | Residual ‖Hc − Ec‖ < 1e−10; eigenvectors orthonormal. |
| F2.4 | **Generic symmetry analysis.** Point groups stored as data (generator matrices plus character table), generated offline with an established library (MolSym, libmsym or pgsym; tool and version recorded). Characters of crystallographic groups are cross-checked against the Bilbao Crystallographic Server; D8h, which Bilbao does not cover, is checked by closure and orthogonality only. The basis representation is built automatically from atom permutations and orbital transformations (s: 1, p: the 3×3 rotation, d: its real 5×5 representation; f: its real 7×7 representation, for sandwiches in Release 4). Reduction and SALCs by projection. v1 groups: O\_h, T\_d, D4h, D2d, D3h, C2v, C2h (butadiene), D6h (benzene, arene sandwiches), D8h (uranocene), plus D2h and C2v as working subgroups for linear molecules (the diatomics, linear H2O), with π pairs labelled by their exact decomposition and mapped to σg, πu and so on, and C1 as the fallback (F2.14). | Irrep counts match reference tables; each group file passes closure and character-orthogonality checks. |
| F2.5 | **Label MOs** by projecting each degenerate eigenspace onto each irrep. The C2v axis convention (molecule in the yz plane, so out-of-plane p is b1) is fixed in CONVENTIONS.md, because b1 and b2 swap with it. | Every level gets one irrep with projection weight > 0.999 in ideal geometry. |
| F2.6 | **Composition.** Weight of each atomic orbital and each SALC block per level, summed over degenerate partners. | Weights sum to 1 per MO; summed values unchanged under rotation within the degenerate space. |
| F2.7 | **Bonding character** per level and per atom pair, from the sign of 2Σ c\_aμ c\_bν H\_aμ,bν: negative = bonding, positive = antibonding, near zero = nonbonding. | O\_h σ: t2g nonbonding. Allyl: middle MO nonbonding with zero coefficient on the central carbon. |
| F2.8 | **Filling** from the F1.7 count; aufbau or a user-specified number of unpaired electrons. Fixed warning that orbital energies alone do not settle spin state; extra warning when a degenerate level is partly filled (e.g. neutral cyclic H3, e′¹), since the ideal geometry is then Jahn–Teller unstable. | Total electrons conserved; both warnings fire in their test cases. |
| F2.9 | **Result object.** Config, energies, eigenvectors, labels, compositions, filling, version and a conventions hash, all as JSON. | Round-trip JSON → Result → JSON is lossless. |

Two functions added for release 1:

| ID | Function | Done when |
| --- | --- | --- |
| F2.10 | **Frontier orbitals.** Mark HOMO and LUMO; flag when either is degenerate or nonbonding. | Allyl cation: HOMO 1b1, LUMO 1a2. BH3: LUMO 1a2″, flagged nonbonding. |
| F2.11 | **Bond orders**, per symmetry-distinct bond, under two named definitions shown side by side. (i) Counting: ½(bonding − antibonding electrons), shared equally over symmetry-equivalent bonds; nonbonding electrons count zero. (ii) Density-matrix: Coulson p\_ab = Σ n\_i c\_ia c\_ib for one orbital per atom, Wiberg index for several. For π-only systems, a separate line adds an assumed σ single bond, labelled as an assumption. | Release 1 oracle values reproduced; every displayed bond order names its definition. |

Four functions added so scope can grow without code changes:

| ID | Function | Done when |
| --- | --- | --- |
| F2.12 | **Parameter scheme.** Default couplings from Harrison's universal law V = η·K/d² (K = 7.62 eV Å²; η: ssσ −1.40, spσ +1.84, ppσ +3.24, ppπ −0.81) and on-site energies from one teaching table (Miessler, Fischer & Tarr, via LibreTexts; DECISIONS.md M7), as frozen in CONVENTIONS.md, with per-molecule curated overrides. The η values are Harrison's 1980 book set (DECISIONS.md M2). Couplings involving H extrapolate Harrison beyond his table and are scaled by 0.6, a pedagogical choice inside the admissible range 0.4235–1.0313 for the required orderings (DECISIONS.md B2). Every parameter records its source. | Any element in the scheme's table runs without hand-tuning; overrides show their source in exports. |
| F2.13 | **Pairing-schedule generator.** For each set of equivalent atoms, pair atoms related by a symmetry operation, then pair the pairs, recursively; output the construction-slider schedule (P1) and lecture-vocabulary narration. | Reproduces the lecture's CH4 and SF6 constructions; generates a valid schedule for every template without hand-authoring. |
| F2.14 | **C1 fallback.** If the point group is missing from the library, solve and draw anyway with labels shown as "—"; method mode runs on the equivalent sets it can detect. | Any molecule with a structure and parameters can be drawn. |
| F2.15 | **Matcher–engine consistency test.** Every matcher fixture and every template member is run through the engine in CI. | A template member whose group or parameters are missing fails CI, never silently in the app. |

## 10. Diagram rendering

The engine emits a DiagramSpec (JSON: levels, labels, lines, weights); renderers only draw it and never compute. This keeps every plotted value traceable to the Result.

| ID | Function | Done when |
| --- | --- | --- |
| F3.1 | **Three-column diagram.** Central-atom orbitals left, MOs centre, ligand LGOs right, labelled in lecture style (LGO1′ numbering; bonding / nonbonding / antibonding) with irreps as an optional overlay (systems with no central atom, such as cyclic H3 and allyl, use two columns: atomic levels and MOs); degenerate levels drawn as side-by-side lines with a multiplicity tag (e.g. eg ×2). | Every level in the Result appears once, at its stored energy. |
| F3.2 | **Two line layers.** Dashed thin lines = symmetry-allowed interaction; solid lines with thickness ∝ composition weight = actual mixing, hidden below 5%. | Each layer toggles independently; no line encodes both meanings. |
| F3.3 | **Energy axis** labelled "model units" unless a calibration is loaded. | Label present in the render and every export. |
| F3.4 | **Electron filling** as arrows on levels, following F2.8. | Arrow count equals the electron count. |
| F3.5 | **Frontier bracket** marking the HOMO–LUMO gap and any nonbonding set (in advanced complex mode, the d-manifold splitting such as Δ\_o). | Bracketed value equals the eigenvalue difference in the Result. |
| F3.6 | **Symmetry-only view.** Collapses energies to an ordering of irreps and connectivity, showing only what symmetry establishes. | Switching views changes no labels or connections. |
| F3.7 | **3D viewer** (WebGL2, per the pinned stack) of the complex: ligand orbitals as lobes; a selected SALC or MO drawn with lobe size ∝ coefficient and colour = sign. | Selecting an item in the diagram highlights the same item in 3D. |
| F3.8 | **Correlation (Walsh-type) plot.** Level energies against a distortion coordinate or one parameter, labelled by irreps of the lower-symmetry group. | Levels are tracked by irrep, not energy rank, falling back to eigenvector overlap where the path has no symmetry: different irreps may cross, same irreps avoid crossing. |
| F3.9 | **Responsive layout.** Usable at 360 px width; diagram above inspector on phones. | No horizontal scrolling at 360 px. |

## 11. Explanation, provenance and export

Every displayed quantity carries a source tag, and everything shown can be exported and reloaded to reproduce the same diagram.

| ID | Function | Done when |
| --- | --- | --- |
| F5.1 | **Provenance on hover.** Each level, label and line is tagged geometry, symmetry, parameter (named) or convention. | No displayed quantity lacks a tag. |
| F5.2 | **Generated explanation.** A deterministic text account of how this diagram was built, filled from templates with Result values. No language model at runtime. | Same config always yields identical text. |
| F5.3 | **Figure export** as SVG and PNG; PDF via SVG. | Exported figure matches the on-screen DiagramSpec. |
| F5.4 | **Data export** as JSON: config, Result, software version, conventions hash. | Re-importing reproduces the diagram exactly. |
| F5.5 | **Config import** from JSON or file upload, validated by F1.8. | Invalid files are rejected with the failing field named. |
| F5.6 | **Shareable link.** Compressed config in the URL fragment; no server storage. | Opening the link reproduces the diagram on another device. |
| F5.7 | **Scope banner.** Persistent short note: model units, which preset (lecture or model) is loaded, not a spectral prediction. | Visible on every page and in exports. |

## 12. Verification hooks for the checker AI

The checker must derive expected values blind from this spec and CONVENTIONS.md before reading the code. A checker that reads the code first tends to confirm it, and two models can share the same convention error.

**Analytic oracles.** Each block is a small matrix whose eigenvalues are known in closed form. The values below are my own derivation from the angular factors in section 9; the checker (or the owner) re-derives them independently before they are used. Rows with d orbitals gate Release 4 (advanced complex mode); the O\_h a1g and t1u rows also define SF6 and so gate Release 1.

| Case | Block (per partner) | Matrix in basis (central…, SALC) |
| --- | --- | --- |
| O\_h σ | a1g | \[\[ε\_s, √6β\_s\], \[√6β\_s, ε\_L\]\] |
| O\_h σ | t1u | \[\[ε\_p, √2β\_p\], \[√2β\_p, ε\_L\]\] |
| O\_h σ | eg | \[\[ε\_d, √3β\_d\], \[√3β\_d, ε\_L\]\] |
| O\_h σ | t2g | ε\_d, nonbonding |
| SF6 σ, d off | a1g, t1u | as the O\_h rows with S 3s, 3p |
| SF6 σ, d off | eg | ε\_F, nonbonding, fluorine only (F–F coupling off in this generic row; with F–F couplings on, as in both SF6 presets, it sits at ε\_F − 2t\_c + t\_t, CONVENTIONS A.3) |
| T\_d σ | a1 | \[\[ε\_s, 2β\_s\], \[2β\_s, ε\_L\]\] |
| T\_d σ | t2 | \[\[ε\_p, 0, 2β\_p/√3\], \[0, ε\_d, 2β\_d/√3\], \[2β\_p/√3, 2β\_d/√3, ε\_L\]\] |
| T\_d σ | e | ε\_d, nonbonding |
| D4h σ | a1g | \[\[ε\_s, 0, 2β\_s\], \[0, ε\_d(z²), −β\_d\], \[2β\_s, −β\_d, ε\_L\]\] |
| D4h σ | b1g | \[\[ε\_d, √3β\_d\], \[√3β\_d, ε\_L\]\] |
| D4h σ | eu | \[\[ε\_p, √2β\_p\], \[√2β\_p, ε\_L\]\] |
| D4h σ | b2g, eg, a2u | central atom only, nonbonding |
| Cyclic H3, D3h (t = ssσ, nearest neighbours) | a1′ | ε\_H + 2t |
| Cyclic H3, D3h | e′ | ε\_H − t (×2) |
| H2O, C2v, HOH angle θ, H–H off | a1 | \[\[ε\_2s, 0, √2·ssσ\], \[0, ε\_2p, √2 cos(θ/2)·spσ\], \[√2·ssσ, √2 cos(θ/2)·spσ, ε\_H\]\] |
| H2O | b2 | \[\[ε\_2p, √2 sin(θ/2)·spσ\], \[√2 sin(θ/2)·spσ, ε\_H\]\] |
| H2O | b1 | ε\_2p, nonbonding (out-of-plane O 2p) |
| Allyl π (t = ppπ) | b1 | \[\[ε\_C, √2t\], \[√2t, ε\_C\]\] in basis (C2, (C1 + C3)/√2), so ε\_C ± √2t |
| Allyl π | a2 | ε\_C, nonbonding, node on C2 |

Here β\_x is shorthand for the relevant σ Slater–Koster parameter (sdσ, spσ or ssσ). Two limit checks on H2O: at θ = 180° the a1 coupling of O 2p\_z vanishes and the b2 coupling reaches √2·spσ, matching linear D∞h; the summed weights must vary smoothly along the bending path.

**Release 1 oracles** (my derivation, to be re-derived blind like the rest):

| Case | Quantity | Expected |
| --- | --- | --- |
| H2 (and any homonuclear σ pair) | σ and σ\* | ε\_H ± ssσ (ssσ < 0, so ε\_H + ssσ is lower) |
| CH4, T\_d, H–H off | a1 block | \[\[ε\_2s, 2·ssσ\], \[2·ssσ, ε\_H\]\] |
| CH4 | t2 block (×3) | \[\[ε\_2p, (2/√3)·spσ\], \[(2/√3)·spσ, ε\_H\]\] |
| Four H in T\_d (CH4's LGOs alone), coupling t\_HH | LGO energies | a1: ε\_H + 3t\_HH; t2 (×3): ε\_H − t\_HH (matches the counting rule R6) |
| BH3, D3h, H–H off | a1′ block | \[\[ε\_2s, √3·ssσ\], \[√3·ssσ, ε\_H\]\] |
| BH3 | e′ block (×2) | \[\[ε\_2p, √(3/2)·spσ\], \[√(3/2)·spσ, ε\_H\]\] |
| BH3 | a2″ | ε\_2p, nonbonding, empty: the LUMO |
| Six H in O\_h, cis t\_c, trans t\_t | LGO energies | a1g: ε\_H + 4t\_c + t\_t; t1u (×3): ε\_H − t\_t; eg (×2): ε\_H − 2t\_c + t\_t |
| Allyl cation | π bond order per C–C | counting 0.5; Coulson 1/√2 ≈ 0.707 |
| H3+ | H–H bond order | counting 1/3; Coulson 2/3 |
| BH3 | B–H bond order | counting 1 |
| Homonuclear diatomic (B2, C2, N2, O2, F2) | σg, σu, πu, πg blocks | σg: D(ε\_s + ssσ, ε\_p − ppσ; −spσ); σu: D(ε\_s − ssσ, ε\_p + ppσ; +spσ); πu: ε\_p + ppπ (×2); πg: ε\_p − ppπ (×2) |
| CO | σ and π blocks | σ: 4×4 in (s\_C, pz\_C, s\_O, pz\_O), CONVENTIONS A.6; π: D(ε\_pC, ε\_pO; ppπ) (×2); counting bond order 3; HOMO σ, carbon-weighted |
| Diatomic orderings | Reference, asserted as a fidelity comparison | B2, C2, N2: πu below σg(2p); O2, F2: σg(2p) below πu; a disagreement is shown, never tuned away |

The six-H energies sum to 6ε\_H, a quick trace check. With t\_c, t\_t < 0 and |t\_c| larger, the order is a1g < t1u < eg.

Signs of off-diagonal elements depend on SALC phase and may flip; eigenvalues may not.

**Property tests** (no reference values needed):

- Zero coupling returns the separate central-atom and ligand levels.
- Number of MOs equals basis size for every configuration.
- A random orthogonal rotation inside any degenerate eigenspace leaves energies, degenerate-summed compositions and bonding classes unchanged.
- Flipping the sign of any single eigenvector changes no displayed quantity.
- Small distortion splits degeneracies exactly as the descent-in-symmetry tables say.
- Weak-coupling d-block splitting matches 3e\_σ (O\_h σ) and (4/3)e\_σ (T\_d σ) to first order (Release 4, advanced complex mode).

**Mutation check.** The checker deliberately flips one angular-factor sign or permutes two ligands, and confirms at least one test fails. A suite that survives the mutation is too weak.

**Process.** One command (`make check`) runs everything. Snapshot tests compare DiagramSpec JSON, not pixels. Disagreements between builder, checker and oracles go to the owner; neither AI may edit CONVENTIONS.md or the oracle file.

**Independent model cross-check (CrystOD).** The checker's environment includes a pinned copy of CrystOD, supplied as a zip; record its version and file hash. It is never a dependency of the shipped engine: it needs pymatgen, phonopy, spglib and SciPy, which the browser build cannot carry, and it uses a different model (extended Hückel with overlap). Use only its default route; its PySCF route is out of scope (section 13).

| Compare | Expectation |
| --- | --- |
| Point group detected from the same geometry | Identical |
| Degeneracy pattern of the MOs | Identical for irreps both bases share (CrystOD may include orbitals this app omits, e.g. all F valence orbitals in SF6) |
| Irrep labels | Identical after mapping axis conventions; C2v b1/b2 depend on the molecular plane (CONVENTIONS.md §5) |
| Qualitative order within each irrep, compared through eigenvector-derived descriptors (bonding character and dominant fragment), after excluding CrystOD levels that lie mostly outside this app's basis (weight below 0.5); HOMO and LUMO irreps | Agreement expected; a mismatch goes to the owner rather than failing automatically, because the models differ |
| Energies, gaps, compositions, bond orders | Not compared: different models and parameters |

Deferred until the programme is finished (owner decision B1); when run (`make check-release`), run it on every Release 1 case, and on the matcher fixtures from Release 3. A label or degeneracy mismatch blocks release until explained. CrystOD's own correctness is not assumed: where it disagrees with the analytic oracles, the oracles win.

## 13. Out of scope

Excluded from every release unless the release plan (section 6) says otherwise, so the site never implies predictive power the model lacks.

- Uploading 3D structure files; the matcher (section 3) takes formulas, names and SMILES and uses curated or VSEPR geometries.
- Transition-metal complexes in method mode; the engine supports them in an advanced mode that the course method does not cover.
- DFT or any ab initio comparison.
- Multiplets, term symbols, Racah parameters and spin–orbit coupling.
- Assigning parameters automatically from an element's or ligand's name beyond the curated presets.
- Nonorthogonal overlap (Hc = ESc); planned as a later, separately validated option.
- Point groups with complex-conjugate irrep pairs (C3, C3h, S4, T and others) until section 14's real-form handling lands.
- Any language-model call at runtime.

## 14. Extensibility

After the Slater–Koster generalisation, most later additions are data, not code. The earlier spec would have locked the tool in: it assumed every system has a central atom and that every coupling is a central–ligand AOM factor. Cyclic H3 and allyl break both assumptions, and Slater–Koster removes both while reproducing the AOM exactly.

| Addition | Cost | Why |
| --- | --- | --- |
| More complexes in existing groups (trans-MA4B2, D4h) | Data only | Preset, parameters and oracle rows |
| New groups with real irreps (C4v, C3v, D3d, D5h) | Data file plus tests | Basis representation is built generically (F2.4); unlocks NH3, square-pyramidal and trigonal-bipyramidal molecules, cis and fac/mer complexes, both conformations of ferrocene |
| Larger π systems | Data only | Hückel is the p\_π-only special case already used for allyl |
| Linear molecules (CO2, HF) | Done for the Release 1 diatomics | D∞h and C∞v run in D2h or C2v, with a lookup mapping labels to σg, πu and so on; polyatomic linear cases reuse it |
| Groups with complex characters (C3, S4, T) | Moderate code | Conjugate irrep pairs must be combined into real two-dimensional representations before projection |
| User-entered geometry | Moderate code | Needs tolerance-based point-group detection; MolSym or libmsym earn their place here |
| Nonorthogonal overlap / extended Hückel | Large | Needs overlap integrals, a generalised eigensolver, and a redefinition of composition (Mulliken vs Löwdin) |
| Many-electron states | Large, separate module | Multiplets and spin states are a different model layer |
| Periodic solids | Not this project | Band structure needs k-space machinery |

**Adding molecules: molecule packs.** Every molecule, template and point group lives in a small data file, and one command validates it. Adding a molecule never touches engine code.

A molecule pack declares:

- **Identity:** name, formula, charge, SMILES, aliases the matcher should recognise.
- **Structure:** coordinates, or a template id plus bond lengths and angles.
- **Basis and counting:** orbitals per atom; electron-counting convention.
- **Parameters:** "use scheme" (F2.12) by default, or curated overrides with sources.
- **Point group:** declared, or detected and cross-checked; C1 if absent from the library.
- **Optional extras:** lecture-simplification switches, narration overrides, morph paths to existing cases (new case-map edges), curated fidelity checks with citations.

`make add-molecule <pack>` runs the gate: schema check, basis count (R4), engine run, ledger totals against eigenvalues, invariant checks along any declared morphs, the generated pairing schedule (F2.13), and the matcher fixture. A pack that passes appears in the app; one that fails names the failing check.

| Kind of addition | Effort |
| --- | --- |
| New member of an existing template | One pack, no code |
| New template of equivalent ligands, group in library | One pack plus a template entry, no code |
| New point group with real characters | One group data file (passes closure and orthogonality checks), no code |
| New element | One row in the parameter table with its source |
| Groups with complex characters, new orbital types beyond f, new coupling kinds | Engine code |

Packs are plain text, so an AI agent can draft one from a prompt ("add NH3 as a C3v member of a new pyramidal template"); the gate, not the agent, decides whether it ships. New molecules start with amber fidelity chips until someone adds cited checks.
