# What signing off CONVENTIONS v1.4 approves

Signing off (decision B6) freezes `CONVENTIONS.md` v1.4. From then on, any change needs a new version and your approval. This page lists what you would be approving. The reasons for each decision are in `DECISIONS.md`.

## The model

- **Orbital energies:** every element takes its starting orbital energies from Miessler, Fischer & Tarr's table (via LibreTexts, snapshot packaged). *Your decision M7.*
- **Couplings:** how strongly two orbitals interact comes from Harrison's distance rule, V = η·7.62/d² (η = −1.40, +1.84, +3.24, −0.81). *Delegated M2.*
- **No overlap:** bonding and antibonding split symmetrically. The app explains this.
- **Hydrogen factor:** couplings to hydrogen are multiplied by 0.6. *Delegated, then your instruction: follow Miessler (B2).*

## Decisions you made

| | Decision |
|---|---|
| B1 | CrystOD cross-check deferred to the end of the programme |
| B3 | Uranium sandwich: rings flat, above and below U |
| B4 | N₂ included |
| M5 | Period-3 diatomics Al₂–Cl₂ included |
| M6 | One entry per molecule, with a model/lecture toggle |
| M7 | All orbital energies from Miessler |

## Decisions I made under delegation (review if you want to)

| | Decision | Why it matters |
|---|---|---|
| B2 | Hydrogen factor 0.6 | Chosen to reproduce textbook orderings, so those orderings are not independent evidence |
| B3 | Sandwich atoms on ±x within the ring plane | Labelling convention only |
| B5 | Construction steps with zero coupling are skipped, with a note | Display only |
| M1 | B₂, C₂ and N₂ keep their wrong σ/π order, shown in red, not tuned away | The model's known failure, shown honestly |
| M3 | SF₆ fluorine–fluorine couplings at full strength | Spreads the fluorine combinations widely |
| M4 | SF₆ "lecture" preset is identical to the model | The app says so |

## Things resting on something weaker than a packaged source

1. **The upper limit on the hydrogen factor.** The four orderings that bound it give 0.4235 to 0.6668. The upper limit comes from water's LUMO being 3a₁* rather than 2b₂*, and CH₄'s being a₁* rather than t₂*. I recalled both from Miessler's figures; please check them in the book.
2. **CH₄'s second photoelectron band "near 23 eV".** This comes only from your spec. The first band (12.61 eV adiabatic, 13.60 eV vertical) is from NIST CCCBDB.
3. **Al₂'s ground state.** CCCBDB lists it as ³Σg with degeneracy 4, which is not self-consistent. Every candidate state is a triplet, so the spin comparison stands.

## Known model failures, shown in the app, not hidden

- B₂, C₂, N₂: wrong σg/πu order.
- B₂, C₂: wrong spin.
- Al₂: wrong spin.
- Absolute energies: compared with a measured spectrum only for CH₄, whose first band the model places about 7 eV too deep. The app tells students to read patterns, not positions.

## Not covered by this sign-off

- The CrystOD cross-check (deferred, B1).
- Release 2 and later features.
- Matching the "lecture" presets to your actual slides (no slides checked yet).
