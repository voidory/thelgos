# Design proposal for the LGO Method Explorer

**Status:** layout approved by the owner (2026-11-06) and now the default; the old look stays at `/web/?design=classic`. The MO diagram was then redrawn after the owner's critique (section "Diagram redesign").

## Principles

1. **The diagram is the page.** It gets the widest column and the most contrast; everything else supports reading it.
2. **The lecture's procedure is the navigation.** The construction is a row of named steps ("Separate atoms → Fix one, flip partner → Combine equal-energy sets → Put the central atom back"), not an anonymous 0–4 slider.
3. **Red and blue mean phase, and nothing else** (spec L2). The interface uses one neutral accent (teal) for selection, buttons and badges, and amber for warnings.
4. **One inspector, in tabs:** Structure, Contact ledger, Matching, Bond order, Export. This replaces the long scrolling column.
5. **Numbers line up.** Tabular figures everywhere, so energies and ledger columns read down cleanly.

## Layout

| Region | Desktop (≥ 900 px) | Phone |
|---|---|---|
| Header | Title plus a one-line scope banner | Same, wrapped |
| Left rail | Cases grouped by family: central atom with ligands, diatomics, rings and chains, ligand sets | Horizontally scrolling chips |
| Stage | Construction timeline (play, named steps, fine slider), narration and rule badges, then the diagram with an energy axis and gridlines, then "About this diagram" | Full width, in that order |
| Inspector | Sticky card with tabs | Below the stage, same tabs |

## Already in the preview

- Grouped case rail, named construction steps, accent-coloured play button.
- Energy axis with 5 eV gridlines and labels.
- Tabbed inspector.
- Calmer palette, tabular numerals, rounded cards, thicker level strokes.

## Still to do in the design pass, if you approve the direction

1. **Label placement:** close-lying levels (O2's σg and σu*, CO's σ* and π) currently overlap. Add simple collision avoidance in the right-hand label gutter.
2. **Orbital thumbnails:** a small lobe sketch next to each MO, coloured by phase, so the diagram reads like the lecture's drawings. The data already exists in the structure spec.
3. **Phone diagram:** taller aspect and larger type on small screens, with the timeline pinned while scrolling.
4. **Slider animation:** levels glide between plateaus instead of redrawing.
5. **3D view:** proper ellipsoids along arbitrary directions, and highlighting the item selected in the diagram.
6. **Accessibility pass:** focus styles, keyboard access for every control, and contrast checks in both light and dark modes.

## Questions for you

- Does teal as the single interface colour work for you, or should it match your course or department colours?
- Should the lecture-vs-model choice (H2O, SF6) be a toggle on one case instead of two separate entries in the list?

## Diagram redesign (after the owner's critique)

The owner asked for larger type, shorter level lines, and a more "physics" look. Changes:

| Before | After | Why |
|---|---|---|
| Level lines stretched across the whole column | Short bars, 30 px each; a degenerate set is a row of bars | Textbook convention; degeneracy is read by counting bars, so the ×n tags are gone |
| 11 px system font | STIX Two Text (the typeface of scientific publishing) at 17 px, with real subscripts and superscripts: 1a₁g, 2t₁u*, 1πg* | Readable at a glance, and looks like the lecture notes and textbooks |
| Irrep names as typed (A1g) | Mulliken lower case for orbitals, numbered within each symmetry (1a₁g, 2a₁g), valence only, per CONVENTIONS §6 | Standard MO notation |
| Text arrows ↑↓ | Drawn half-headed electron arrows | How electrons are drawn by hand |
| Mixing lines dark | Lighter grey under the bars; thickness still ∝ weight | Lines support the levels instead of competing with them |
| Labels could collide | Labels at least 17 px apart in each gutter, with leader lines; text halos so lines never cross a label | Readability where levels crowd together |
| Near-degenerate levels drawn on top of each other (water's 2b₂*/3a₁*, 0.1 eV apart) | Levels within 7 px share a row side by side, each still at its own height | Both stay visible; the energy scale is untouched |
| HOMO/LUMO not marked | Small accent tags after the label | Frontier orbitals at a glance |
| No axis | A quiet energy axis with ticks, "E / eV (model)" | Keeps the to-scale drawing honest |
| Phone: whole diagram shrunk to half size | A compact geometry below 560 px with shorter bars and tighter columns, so labels stay at reading size | The phone view is now usable |

The energies stay to scale (spec F3.1): the readability problems came from line length and type, not from the scale.
