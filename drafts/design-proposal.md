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

## Design pass: status

1. **Label placement:** done (pixel-proximity rows and gutter collision avoidance).
2. **Orbital sketch:** done, as a single card for the *selected* level rather than a thumbnail on every MO. Thumbnails on all levels crowd the diagram (SF6 has 10 levels; the cards would overlap the correlation lines). The card goes in the emptiest free spot near its level and is left out when none exists.
3. **Phone diagram:** done (compact geometry; the timeline stays pinned while the diagram scrolls; the canvas widens instead of clipping long labels). Mid-size screens (900–1279 px) use a narrower rail and inspector, or put the inspector below the diagram, so they get the full diagram rather than the phone one.
4. **Glide:** done for step jumps and "why does it jump". Every frame is a real engine evaluation of H(s), not a screen interpolation, so the motion is the model's. Reduced-motion users get the jump.
5. **3D view:** ellipsoids now follow any direction. The 3D view opens in the same orientation as the 2D sketch. Ledger contacts are drawn as rods, and the hovered one is highlighted, as in the 2D sketch. Fragment orbitals (LGOs) are selectable too.
6. **Accessibility:** done.
   - Focus rings everywhere.
   - Every MO level is a named button (label, energy, electrons, HOMO/LUMO) and keeps focus across redraws.
   - Screen readers get a text table of the levels.
   - The inspector uses WAI-ARIA tabs (arrow keys, Home and End).
   - Reduced-motion is respected.
   - Contrast: every text token is at least 4.5:1 on the panel and the page background in both themes. The light theme's dashed "allowed" line was 2.9:1 and is now 3.6:1, above the 3:1 needed for graphics.

**Correction made along the way:** the structure sketch and the 3D view drew each p function separately, so a level with p_x and p_y on one atom showed two crossed dumbbells. That is not the shape of the orbital. The p functions on one atom add as a vector, giving a single p orbital along (c_x, c_y, c_z) with size |c|, and the drawings now show that. The ledger and contact lines still list the basis functions separately.

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
