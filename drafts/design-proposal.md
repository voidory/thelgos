# Design proposal for the LGO Method Explorer

**Status:** proposal for the owner. Preview it at `/web/?design=proposal`; the default look is unchanged until you approve. Screenshots were sent in the session.

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
