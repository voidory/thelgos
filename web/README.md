# Web front end (Release 1 subset, provisional)

Static files only. The page runs the Python engine from `engine/` in Pyodide 0.27.7 inside a Web Worker (`worker.js`), reading `params/`, `groups/` and `CONVENTIONS.md` straight from the repository. So the browser runs exactly the files `make check` tests. `app.js` draws the DiagramSpec it receives and computes no chemistry.

To run it, serve the repository root and open `/web/`:

```bash
python3 -m http.server 8000
```

## Look

The approved layout (case rail, construction timeline, tabbed inspector) is the default; `?design=classic` shows the first scaffold. The MO diagram uses short bars, STIX Two Text labels in standard notation, drawn electron arrows, an energy axis and collision-free labels; see `drafts/design-proposal.md`.

## What is built

| Spec | Built |
|---|---|
| L2 | Structure view: atoms in a fixed 2D projection, orbital lobes coloured by phase (red positive, blue negative), one contact line per ledger contact row (solid stabilising, dashed destabilising, width by size); hovering a row or a line highlights both |
| L3 | Rule badges (R1–R10) on each construction stage, the ledger, matching and bond order; tapping one shows the rule, its meaning in the model and its badge |
| P3 | Why-jump: "Show where this level appeared" moves the slider to the checkpoint at which the selected level first reached its final energy, with that moment's ledger and structure |
| P6 | Matching table: ⟨central orbital\|H\|LGO⟩ tallies; zero means no interaction, and the zero pattern follows the irreps |
| F3.7 | 3D view (WebGL2): atoms and the selected orbital's lobes (p functions on an atom merged into one oriented lobe); opens in the 2D sketch's orientation; drag to rotate |
| F5.1–F5.6 | Provenance tooltips; generated explanation; SVG, PNG and JSON export; JSON import with validation; share links (compressed config in the URL fragment, nothing stored on a server) |
| L1 | Contact ledger for the selected level: every on-site and contact row, phase, contact weight, energy contribution, running total. The total equals the eigenvalue |
| P1–P2 | Construction slider over H(s), with stage ticks and narration, and B5 notes for skipped stages |
| P7 | Filling arrows; counting and density-matrix bond orders, side by side |
| Watch mode | Plays the construction stage by stage |
| F3.1–F3.4 | Three-column diagram: central-atom orbitals, MOs, ligand group orbitals; diatomics in the textbook layout (atom 1 orbitals, MOs, atom 2 orbitals, with σg/πu* labels); two columns for other systems without a central atom; symmetry-adapted fragment levels; dashed allowed and solid mixing layers, each with its own toggle; axis labelled in model units; electron arrows |
| F3.9 | Usable at 360 px with no horizontal scroll |
| — | "About this diagram" panel: reading notes (`params/case_notes.json`, explanatory text written for the app) and comparisons with curated experimental data (`params/curated_checks.json`, cited), computed by `engine/lgo_engine/fidelity.py` with the checker's rules. Groundwork for the spec §7 fidelity bar |
| F5.7 | Scope banner: model units, preset tag, not a spectral prediction, engine and CONVENTIONS versions |

## Browser tests

`make test-web` (needs Node, Playwright and network access to cdn.jsdelivr.net) runs `tests/web_test.js` in headless Chromium against the real Pyodide runtime. It checks:
- the runtime is Python 3.12.7 with NumPy 2.0.2;
- every fixture's energies equal `evidence/reference_spectra.json` to 1e−9 eV;
- warm recompute stays under 100 ms. This is measured on the build machine, not the reference phone;
- the ledger total equals the eigenvalue;
- the slider and its four SF6 stages, and the CH4 B5 skip note;
- L2 one-to-one ledger rows and structure lines, L3 badges, F5.1 tooltips, P6 zero tallies, P3 why-jump for SF6 eg (checkpoint 3 of 4);
- JSON export, re-import and rejection of an invalid file (failing field named); a share link opened in a new page; the 3D view painting;
- the reading panel: N2's ordering disagreement and CO's note on its starred HOMO; Si2 keeps its element symbol in the column titles and compares spin only; Al2 shows its spin disagreement and the on-site table note;
- the rail lists H2O and SF6 once each; the model/lecture toggle switches presets, is remembered, and hides for single-preset cases;
- Watch mode plays through; step jumps glide, and with reduced motion they do not;
- the selected level's sketch card sits clear of bars and labels;
- accessibility: every level is a named button and is listed in the screen-reader table; the inspector tabs answer to the arrow keys; Enter selects a focused level and keeps focus on it;
- there is no horizontal scroll at 360 px, and no page errors.

`make check` does not run these, because they need the network.

## Not built yet

The Release 2+ features (L4–L5, P4–P5, Practice and Exam modes, D1–D10), and from the design pass: highlighting ledger contacts in the 3D view, and pinning the timeline on phones.
