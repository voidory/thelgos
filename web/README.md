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
| F3.7 | 3D view (WebGL2): atoms, ledger contacts as rods, and the selected orbital's lobes (p functions on an atom merged into one oriented lobe); opens in the 2D sketch's orientation; drag to rotate. Any MO or fragment orbital (LGO, central-atom or atom orbital) can be selected in the diagram; the ledger, 2D structure, sketch card and 3D view follow it; hovering a ledger row highlights that contact in 2D and 3D |
| D2, D4, D5, D6, F3.6, D10 | Release 2 Milestone A: approximation slider (lecture ↔ model), rotation within a degenerate set, water's lone-pair dial, electron stepper with charge, symmetry-only view, and an in-app speed test (footer) |
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
- the reading panel: N2's ordering disagreement and CO's note on its starred HOMO; Si2 keeps its element symbol in the column titles and compares spin only; Al2 shows its spin disagreement and its reading note;
- the rail lists H2O and SF6 once each; the model/lecture toggle switches presets, is remembered, and hides for single-preset cases;
- gate C4, on every case: the s = 0 end point equals the on-site table (1e−12 eV); the ledger identity holds at s = 1 and s = 0.5 (1e−10 eV); counting bond orders equal CONVENTIONS A.5;
- Watch mode plays through; step jumps glide, and with reduced motion they do not;
- the selected level's sketch card sits clear of bars and labels;
- accessibility: every level is a named button and is listed in the screen-reader table; the inspector tabs answer to the arrow keys; Enter selects a focused level and keeps focus on it;
- F3.7: a ligand group orbital is selectable with its own ledger; hovering a ledger contact changes the 3D view;
- Release 2 Milestone A: D6 electron stepper (H3+ + e⁻ is neutral H3 with the Jahn–Teller warning); D2 approximation slider passes lecture → model; D5 lone pairs equivalent at 45°; D4 dial keeps the energy; F3.6 symmetry-only view keeps labels and lines; D10 speed test runs;
- the CH4 photoelectron spectrum: two bands 3:1, badged approximate, measured positions cited;
- the banner shows the CONVENTIONS version read from the file;
- layout: the full diagram from 1024 px up, the phone geometry at 360 px, and nothing clipped or scrolling sideways;
- there is no horizontal scroll at 360 px, and no page errors.

`make check` does not run these, because they need the network.

## Not built yet

The Release 2+ features (L4–L5, P4–P5, Practice and Exam modes, D1–D10), from Release 1: nothing. Share links and exports carry a selected MO but not a selected fragment orbital.
