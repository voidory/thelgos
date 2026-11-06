# Web front end (Release 1 subset, provisional)

Static files only. The page runs the Python engine from `engine/` in Pyodide 0.27.7 inside a Web Worker (`worker.js`), reading `params/`, `groups/` and `CONVENTIONS.md` straight from the repository. So the browser runs exactly the files `make check` tests. `app.js` draws the DiagramSpec it receives and computes no chemistry.

To run it, serve the repository root and open `/web/`:

```bash
python3 -m http.server 8000
```

## What is built

| Spec | Built |
|---|---|
| L1 | Contact ledger for the selected level: every on-site and contact row, phase, contact weight, energy contribution, running total. The total equals the eigenvalue |
| P1–P2 | Construction slider over H(s), with stage ticks and narration, and B5 notes for skipped stages |
| P7 | Filling arrows; counting and density-matrix bond orders, side by side |
| Watch mode | Plays the construction stage by stage |
| F3.1–F3.4 | Three-column diagram: central-atom orbitals, MOs, ligand group orbitals; diatomics in the textbook layout (atom 1 orbitals, MOs, atom 2 orbitals, with σg/πu* labels); two columns for other systems without a central atom; symmetry-adapted fragment levels; dashed allowed and solid mixing layers, each with its own toggle; axis labelled in model units; electron arrows |
| F3.9 | Usable at 360 px with no horizontal scroll |
| F5.7 | Scope banner: model units, preset tag, not a spectral prediction, engine and CONVENTIONS versions |

## Browser tests

`make test-web` (needs Node, Playwright and network access to cdn.jsdelivr.net) runs `tests/web_test.js` in headless Chromium against the real Pyodide runtime. It checks:
- the runtime is Python 3.12.7 with NumPy 2.0.2;
- every fixture's energies equal `evidence/reference_spectra.json` to 1e−9 eV;
- warm recompute stays under 100 ms. This is measured on the build machine, not the reference phone;
- the ledger total equals the eigenvalue;
- the slider and its four SF6 stages, and the CH4 B5 skip note;
- Watch mode plays through;
- there is no horizontal scroll at 360 px, and no page errors.

`make check` does not run these, because they need the network.

## Not built yet

L2–L3 contact lines and rule badges on the structure, P3 why-jump, P6 per-lobe matching view, the F3.7 3D viewer, export (section 11), and the Release 2+ features.
