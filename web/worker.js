// Web Worker: runs the Python engine in Pyodide 0.27.7 (CONVENTIONS tech stack pin).
// The engine and data are fetched from the repository itself, so the browser runs the same files
// that `make check` tests. Nothing is computed in JavaScript.
importScripts("https://cdn.jsdelivr.net/pyodide/v0.27.7/full/pyodide.js");

const ENGINE = ["__init__.py", "model.py", "params.py", "hamiltonian.py", "solve.py", "symmetry.py",
  "analysis.py", "ledger.py", "result.py", "schedule.py", "cases.py", "diagram.py", "fidelity.py", "views.py", "app_api.py"];
const GROUPS = ["C2h", "C2v", "D2d", "D2h", "D3h", "D4h", "D6h", "D8h", "Oh", "Td"];
const DATA = ["params/onsite.csv", "params/pairs.csv", "params/presets.json", "params/channel_defaults.json",
  "params/case_notes.json", "params/curated_checks.json",
  "CONVENTIONS.md", "engine/load.py", ...GROUPS.map(g => `groups/${g}.json`),
  ...ENGINE.map(f => `engine/lgo_engine/${f}`)];

let py = null;
const ready = (async () => {
  const t0 = performance.now();
  py = await loadPyodide();
  await py.loadPackage("numpy");
  const base = new URL("../", self.location.href);
  for (const path of DATA) {
    const r = await fetch(new URL(path, base));
    if (!r.ok) throw new Error(`fetch ${path}: ${r.status}`);
    const bytes = new Uint8Array(await r.arrayBuffer());
    const dir = "/bundle/" + path.split("/").slice(0, -1).join("/");
    py.FS.mkdirTree(dir);
    py.FS.writeFile("/bundle/" + path, bytes);
  }
  py.runPython(`
import sys, json
sys.path[:0] = ["/bundle/engine"]
from load import load_params, load_groups, conventions_sha256, conventions_version, load_extras
from lgo_engine import app_api, __version__
P = load_params("/bundle"); G = load_groups("/bundle"); CSHA = conventions_sha256("/bundle"); X = load_extras("/bundle")
`);
  return { load_ms: performance.now() - t0, version: py.runPython("__version__"),
           python: py.runPython("sys.version.split()[0]"), conventions: JSON.parse(py.runPython("json.dumps(conventions_version('/bundle'))")), numpy: py.runPython("import numpy; numpy.__version__") };
})();

self.onmessage = async (ev) => {
  const { id, cmd, args } = ev.data;
  try {
    const info = await ready;
    let out;
    if (cmd === "info") out = info;
    else if (cmd === "cases") out = JSON.parse(py.runPython("json.dumps(app_api.case_list(P))"));
    else if (cmd === "view") {
      py.globals.set("ARGS", JSON.stringify(args));
      const t0 = performance.now();
      out = JSON.parse(py.runPython(`
a = json.loads(ARGS)
json.dumps(app_api.view(P, G, a["case"], s=a.get("s", 1.0), level=a.get("level"), conventions_sha256=CSHA, extras=X, fragment=a.get("fragment"), electrons=a.get("electrons"), approx=a.get("approx"), dial=a.get("dial") or 0.0, lone_pair_angle=a.get("lone_pair_angle")), allow_nan=False)`));
      out.compute_ms = performance.now() - t0;
    } else if (cmd === "export") {
      py.globals.set("ARGS", JSON.stringify(args));
      out = JSON.parse(py.runPython(`
a = json.loads(ARGS)
json.dumps(app_api.export(P, G, a["case"], s=a.get("s", 1.0), level=a.get("level"), conventions_sha256=CSHA), allow_nan=False)`));
    } else if (cmd === "import") {
      py.globals.set("ARGS", JSON.stringify(args));
      out = JSON.parse(py.runPython(`json.dumps(app_api.validate_import(P, json.loads(ARGS), conventions_sha256=CSHA))`));
    } else throw new Error(`unknown command ${cmd}`);
    self.postMessage({ id, ok: true, out });
  } catch (e) {
    self.postMessage({ id, ok: false, error: String(e) });
  }
};
