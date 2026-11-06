// Main thread: draws what the engine returns and computes no chemistry (spec section 10).
"use strict";
const worker = new Worker("worker.js");
let seq = 0;
const pending = new Map();
worker.onmessage = (ev) => {
  const { id, ok, out, error } = ev.data;
  const p = pending.get(id);
  pending.delete(id);
  ok ? p.resolve(out) : p.reject(new Error(error));
};
const call = (cmd, args) => new Promise((resolve, reject) => {
  const id = ++seq;
  pending.set(id, { resolve, reject });
  worker.postMessage({ id, cmd, args });
});

const $ = (id) => document.getElementById(id);
const NS = "http://www.w3.org/2000/svg";
const state = { case: null, s: 1, level: null, view: null, info: null, watching: false };

function el(tag, attrs = {}, text) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (text !== undefined) e.textContent = text;
  return e;
}

const fmt = (x, d = 3) => (x >= 0 ? "+" : "−") + Math.abs(x).toFixed(d);
// D∞h/C∞v names in the textbook form (σg, πu*, …); the star marks an antibonding level as the engine classifies it.
const GREEK = { sigma: "σ", pi: "π" };
const irrepText = (lv) => {
  if (!lv.linear_label) return lv.irrep;
  const [base, par] = lv.linear_label.split("_");
  return `${GREEK[base] || base}${par || ""}${lv.character === "antibonding" ? "*" : ""}`;
};

function draw(view) {
  const svg = $("diagram");
  svg.replaceChildren();
  const d = view.diagram;
  const W = 720, H = 540, top = 36, bottom = 30;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const atEnd = view.schedule.s >= 1;
  const centreLevels = atEnd ? d.centre : d.slider_levels;
  const all = [...d.left, ...d.right, ...d.centre, ...d.slider_levels].map((x) => x.energy_eV);
  let lo = Math.min(...all), hi = Math.max(...all);
  const pad = (hi - lo) * 0.06 || 1;
  lo -= pad; hi += pad;
  const y = (e) => top + (hi - e) / (hi - lo) * (H - top - bottom);
  const cols = d.columns;
  const colX = cols.length === 3 ? { left: [95, 195], centre: [295, 425], right: [515, 615] }
                                 : { right: [80, 220], centre: [380, 520] };
  for (const c of cols) svg.append(el("text", { x: (colX[c][0] + colX[c][1]) / 2, y: 18, "text-anchor": "middle", class: "muted" }, d.column_titles[c]));
  svg.append(el("text", { x: 12, y: H / 2, transform: `rotate(-90 12 ${H / 2})`, "text-anchor": "middle", class: "muted" }, d.axis_label));

  const pos = {};
  const showIrreps = $("t-irreps").checked;
  // Layout only: levels of one column at the same energy (e.g. CH4's a1 and t2 H combinations with
  // H-H coupling off) share the row side by side instead of being drawn on top of each other.
  const slots = {};
  function slotFor(col, list) {
    const rows = [];
    list.forEach((lv) => {
      const r = rows.find((r) => Math.abs(r.e - lv.energy_eV) < 1e-9);
      if (r) { r.items.push(lv); r.m += lv.multiplicity; } else rows.push({ e: lv.energy_eV, items: [lv], m: lv.multiplicity });
    });
    for (const r of rows) { let k = 0; for (const lv of r.items) { slots[col + ":" + (lv.id || lv.energy_eV)] = { start: k, total: r.m, last: lv === r.items[r.items.length - 1], row: r }; k += lv.multiplicity; } }
  }
  slotFor("left", d.left); slotFor("right", d.right);
  function levelRow(lv, col, cls, key) {
    const [x0, x1] = colX[col];
    const sl = slots[col + ":" + (lv.id || lv.energy_eV)] || { start: 0, total: lv.multiplicity, last: true, row: { items: [lv] } };
    const m = lv.multiplicity, M = sl.total, gap = 6, w = (x1 - x0 - gap * (M - 1)) / M, yy = y(lv.energy_eV);
    const g = el("g");
    for (let i = 0; i < m; i++) {
      const a = x0 + (sl.start + i) * (w + gap);
      g.append(el("line", { x1: a, x2: a + w, y1: yy, y2: yy, class: `level ${cls}` }));
      if (lv.electrons !== undefined && atEnd) {
        // One electron per line first, then pair up (the usual drawing; it does not assert a spin state).
        const n = (i < lv.electrons ? 1 : 0) + (i < lv.electrons - m ? 1 : 0);
        const arrows = n === 2 ? "↑↓" : n === 1 ? "↑" : "";
        if (arrows) g.append(el("text", { x: a + w / 2, y: yy - 3, "text-anchor": "middle" }, arrows));
      }
    }
    const segX0 = x0 + sl.start * (w + gap), segX1 = segX0 + m * w + (m - 1) * gap;
    if (sl.last) {
      const items = sl.row.items;
      const mult = items.map((it) => (it.multiplicity > 1 ? `×${it.multiplicity}` : "")).filter(Boolean).join(" ");
      const irr = showIrreps && cls !== "moving" ? items.map((it) => (it.irrep && it.irrep !== "—" ? irrepText(it) : "")).filter(Boolean).join(" + ") : "";
      const txt = [mult, irr].filter(Boolean).join(" ");
      if (txt) g.append(el("text", { x: col === "left" ? x0 - 4 : x1 + 4, y: yy + 4, "text-anchor": col === "left" ? "end" : "start", class: irr ? "" : "muted" }, txt));
    }
    if (lv.label && col !== "centre") g.append(el("text", { x: (segX0 + segX1) / 2, y: yy + 14, "text-anchor": "middle", class: "muted" }, lv.label));
    pos[key] = { x0: segX0, x1: segX1, y: yy };
    return g;
  }
  for (const lv of d.left) svg.append(levelRow(lv, "left", "frag", lv.id));
  for (const lv of d.right) svg.append(levelRow(lv, "right", "frag", lv.id));
  centreLevels.forEach((lv, k) => {
    const g = levelRow(lv, "centre", atEnd ? (state.level === k ? "selected" : "") : "moving", atEnd ? `M${k}` : `S${k}`);
    if (atEnd) {
      const [x0, x1] = colX.centre, yy = y(lv.energy_eV);
      const hit = el("rect", { x: x0, y: yy - 7, width: x1 - x0, height: 14, class: "mo-hit", "data-mo": k, "data-energy": lv.energy_eV, tabindex: 0 });
      hit.addEventListener("click", () => select(k));
      hit.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") select(k); });
      g.append(hit);
      if (lv.warnings && lv.warnings.length) g.append(el("text", { x: x0 - 6, y: yy + 4, "text-anchor": "end", class: "warn" }, "!"));
    }
    svg.append(g);
  });
  if (atEnd) {
    for (const ln of d.lines) {
      const layerOn = ln.layer === "mixing" ? $("t-mixing").checked : $("t-allowed").checked;
      if (!layerOn || !pos[ln.from] || !pos[ln.to]) continue;
      const f = pos[ln.from], t = pos[ln.to];
      const fromRight = f.x0 > t.x0;
      const xa = fromRight ? f.x0 : f.x1, xb = fromRight ? t.x1 : t.x0;
      svg.append(el("line", { x1: xa, y1: f.y, x2: xb, y2: t.y, class: ln.layer === "mixing" ? "mix" : "allowed",
                              "stroke-width": ln.layer === "mixing" ? (0.6 + 4 * ln.weight).toFixed(2) : 1 }));
    }
  }
}

function narrate(view) {
  const sc = view.schedule;
  const K = sc.stages.length;
  let text = "";
  if (sc.s <= 0) text = "Separate atoms: no couplings switched on.";
  else if (sc.active_stage !== null) text = `Stage ${sc.active_stage + 1} of ${K}: ${sc.stages[sc.active_stage].narration}`;
  else {
    const done = Math.round(sc.s * K);
    text = done >= K ? "Finished diagram." : `After stage ${done} of ${K}: ${sc.stages[done - 1].narration}`;
  }
  const n = $("narration");
  n.textContent = text;
  for (const sk of sc.skipped) {
    const s = document.createElement("span");
    s.className = "skip";
    s.textContent = `Skipped: ${sk.narration} ${sk.note}`;
    n.append(s);
  }
  const ticks = $("ticks");
  ticks.replaceChildren();
  sc.checkpoints.forEach((c, i) => {
    const t = document.createElement("span");
    t.style.left = `${c.s * 100}%`;
    t.textContent = i === 0 ? "0" : String(i);
    t.title = c.label;
    ticks.append(t);
  });
}

function ledger(view) {
  const box = $("ledger");
  box.replaceChildren();
  const L = view.ledger;
  if (!L) { $("ledger-hint").textContent = view.schedule.s >= 1 ? "Tap a molecular orbital to see every contact that sets its energy." : "Finish the construction (slider at the end) to inspect a level."; return; }
  $("ledger-hint").textContent = `Level ${L.level + 1}: ledger total ${fmt(L.total_eV, 6)} eV; eigenvalue ${fmt(L.eigenvalue_eV, 6)} eV.`;
  const t = document.createElement("table");
  t.dataset.total = L.total_eV;
  t.dataset.eigen = L.eigenvalue_eV;
  t.innerHTML = "<thead><tr><th>Contact</th><th>Phase</th><th>c·c·H</th><th>Energy</th><th>Running</th></tr></thead>";
  const tb = document.createElement("tbody");
  let run = 0;
  for (const r of L.rows) {
    if (Math.abs(r.energy_eV) < 1e-9) continue;
    run += r.energy_eV;
    const tr = document.createElement("tr");
    const name = r.kind === "onsite" ? L.basis_labels[r.mu] : `${L.basis_labels[r.mu]}–${L.basis_labels[r.nu]}`;
    tr.innerHTML = `<td>${name}</td><td>${r.kind === "onsite" ? "on-site" : (r.phase === "in" ? "in-phase" : "out-of-phase")}</td>` +
      `<td>${r.kind === "onsite" ? "" : fmt(r.contact_weight_eV)}</td><td>${fmt(r.energy_eV)}</td><td>${fmt(run)}</td>`;
    tb.append(tr);
  }
  t.append(tb);
  box.append(t);
}

function bonds(view) {
  const bo = view.result.bond_orders;
  const c = bo.counting, dm = bo.density_matrix;
  const vals = Object.values(dm.per_bond);
  const dmText = vals.length ? `${Math.min(...vals).toFixed(3)}${Math.max(...vals) - Math.min(...vals) > 1e-9 ? `–${Math.max(...vals).toFixed(3)}` : ""}` : "—";
  $("bonds").innerHTML = `<p>Counting: <b>${c.per_bond === null ? "—" : c.per_bond.toFixed(3)}</b> per bond ` +
    `(${c.bonding_electrons} bonding, ${c.antibonding_electrons} antibonding electrons over ${c.bonds} bonds).<br>` +
    `${dm.definition}: ${dmText}. These are different definitions and need not agree.</p>`;
}

function banner(view) {
  const tag = view.preset.provenance_tag;
  $("banner").textContent = `Model units (eV scale) · preset: ${tag} (${view.preset.id}) · not a spectral prediction · ` +
    `engine ${view.result.version} · CONVENTIONS v1.3, not frozen` + (view.result.warnings.length ? ` · warnings: ${view.result.warnings.join(", ")}` : "");
}

let latest = 0;
async function refresh() {
  const mine = ++latest;
  const view = await call("view", { case: state.case, s: state.s, level: state.s >= 1 ? state.level : null });
  if (mine !== latest) return state.view;  // a newer request superseded this one
  state.view = view;
  draw(view); narrate(view); ledger(view); bonds(view); banner(view);
  document.body.dataset.computeMs = view.compute_ms.toFixed(1);
  document.body.dataset.rendered = `${view.case}|${view.schedule.s}|${view.ledger ? view.ledger.level : ""}`;
  return view;
}

function select(k) { state.level = k; refresh(); }

async function watch() {
  if (state.watching) { state.watching = false; return; }
  state.watching = true;
  $("watch").textContent = "Stop";
  const K = state.view.schedule.stages.length;
  state.level = null;
  for (let k = 0; k < K && state.watching; k++) {
    for (let f = 0; f <= 20 && state.watching; f++) {
      state.s = (k + f / 20) / K;
      $("s").value = state.s;
      await refresh();
      await new Promise((r) => setTimeout(r, 30));
    }
    await new Promise((r) => setTimeout(r, 700));
  }
  state.s = 1; $("s").value = 1;
  await refresh();
  state.watching = false;
  $("watch").textContent = "Watch";
  document.body.dataset.watchDone = "1";
}

(async () => {
  try {
    const info = await call("info");
    state.info = info;
    $("footer").textContent = `Pyodide · Python ${info.python} · NumPy ${info.numpy} · first load ${(info.load_ms / 1000).toFixed(1)} s`;
    const cases = await call("cases");
    const sel = $("case");
    for (const c of cases) sel.append(new Option(`${c.name} (${c.group})`, c.name));
    state.case = cases.find((c) => c.name === "SF6 model").name;
    sel.value = state.case;
    sel.addEventListener("change", () => { state.case = sel.value; state.level = null; state.s = 1; $("s").value = 1; refresh(); });
    $("s").addEventListener("input", () => { state.s = Number($("s").value); refresh(); });
    $("watch").addEventListener("click", watch);
    for (const id of ["t-mixing", "t-allowed", "t-irreps"]) $(id).addEventListener("change", () => draw(state.view));
    for (const id of ["case", "s", "watch"]) $(id).disabled = false;
    await refresh();
    document.body.dataset.ready = "1";
  } catch (e) {
    $("banner").textContent = `Engine failed to load: ${e.message}`;
    document.body.dataset.error = e.message;
  }
})();
