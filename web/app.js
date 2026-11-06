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
  for (const lv of [...d.left, ...d.right]) {
    const g = levelRow(lv, d.left.includes(lv) ? "left" : "right", "frag", lv.id);
    g.append(el("title", {}, view.provenance.fragment));
    svg.append(g);
  }
  centreLevels.forEach((lv, k) => {
    const g = levelRow(lv, "centre", state.level === k ? "selected" : (atEnd ? "" : "moving"), atEnd ? `M${k}` : `S${k}`);
    g.append(el("title", {}, `${view.provenance.level_energy} ${atEnd ? view.provenance.level_label : ""}`.trim()));
    {
      const [x0, x1] = colX.centre, yy = y(lv.energy_eV);
      const hit = el("rect", { x: x0, y: yy - 7, width: x1 - x0, height: 14, class: "mo-hit", tabindex: 0,
                               [atEnd ? "data-mo" : "data-slider-level"]: k, "data-energy": lv.energy_eV });
      hit.addEventListener("click", () => select(k));
      hit.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") select(k); });
      g.append(hit);
      if (atEnd && lv.warnings && lv.warnings.length) g.append(el("text", { x: x0 - 6, y: yy + 4, "text-anchor": "end", class: "warn" }, "!"));
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
      const line = el("line", { x1: xa, y1: f.y, x2: xb, y2: t.y, class: ln.layer === "mixing" ? "mix" : "allowed",
                                "stroke-width": ln.layer === "mixing" ? (0.6 + 4 * ln.weight).toFixed(2) : 1 });
      line.append(el("title", {}, ln.layer === "mixing" ? `${view.provenance.mixing_line} Weight ${(100 * ln.weight).toFixed(0)} %.` : view.provenance.allowed_line));
      svg.append(line);
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
  const stage = sc.active_stage !== null ? sc.stages[sc.active_stage] : (sc.s > 0 ? sc.stages[Math.max(0, Math.round(sc.s * K) - 1)] : null);
  badgesFor($("badges"), stage ? stage.rules : [], view.rules);
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
  const wj = $("whyjump");
  wj.hidden = !(L && view.schedule.s >= 1);
  if (!L) { $("ledger-hint").textContent = "Tap a molecular orbital (at any point of the construction) to see every contact that sets its energy."; return; }
  const when = view.schedule.s >= 1 ? "" : ` at this point of the construction (s = ${view.schedule.s.toFixed(2)})`;
  $("ledger-hint").textContent = `Level ${L.level + 1}${when}: ledger total ${fmt(L.total_eV, 6)} eV; eigenvalue ${fmt(L.eigenvalue_eV, 6)} eV.`;
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
    tr.dataset.row = L.rows.indexOf(r);
    if (r.kind === "contact") { tr.addEventListener("mouseenter", () => hot(tr.dataset.row, true)); tr.addEventListener("mouseleave", () => hot(tr.dataset.row, false)); }
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

function about(view) {
  const ul = $("comparisons");
  ul.replaceChildren();
  for (const c of view.comparisons || []) {
    const li = document.createElement("li");
    li.className = c.agrees ? "agrees" : "disagrees";
    li.dataset.agrees = c.agrees ? "1" : "0";
    const v = document.createElement("span");
    v.className = "verdict";
    v.textContent = `${c.category}: ${c.agrees ? "agrees with experiment" : "disagrees with experiment"}. `;
    const t = document.createElement("span");
    t.textContent = `${c.model}. ${c.reference}.`;
    const src = document.createElement("span");
    src.className = "src";
    src.textContent = `Source: ${c.source}`;
    li.append(v, t, src);
    ul.append(li);
  }
  $("explanation").textContent = view.explanation;
  const box = $("notes");
  box.replaceChildren();
  for (const n of view.notes || []) {
    const h = document.createElement("h3");
    h.textContent = n.title;
    const p = document.createElement("p");
    p.textContent = n.text;
    box.append(h, p);
  }
}

function banner(view) {
  const tag = view.preset.provenance_tag;
  $("banner").textContent = `Model units (eV scale) · preset: ${tag} (${view.preset.id}) · not a spectral prediction · ` +
    `engine ${view.result.version} · CONVENTIONS v1.3, not frozen` + (view.result.warnings.length ? ` · warnings: ${view.result.warnings.join(", ")}` : "");
}

let latest = 0;
async function refresh() {
  const mine = ++latest;
  const view = await call("view", { case: state.case, s: state.s, level: state.level });
  if (mine !== latest) return state.view;  // a newer request superseded this one
  state.view = view;
  draw(view); narrate(view); ledger(view); bonds(view); banner(view); about(view);
  structure(view); matching(view); slotBadges(view);
  if ($("t-3d").checked) view3d.render(view);
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


// ---- Rule badges (L3): each step carries the rule it uses and the badge that rule earns. ----
function badgesFor(box, ids, rules) {
  box.replaceChildren();
  for (const id of ids || []) {
    const r = rules[id];
    if (!r) continue;
    const b = document.createElement("button");
    b.type = "button";
    b.className = "badge";
    b.dataset.rule = id;
    b.innerHTML = `${id} <span class="kind">${r.badge.split(":")[0].split(";")[0]}</span>`;
    b.title = `${r.rule}. ${r.badge}`;
    b.addEventListener("click", () => {
      const d = $("rule-detail");
      const same = !d.hidden && d.dataset.rule === id;
      d.hidden = same;
      d.dataset.rule = id;
      d.textContent = `${id}. ${r.rule}. In the model: ${r.model}. Badge: ${r.badge}. Switches that remove approximations arrive in Release 2.`;
    });
    box.append(b);
  }
}
function slotBadges(view) {
  for (const slot of document.querySelectorAll(".badge-slot")) badgesFor(slot, slot.dataset.rules.split(" "), view.rules);
}

// ---- Linked contacts: a ledger row and its structure line highlight together (L2). ----
function hot(row, on) {
  for (const e of document.querySelectorAll(`[data-row="${row}"]`)) e.classList.toggle("hot", on);
}

// ---- Structure view (L2): atoms, orbital lobes coloured by phase, contact lines. ----
function structure(view) {
  const svg = $("structure");
  svg.replaceChildren();
  const S = view.structure;
  const W = 600, H = 420, pad = 60;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  if (!S) {
    svg.append(el("text", { x: W / 2, y: H / 2, "text-anchor": "middle", class: "muted" }, "Select a level to see its orbital and contacts."));
    return;
  }
  const xs = S.atoms.map((a) => a.x), ys = S.atoms.map((a) => a.y);
  const span = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys), 1);
  const k = Math.min((W - 2 * pad), (H - 2 * pad)) / span;
  const cx = (Math.max(...xs) + Math.min(...xs)) / 2, cy = (Math.max(...ys) + Math.min(...ys)) / 2;
  const X = (x) => W / 2 + (x - cx) * k, Y = (y) => H / 2 - (y - cy) * k;
  const at = Object.fromEntries(S.atoms.map((a) => [a.index, a]));
  const pairCount = {};
  for (const c of S.contacts) {
    const key = c.atoms.slice().sort().join("-");
    const n = (pairCount[key] = (pairCount[key] || 0) + 1) - 1;
    const a = at[c.atoms[0]], b = at[c.atoms[1]];
    let x1 = X(a.x), y1 = Y(a.y), x2 = X(b.x), y2 = Y(b.y);
    const len = Math.hypot(x2 - x1, y2 - y1) || 1, nx = -(y2 - y1) / len, ny = (x2 - x1) / len;
    const off = (n % 2 ? -1 : 1) * Math.ceil(n / 2) * 5;
    const line = el("line", { x1: x1 + nx * off, y1: y1 + ny * off, x2: x2 + nx * off, y2: y2 + ny * off,
                              class: `contact ${c.style}`, "stroke-width": (0.8 + 2.2 * Math.min(1, Math.abs(c.energy_eV) / 4)).toFixed(2),
                              "data-row": c.ledger_row });
    line.append(el("title", {}, `${view.ledger.basis_labels[c.mu]}–${view.ledger.basis_labels[c.nu]}: ${c.phase}-phase, ${fmt(c.energy_eV)} eV (${c.style === "solid" ? "stabilising" : "destabilising"})`));
    line.addEventListener("mouseenter", () => hot(c.ledger_row, true));
    line.addEventListener("mouseleave", () => hot(c.ledger_row, false));
    svg.append(line);
  }
  const maxc = Math.max(...S.lobes.map((l) => Math.abs(l.coefficient)), 1e-9);
  const order = S.atoms.slice().sort((a, b) => a.depth - b.depth);
  for (const a of order) {
    for (const l of S.lobes.filter((l) => l.atom === a.index)) {
      const r = 8 + 30 * Math.abs(l.coefficient) / maxc;
      const x = X(a.x), y = Y(a.y);
      if (l.kind === "s") {
        svg.append(el("circle", { cx: x, cy: y, r, class: l.coefficient > 0 ? "lobe-pos" : "lobe-neg" }));
      } else {
        const [dx, dy] = l.dir2d, n = Math.hypot(dx, dy);
        if (n < 1e-6) {  // p orbital pointing at the viewer: draw it face-on
          svg.append(el("circle", { cx: x, cy: y, r: r * 0.6, class: l.coefficient > 0 ? "lobe-pos" : "lobe-neg" }));
          continue;
        }
        const ang = Math.atan2(-dy, dx) * 180 / Math.PI, len = r * Math.max(0.35, n);
        for (const sgn of [1, -1]) {
          const positive = (sgn > 0) === (l.coefficient > 0);
          svg.append(el("ellipse", { cx: x + sgn * len * Math.cos(ang * Math.PI / 180), cy: y + sgn * len * Math.sin(ang * Math.PI / 180),
                                     rx: len, ry: r * 0.42, transform: `rotate(${ang} ${x + sgn * len * Math.cos(ang * Math.PI / 180)} ${y + sgn * len * Math.sin(ang * Math.PI / 180)})`,
                                     class: positive ? "lobe-pos" : "lobe-neg" }));
        }
      }
    }
    svg.append(el("circle", { cx: X(a.x), cy: Y(a.y), r: 9, class: "atom" }));
    svg.append(el("text", { x: X(a.x), y: Y(a.y) + 4, "text-anchor": "middle" }, a.element));
  }
  svg.append(el("text", { x: 8, y: H - 8, class: "muted" }, "Red = positive phase, blue = negative. Solid line: stabilising contact; dashed: destabilising."));
}

// ---- Matching view (P6): coupling tallies between central orbitals and LGO components. ----
function matching(view) {
  const box = $("matching");
  box.replaceChildren();
  const m = view.diagram.matching;
  if (!m) { box.textContent = "No central atom: nothing to match."; return; }
  const t = document.createElement("table");
  t.innerHTML = `<thead><tr><th></th>${m.cols.map((c) => `<th title="${c.irrep}">${c.label}${c.irrep !== "—" ? `<br><small>${c.irrep}</small>` : ""}</th>`).join("")}</tr></thead>`;
  const tb = document.createElement("tbody");
  m.rows.forEach((r, i) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${r.label}${r.irrep !== "—" ? ` <small>${r.irrep}</small>` : ""}</td>` +
      m.tally_eV[i].map((v) => (v === 0 ? `<td class="zero">0</td>` : `<td>${fmt(v, 2)}</td>`)).join("");
    tb.append(tr);
  });
  t.append(tb);
  const note = document.createElement("p");
  note.className = "muted-text";
  note.textContent = m.rule;
  box.append(t, note);
}

// ---- Export and share (F5.3-F5.6). ----
function download(name, blob) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = name;
  document.body.append(a);
  a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
}
function standaloneSVG() {
  // Inline the computed styles so the exported figure matches the screen without the stylesheet.
  const src = $("diagram");
  const copy = src.cloneNode(true);
  copy.setAttribute("xmlns", NS);
  const all = src.querySelectorAll("*"), out = copy.querySelectorAll("*");
  all.forEach((e, i) => {
    const cs = getComputedStyle(e);
    for (const p of ["stroke", "stroke-width", "stroke-dasharray", "fill", "fill-opacity", "font-size", "font-family", "font-weight"]) out[i].style.setProperty(p, cs.getPropertyValue(p));
  });
  const bg = el("rect", { x: 0, y: 0, width: "100%", height: "100%", fill: getComputedStyle(src).backgroundColor });
  copy.insertBefore(bg, copy.firstChild);
  copy.append(el("text", { x: 10, y: copy.viewBox.baseVal.height - 8, style: "font-size:10px;fill:#888" },
                 `${state.case} · ${$("banner").textContent}`));
  return new XMLSerializer().serializeToString(copy);
}
const fileBase = () => `lgo-${state.case.replace(/[^A-Za-z0-9+]+/g, "_")}`;
async function exportSVG() { download(`${fileBase()}.svg`, new Blob([standaloneSVG()], { type: "image/svg+xml" })); status("SVG downloaded."); }
async function exportPNG() {
  const svgText = standaloneSVG();
  const img = new Image();
  const url = URL.createObjectURL(new Blob([svgText], { type: "image/svg+xml" }));
  await new Promise((r, j) => { img.onload = r; img.onerror = j; img.src = url; });
  const vb = $("diagram").viewBox.baseVal, scale = 2;
  const c = document.createElement("canvas");
  c.width = vb.width * scale; c.height = vb.height * scale;
  c.getContext("2d").drawImage(img, 0, 0, c.width, c.height);
  URL.revokeObjectURL(url);
  c.toBlob((b) => { download(`${fileBase()}.png`, b); status("PNG downloaded."); }, "image/png");
}
async function exportJSON() {
  const data = await call("export", { case: state.case, s: state.s, level: state.level });
  download(`${fileBase()}.json`, new Blob([JSON.stringify(data, null, 2)], { type: "application/json" }));
  status("JSON downloaded: config, full result, engine version and CONVENTIONS hash.");
  return data;
}
async function importJSON(text) {
  let data;
  try { data = JSON.parse(text); } catch (e) { status(`Import rejected: not valid JSON (${e.message}).`); return; }
  try {
    const ok = await call("import", data);
    Object.assign(state, { case: ok.config.case, s: ok.config.s, level: ok.config.level });
    $("case").value = state.case; $("s").value = state.s;
    await refresh();
    status(["Imported.", ...ok.notes].join(" "));
  } catch (e) {
    status(`Import rejected: ${String(e.message).split("ValueError: ").pop().trim()}`);
  }
}
const b64url = (bytes) => btoa(String.fromCharCode(...bytes)).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
const unb64url = (s) => Uint8Array.from(atob(s.replace(/-/g, "+").replace(/_/g, "/")), (c) => c.charCodeAt(0));
async function pack(obj) {
  const stream = new Blob([JSON.stringify(obj)]).stream().pipeThrough(new CompressionStream("deflate-raw"));
  return b64url(new Uint8Array(await new Response(stream).arrayBuffer()));
}
async function unpack(s) {
  const stream = new Blob([unb64url(s)]).stream().pipeThrough(new DecompressionStream("deflate-raw"));
  return JSON.parse(await new Response(stream).text());
}
async function shareLink() {
  const cfg = { v: 1, case: state.case, s: state.s, level: state.level,
                show: { mixing: $("t-mixing").checked, allowed: $("t-allowed").checked, irreps: $("t-irreps").checked } };
  const url = `${location.origin}${location.pathname}#c=${await pack(cfg)}`;
  history.replaceState(null, "", url);
  try { await navigator.clipboard.writeText(url); status("Link copied. Opening it reproduces this view; nothing is stored on a server."); }
  catch { status(`Link: ${url}`); }
  return url;
}
async function fromLink() {
  const m = location.hash.match(/#c=([A-Za-z0-9_-]+)/);
  if (!m) return false;
  try {
    const cfg = await unpack(m[1]);
    Object.assign(state, { case: cfg.case, s: cfg.s, level: cfg.level });
    $("t-mixing").checked = cfg.show.mixing; $("t-allowed").checked = cfg.show.allowed; $("t-irreps").checked = cfg.show.irreps;
    return true;
  } catch { status("The link could not be read; showing the default view."); return false; }
}
const status = (t) => { $("export-status").textContent = t; };

// ---- 3D viewer (F3.7): WebGL2, atoms and the selected orbital's lobes; drag to rotate. ----
const view3d = (() => {
  let gl, prog, sphere, rot = [0.5, -0.4], lastView = null, buf = {};
  const VS = `#version 300 es
  in vec3 p; uniform mat4 M; uniform mat4 P; uniform vec3 c; uniform vec3 s; out vec3 n; out vec3 col;
  void main(){ vec4 w = M * vec4(p * s + c, 1.0); n = normalize(mat3(M) * (p / s)); col = vec3(0); gl_Position = P * w; }`;
  const FS = `#version 300 es
  precision mediump float; in vec3 n; uniform vec4 k; out vec4 o;
  void main(){ float d = max(dot(normalize(n), normalize(vec3(0.4,0.6,1.0))), 0.0); o = vec4(k.rgb * (0.35 + 0.65*d), k.a); }`;
  function icosphere() {
    const pts = [], idx = [];
    const st = 16, sl = 24;
    for (let i = 0; i <= st; i++) for (let j = 0; j <= sl; j++) {
      const t = Math.PI * i / st, f = 2 * Math.PI * j / sl;
      pts.push(Math.sin(t) * Math.cos(f), Math.cos(t), Math.sin(t) * Math.sin(f));
    }
    for (let i = 0; i < st; i++) for (let j = 0; j < sl; j++) {
      const a = i * (sl + 1) + j, b = a + sl + 1;
      idx.push(a, b, a + 1, b, b + 1, a + 1);
    }
    return { pts: new Float32Array(pts), idx: new Uint16Array(idx) };
  }
  function init(canvas) {
    gl = canvas.getContext("webgl2", { antialias: true, preserveDrawingBuffer: true });
    if (!gl) return false;
    const sh = (t, src) => { const s = gl.createShader(t); gl.shaderSource(s, src); gl.compileShader(s); return s; };
    prog = gl.createProgram();
    gl.attachShader(prog, sh(gl.VERTEX_SHADER, VS)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, FS)); gl.linkProgram(prog);
    sphere = icosphere();
    const vao = gl.createVertexArray(); gl.bindVertexArray(vao);
    const vb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, vb); gl.bufferData(gl.ARRAY_BUFFER, sphere.pts, gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(prog, "p"); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
    const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, sphere.idx, gl.STATIC_DRAW);
    let drag = null;
    canvas.addEventListener("pointerdown", (e) => { drag = [e.clientX, e.clientY]; canvas.setPointerCapture(e.pointerId); });
    canvas.addEventListener("pointermove", (e) => {
      if (!drag) return;
      rot[0] += (e.clientX - drag[0]) * 0.01; rot[1] += (e.clientY - drag[1]) * 0.01; drag = [e.clientX, e.clientY];
      if (lastView) render(lastView);
    });
    canvas.addEventListener("pointerup", () => { drag = null; });
    return true;
  }
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const rgb = (hex) => { const h = hex.replace("#", ""); return [0, 2, 4].map((i) => parseInt(h.slice(i, i + 2), 16) / 255); };
  function matMul(a, b) { const o = new Float32Array(16); for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) { let s = 0; for (let k = 0; k < 4; k++) s += a[k * 4 + j] * b[i * 4 + k]; o[i * 4 + j] = s; } return o; }
  function render(view) {
    lastView = view;
    const canvas = $("view3d");
    if (!gl && !init(canvas)) { canvas.replaceWith(Object.assign(document.createElement("p"), { textContent: "3D view needs WebGL2." })); return; }
    const S = view.structure;
    gl.viewport(0, 0, canvas.width, canvas.height);
    const bg = rgb(css("--panel").length === 7 ? css("--panel") : "#ffffff");
    gl.clearColor(bg[0], bg[1], bg[2], 1); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST); gl.useProgram(prog);
    if (!S) return;
    const ext = Math.max(...S.atoms.map((a) => Math.hypot(...a.xyz)), 1) + 1.2;
    const [ay, ax] = rot, cy = Math.cos(ay), sy = Math.sin(ay), cx = Math.cos(ax), sx = Math.sin(ax);
    const Ry = new Float32Array([cy, 0, -sy, 0, 0, 1, 0, 0, sy, 0, cy, 0, 0, 0, 0, 1]);
    const Rx = new Float32Array([1, 0, 0, 0, 0, cx, sx, 0, 0, -sx, cx, 0, 0, 0, 0, 1]);
    // Engine coordinates have z up; map z to screen y.
    const Up = new Float32Array([1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1]);
    const M = matMul(Rx, matMul(Ry, Up));
    const a = canvas.width / canvas.height, sc = 1 / ext;
    const P = new Float32Array([sc / a, 0, 0, 0, 0, sc, 0, 0, 0, 0, -sc * 0.5, 0, 0, 0, 0, 1]);
    gl.uniformMatrix4fv(gl.getUniformLocation(prog, "M"), false, M);
    gl.uniformMatrix4fv(gl.getUniformLocation(prog, "P"), false, P);
    const draw = (c, s, col) => {
      gl.uniform3fv(gl.getUniformLocation(prog, "c"), c); gl.uniform3fv(gl.getUniformLocation(prog, "s"), s);
      gl.uniform4fv(gl.getUniformLocation(prog, "k"), col);
      gl.drawElements(gl.TRIANGLES, sphere.idx.length, gl.UNSIGNED_SHORT, 0);
    };
    const pos = rgb(css("--pos")), neg = rgb(css("--neg")), atom = rgb(css("--muted"));
    const maxc = Math.max(...S.lobes.map((l) => Math.abs(l.coefficient)), 1e-9);
    for (const at of S.atoms) draw(at.xyz, [0.18, 0.18, 0.18], [...atom, 1]);
    for (const l of S.lobes) {
      const at = S.atoms[l.atom].xyz, r = 0.18 + 0.55 * Math.abs(l.coefficient) / maxc;
      if (l.kind === "s") { draw(at, [r, r, r], [...(l.coefficient > 0 ? pos : neg), 1]); continue; }
      const d = l.dir3d;
      for (const sgn of [1, -1]) {
        const positive = (sgn > 0) === (l.coefficient > 0);
        const c = at.map((v, i) => v + sgn * d[i] * r);
        // Ellipsoid elongated along d: scale is axis-aligned, so use a sphere of radius r along the
        // dominant axis. Adequate for the Cartesian and axial directions used in Release 1.
        const s = [0, 1, 2].map((i) => (Math.abs(d[i]) > 0.5 ? r : r * 0.55));
        draw(c, s, [...(positive ? pos : neg), 1]);
      }
    }
  }
  return { render };
})();

(async () => {
  try {
    const info = await call("info");
    state.info = info;
    $("footer").textContent = `Pyodide · Python ${info.python} · NumPy ${info.numpy} · first load ${(info.load_ms / 1000).toFixed(1)} s`;
    const cases = await call("cases");
    const sel = $("case");
    for (const c of cases) sel.append(new Option(`${c.name} (${c.group})`, c.name));
    state.case = cases.find((c) => c.name === "SF6 model").name;
    if (await fromLink()) $("s").value = state.s;
    sel.value = state.case;
    sel.addEventListener("change", () => { state.case = sel.value; state.level = null; state.s = 1; $("s").value = 1; refresh(); });
    $("s").addEventListener("input", () => { state.s = Number($("s").value); state.level = null; refresh(); });
    $("whyjump").addEventListener("click", () => {
      const o = state.view.origins[state.level];
      state.s = o.s; state.level = o.level_index_there; $("s").value = state.s;
      refresh().then(() => status(`Jumped to where this level appeared: ${o.label}`));
    });
    $("ex-svg").addEventListener("click", exportSVG);
    $("ex-png").addEventListener("click", exportPNG);
    $("ex-json").addEventListener("click", exportJSON);
    $("ex-link").addEventListener("click", shareLink);
    $("im-json").addEventListener("change", async (e) => { const f = e.target.files[0]; if (f) await importJSON(await f.text()); e.target.value = ""; });
    $("t-3d").addEventListener("change", () => { $("view3d").hidden = !$("t-3d").checked; $("structure").hidden = $("t-3d").checked; if ($("t-3d").checked) view3d.render(state.view); });
    window.lgo = { exportJSON, importJSON, shareLink, state };
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
