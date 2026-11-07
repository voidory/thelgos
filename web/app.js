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
const state = { case: null, s: 1, level: null, fragment: null, view: null, info: null, watching: false,
  electrons: null, approx: null, dial: 0, lp: 0 };

function el(tag, attrs = {}, text) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (text !== undefined) e.textContent = text;
  return e;
}

const fmt = (x, d = 3) => (x >= 0 ? "+" : "−") + Math.abs(x).toFixed(d);
// ---- Orbital labels in textbook form, set in STIX Two with real subscripts. ----
// Mulliken symbols are written in lower case for orbitals (1a₁g, 2t₁u); linear molecules use σ/π with
// g/u subscripts. The star marks a level the engine classifies as antibonding.
const GREEK = { sigma: "σ", pi: "π" };
function labelParts(lv, { number = true } = {}) {
  // Returns [[text, "base"|"sub"|"sup"], ...]
  const out = [];
  if (number && lv.number) out.push([String(lv.number), "base"]);
  if (lv.linear_label) {
    const [base, par] = lv.linear_label.split("_");
    out.push([GREEK[base] || base, "base"]);
    if (par) out.push([par, "sub"]);
  } else if (lv.irrep && lv.irrep !== "—") {
    const m = lv.irrep.match(/^([A-Z])(\d?)([gu]?)('*)$/);
    if (m) {
      out.push([m[1].toLowerCase(), "base"]);
      if (m[2] || m[3]) out.push([m[2] + m[3], "sub"]);
      if (m[4]) out.push([m[4].length === 1 ? "′" : "″", "base"]);
    } else out.push([lv.irrep, "base"]);
  }
  if (lv.character === "antibonding") out.push(["*", "sup"]);
  return out;
}
const irrepText = (lv) => labelParts(lv).map((p) => p[0]).join("");
function richText(parts, attrs) {
  // Sub- and superscripts by explicit dy offsets (baseline-shift is unreliable across browsers),
  // so a star after a subscript (σg*) sits up as a superscript.
  const t = el("text", attrs);
  let shift = 0;
  for (const [txt, kind] of parts) {
    const target = kind === "sub" ? 5 : kind === "sup" ? -8 : 0;
    const a = { dy: target - shift };
    if (kind === "sub" || kind === "sup") a["font-size"] = "68%";
    if (kind === "tag") { a.class = "frontier"; a.dx = 8; }
    if (kind === "warn") { a.class = "warn-mark"; a.dx = 4; }
    t.append(el("tspan", a, txt));
    shift = target;
  }
  return t;
}
function orbitalNameParts(label) {
  // "O 2p" -> O, 2p ; "S 3px" -> S, 3p, sub x
  const m = label.match(/^([A-Z][a-z]?) (\d)([spdf])([xyz]?)$/);
  if (!m) return [[label, "base"]];
  const parts = [[`${m[1]} `, "base"], [`${m[2]}${m[3]}`, "base"]];
  if (m[4]) parts.push([m[4], "sub"]);
  return parts;
}

// Electron pair drawn as half-headed arrows (harpoons), the way they are drawn by hand.
function harpoon(x, y, up) {
  const h = 13, top = y - h / 2, bot = y + h / 2;
  const d = up ? `M${x},${bot} L${x},${top} L${x - 3.5},${top + 6}` : `M${x},${top} L${x},${bot} L${x + 3.5},${bot - 6}`;
  return el("path", { d, class: "electron" });
}

// Plain-text name of a level for screen readers ("1t1u*, −27.80 eV, 6 electrons, HOMO").
function levelName(lv, k, atEnd, view) {
  const lab = atEnd ? labelParts(lv).map(([t]) => t).join("") : `Level ${k + 1}`;
  const bits = [lab || `Level ${k + 1}`, `${fmt(lv.energy_eV, 2)} eV`];
  if (atEnd && lv.electrons !== undefined) bits.push(`${lv.electrons} electron${lv.electrons === 1 ? "" : "s"}`);
  if (atEnd && k === view.result.frontier.HOMO) bits.push("HOMO");
  if (atEnd && k === view.result.frontier.LUMO) bits.push("LUMO");
  if (lv.multiplicity > 1) bits.push(`${lv.multiplicity}-fold degenerate`);
  return bits.join(", ");
}

// Text alternative to the diagram: every level in a table, read by screen readers, hidden on screen.
function levelTable(view) {
  const box = $("levels-sr");
  if (!box) return;
  const atEnd = view.schedule.s >= 1;
  const list = atEnd ? view.diagram.centre : view.diagram.slider_levels;
  const t = document.createElement("table");
  t.innerHTML = `<caption>Molecular-orbital levels of ${view.case}${atEnd ? "" : ` at construction point s = ${view.schedule.s.toFixed(2)}`}, highest first</caption>
    <thead><tr><th scope="col">Level</th><th scope="col">Energy / eV</th><th scope="col">Degeneracy</th><th scope="col">Electrons</th></tr></thead>`;
  const tb = document.createElement("tbody");
  list.map((lv, k) => [lv, k]).sort((a, b) => b[0].energy_eV - a[0].energy_eV).forEach(([lv, k]) => {
    const tr = document.createElement("tr");
    for (const v of [levelName(lv, k, atEnd, view).split(", ")[0] + (atEnd && k === view.result.frontier.HOMO ? " (HOMO)" : atEnd && k === view.result.frontier.LUMO ? " (LUMO)" : ""),
                     fmt(lv.energy_eV, 3), String(lv.multiplicity), lv.electrons === undefined ? "–" : String(lv.electrons)]) {
      const td = document.createElement("td"); td.textContent = v; tr.append(td);
    }
    tb.append(tr);
  });
  t.append(tb);
  box.replaceChildren(t);
}

// Column titles are set in capitals, except a leading element symbol ("Si (atom 1)" -> "Si (ATOM 1)", not "SI").
function capsTitle(t) {
  const m = t.match(/^([A-Z][a-z]?)(\s.*)$/);
  return m ? m[1] + m[2].toUpperCase() : t.toUpperCase();
}

function draw(view) {
  const svg = $("diagram");
  // Keep keyboard focus on the same level across redraws (selecting a level redraws the diagram).
  const fe = document.activeElement, refocus = fe && svg.contains(fe) && fe.classList.contains("mo-hit") ? fe.dataset.mo ?? fe.dataset.sliderLevel ?? fe.dataset.frag : null;
  svg.replaceChildren();
  const d = view.diagram;
  const atEnd = view.schedule.s >= 1;
  const centreLevels = atEnd ? d.centre : d.slider_levels;
  const three = d.columns.length === 3;
  // Compact geometry on narrow screens: a narrower canvas so the text is not scaled down to nothing.
  const compact = (svg.parentElement.clientWidth || 640) < 560;
  svg.classList.toggle("compact", compact);
  const W = compact ? 420 : 640, H = compact ? 520 : 600, top = compact ? 44 : 52, bottom = compact ? 28 : 34;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const all = [...d.left, ...d.right, ...d.centre, ...d.slider_levels].map((x) => x.energy_eV);
  let lo = Math.min(...all), hi = Math.max(...all);
  const pad = (hi - lo) * 0.05 || 1;
  lo -= pad; hi += pad;
  // F3.6 symmetry-only view: heights show order only (energies collapsed to evenly spaced ranks), since
  // symmetry fixes labels and connections but not energies. Labels and lines are unchanged.
  const symOnly = $("t-symonly") && $("t-symonly").checked;
  const ranks = [...all].sort((a, b) => a - b).filter((e, i, a) => i === 0 || e - a[i - 1] > 1e-6);
  const rankOf = (e) => ranks.findIndex((r) => Math.abs(r - e) <= 1e-6);
  const y = symOnly
    ? (e) => H - bottom - 14 - (rankOf(e) / Math.max(ranks.length - 1, 1)) * (H - top - bottom - 28)
    : (e) => top + (hi - e) / (hi - lo) * (H - top - bottom);
  svg.classList.toggle("symonly", !!symOnly);
  // Column centres: short bars sit at these; labels go in the gutters beside them.
  const cx = compact ? (three ? { left: 120, centre: 232, right: 340 } : { right: 150, centre: 295 })
                     : (three ? { left: 165, centre: 335, right: 505 } : { right: 200, centre: 420 });
  const BAR = compact ? 18 : 30, GAP = compact ? 4 : 7;

  // Energy axis with quiet ticks.
  const step = (hi - lo) > 60 ? 10 : 5;
  const AX = compact ? 32 : 52;
  svg.append(el("line", { x1: AX, x2: AX, y1: top - 6, y2: H - bottom + 4, class: "axis-line" }));
  for (let e = Math.ceil(lo / step) * step; e <= hi && !symOnly; e += step) {
    svg.append(el("line", { x1: AX - 4, x2: AX, y1: y(e), y2: y(e), class: "axis-line" }));
    svg.append(el("text", { x: AX - 7, y: y(e) + 4, "text-anchor": "end", class: "axis-num" }, e === 0 ? "0" : String(e).replace("-", "−")));
  }
  const axisTitle = symOnly ? "order only (symmetry)" : compact ? "E / eV" : "E / eV (model)";
  if (compact) svg.append(el("text", { x: 4, y: top - 14, class: "axis-title" }, axisTitle));
  else svg.append(el("text", { x: 12, y: (top + H - bottom) / 2, transform: `rotate(-90 12 ${(top + H - bottom) / 2})`, "text-anchor": "middle", class: "axis-title" }, axisTitle));
  const short = { "Molecular orbitals": "MOs", "Ligand group orbitals": "LGOs" };
  for (const c of d.columns) svg.append(el("text", { x: cx[c], y: 24, "text-anchor": "middle", class: "col-title" },
                                           compact ? (short[d.column_titles[c]] || d.column_titles[c].replace(" orbitals", "")) : capsTitle(d.column_titles[c])));

  const pos = {};
  const labels = { left: [], right: [], centre: [] };
  // Lay out one column: rows of levels at the same energy share a centred group of short bars.
  function layout(col, list, isMO) {
    // Rows: levels closer than 7 px on screen share a row and sit side by side (each still at its own
    // height), so near-degenerate levels such as water's 2b2*/3a1* stay distinguishable. Groups that
    // belong to different levels are separated by a wider gap.
    const rows = [];
    list.forEach((lv, k) => {
      const r = rows.find((r) => Math.abs(y(r.e) - y(lv.energy_eV)) < 7);
      if (r) { r.items.push([lv, k]); r.m += lv.multiplicity; } else rows.push({ e: lv.energy_eV, items: [[lv, k]], m: lv.multiplicity });
    });
    const GROUP = 16;
    for (const r of rows) {
      const width = r.m * BAR + (r.m - r.items.length) * GAP + (r.items.length - 1) * GROUP;
      let x = cx[col] - width / 2;
      for (const [lv, k] of r.items) {
        const yy = y(lv.energy_eV);
        const g = el("g");
        const x0 = x;
        for (let i = 0; i < lv.multiplicity; i++) {
          const a = x + i * (BAR + GAP);
          const cls = isMO ? (state.level === k ? "level selected" : atEnd ? "level" : "level moving") : (state.fragment === lv.id ? "level frag selected" : "level frag");
          g.append(el("line", { x1: a, x2: a + BAR, y1: yy, y2: yy, class: cls }));
          if (isMO && atEnd && lv.electrons !== undefined) {
            const n = (i < lv.electrons ? 1 : 0) + (i < lv.electrons - lv.multiplicity ? 1 : 0);
            if (n >= 1) g.append(harpoon(a + BAR / 2 - (n === 2 ? 5 : 0), yy, true));
            if (n === 2) g.append(harpoon(a + BAR / 2 + 5, yy, false));
          }
        }
        x += lv.multiplicity * BAR + (lv.multiplicity - 1) * GAP + GROUP;
        const x1 = x - GROUP;
        const key = isMO ? (atEnd ? `M${k}` : `S${k}`) : lv.id;
        pos[key] = { x0, x1, y: yy };
        if (isMO) {
          g.append(el("title", {}, `${view.provenance.level_energy} ${atEnd ? view.provenance.level_label : ""} E = ${fmt(lv.energy_eV, 3)} eV`.trim()));
          const hit = el("rect", { x: x0 - 4, y: yy - 9, width: x1 - x0 + 8, height: 18, class: "mo-hit", tabindex: 0,
                                   role: "button", "aria-pressed": state.level === k, "aria-label": levelName(lv, k, atEnd, view),
                                   [atEnd ? "data-mo" : "data-slider-level"]: k, "data-energy": lv.energy_eV, "data-multiplicity": lv.multiplicity });
          hit.addEventListener("click", () => select(k));
          hit.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") select(k); });
          g.append(hit);
        } else {
          g.append(el("title", {}, view.provenance.fragment));
          // Fragment orbitals (LGOs, central-atom or atom orbitals) are selectable too: ledger, structure and 3D follow.
          const fname = `${lv.label}${lv.irrep && lv.irrep !== "—" ? ` (${lv.irrep})` : ""}, ${fmt(lv.energy_eV, 2)} eV${lv.multiplicity > 1 ? `, ${lv.multiplicity}-fold degenerate` : ""}`;
          const fh = el("rect", { x: x0 - 4, y: yy - 9, width: x1 - x0 + 8, height: 18, class: "mo-hit frag-hit", tabindex: 0, role: "button",
                                  "aria-pressed": state.fragment === lv.id, "aria-label": fname, "data-frag": lv.id });
          fh.addEventListener("click", () => selectFragment(lv.id));
          fh.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") selectFragment(lv.id); });
          g.append(fh);
        }
        svg.append(g);
      }
      // One label per row; the gutter is to the right of MO and right-hand bars, left of left-hand bars.
      const rightmost = cx[col] + width / 2, leftmost = cx[col] - width / 2;
      const yy = y(r.e);
      const parts = [];
      const bases = r.items.map(([lv]) => (lv.label || "").replace(/[xyz]$/, ""));
      const oneBase = !isMO && !(col === "right" && d.column_titles.right === "Ligand group orbitals") && bases.every((b) => b === bases[0]) && bases[0];
      if (oneBase) parts.push(...orbitalNameParts(bases[0]));
      else r.items.forEach(([lv], i) => {
        if (i) parts.push([", ", "base"]);
        if (isMO) { if (atEnd && $("t-irreps").checked) parts.push(...labelParts(lv)); }
        else if (col === "right" && d.column_titles.right === "Ligand group orbitals") parts.push(...labelParts(lv, { number: false }));
        else parts.push(...orbitalNameParts(lv.label));
      });
      const homo = isMO && atEnd && r.items.some(([, k]) => k === view.result.frontier.HOMO);
      const lumo = isMO && atEnd && r.items.some(([, k]) => k === view.result.frontier.LUMO);
      const warn = isMO && atEnd && r.items.some(([lv]) => lv.warnings && lv.warnings.length);
      if (homo || lumo) parts.push([homo && lumo ? "HOMO = LUMO" : homo ? "HOMO" : "LUMO", "tag"]);
      if (warn) parts.push(["◆", "warn"]);
      if (parts.length) labels[col].push({ y: yy, x: col === "left" ? leftmost - 10 : rightmost + 10, anchor: col === "left" ? "end" : "start", parts, ox: col === "left" ? leftmost : rightmost });
    }
  }
  layout("left", d.left, false);
  layout("right", d.right, false);
  layout("centre", centreLevels, true);

  // Correlation lines, under the labels.
  if (atEnd) {
    const g = el("g");
    for (const ln of d.lines) {
      const on = ln.layer === "mixing" ? $("t-mixing").checked : $("t-allowed").checked;
      if (!on || !pos[ln.from] || !pos[ln.to]) continue;
      const f = pos[ln.from], t = pos[ln.to];
      const fromRight = f.x0 > t.x0;
      const line = el("line", { x1: fromRight ? f.x0 : f.x1, y1: f.y, x2: fromRight ? t.x1 : t.x0, y2: t.y,
                                class: ln.layer === "mixing" ? "mix" : "allowed",
                                "stroke-width": ln.layer === "mixing" ? (0.5 + 2.2 * ln.weight).toFixed(2) : 0.9 });
      line.append(el("title", {}, ln.layer === "mixing" ? `${view.provenance.mixing_line} Weight ${(100 * ln.weight).toFixed(0)} %.` : view.provenance.allowed_line));
      g.append(line);
    }
    svg.insertBefore(g, svg.firstChild.nextSibling);
  }

  // Labels: keep at least 17 px apart in each gutter; a short leader line shows where a label was moved from.
  for (const col of Object.keys(labels)) {
    const L = labels[col].sort((a, b) => a.y - b.y);
    for (let it = 0; it < 30; it++) {
      let moved = false;
      for (let i = 1; i < L.length; i++) {
        const need = 17 - (L[i].ly ?? L[i].y) + (L[i - 1].ly ?? L[i - 1].y);
        if (need > 0) { L[i - 1].ly = (L[i - 1].ly ?? L[i - 1].y) - need / 2; L[i].ly = (L[i].ly ?? L[i].y) + need / 2; moved = true; }
      }
      if (!moved) break;
    }
    for (const lb of L) {
      const ly = lb.ly ?? lb.y;
      if (Math.abs(ly - lb.y) > 1.5) {
        const xa = lb.anchor === "start" ? lb.ox + 2 : lb.ox - 2, xb = lb.anchor === "start" ? lb.x - 2 : lb.x + 2;
        svg.append(el("polyline", { points: `${xa},${lb.y} ${(xa + xb) / 2},${lb.y} ${xb},${ly}`, class: "leader" }));
      }
      svg.append(richText(lb.parts, { x: lb.x, y: ly + 5, "text-anchor": lb.anchor, class: "orb-label" }));
    }
  }

  // Sketch card for the selected level: the orbital drawn as in the lecture, placed in the emptiest
  // spot near its level (fewest crossed correlation lines, no overlap with bars or labels).
  // Smaller in the compact geometry; left out when no spot is free of bars and labels.
  const pk = state.fragment ?? (atEnd ? `M${state.level}` : `S${state.level}`);
  if (view.structure && (state.level !== null || state.fragment) && pos[pk]) sketchCard(svg, view.structure, pos[pk], { W, H, top, bottom, AX, size: compact ? 84 : 108 });
  // Never clip: if labels reach past the canvas (long labels in the compact geometry), widen the viewBox to fit.
  const bb = svg.getBBox();
  if (bb.width > 0) {
    const x0 = bb.x < 0 ? bb.x - 4 : 0, x1 = bb.x + bb.width > W ? bb.x + bb.width + 4 : W;
    if (x0 < 0 || x1 > W) svg.setAttribute("viewBox", `${x0} 0 ${x1 - x0} ${H}`);
  }
  if (refocus != null) { const h = svg.querySelector(`[data-mo="${refocus}"], [data-slider-level="${refocus}"], [data-frag="${refocus}"]`); if (h) h.focus({ preventScroll: true }); }
}

function sketchCard(svg, S, p, { W, H, top, bottom, AX, size }) {
  const m = 6;
  const boxes = [], segs = [];
  for (const t of svg.querySelectorAll("text")) { const b = t.getBBox(); boxes.push([b.x, b.y, b.x + b.width, b.y + b.height]); }
  for (const l of svg.querySelectorAll("line.level")) boxes.push([+l.getAttribute("x1"), +l.getAttribute("y1") - 8, +l.getAttribute("x2"), +l.getAttribute("y1") + 8]);
  for (const l of svg.querySelectorAll("line.mix, line.allowed")) segs.push(["x1", "y1", "x2", "y2"].map((k) => +l.getAttribute(k)));
  const hitsSeg = (x0, y0, x1, y1, [ax, ay, bx, by]) => {
    // Liang–Barsky clip of segment against the rectangle.
    let t0 = 0, t1 = 1;
    const dx = bx - ax, dy = by - ay;
    for (const [pp, q] of [[-dx, ax - x0], [dx, x1 - ax], [-dy, ay - y0], [dy, y1 - ay]]) {
      if (pp === 0) { if (q < 0) return false; continue; }
      const r = q / pp;
      if (pp < 0) { if (r > t1) return false; if (r > t0) t0 = r; } else { if (r < t0) return false; if (r < t1) t1 = r; }
    }
    return true;
  };
  let best = null;
  for (let bx = AX + 10; bx <= W - size - 4; bx += 8) {
    for (let by = top; by <= H - bottom - size; by += 8) {
      const x0 = bx - m, y0 = by - m, x1 = bx + size + m, y1 = by + size + m;
      let cost = 0;
      for (const [a, b, c, d] of boxes) if (a < x1 && c > x0 && b < y1 && d > y0) cost += 100;
      for (const sg of segs) if (hitsSeg(x0, y0, x1, y1, sg)) cost += 4;
      const ex = bx + size < p.x0 ? p.x0 - (bx + size) : bx > p.x1 ? bx - p.x1 : 0;
      const ey = by > p.y ? by - p.y : by + size < p.y ? p.y - by - size : 0;
      cost += Math.hypot(ex, ey) / 25;
      if (!best || cost < best.cost) best = { cost, bx, by };
    }
  }
  if (best.cost >= 100) return;
  const { bx, by } = best;
  const xs = S.atoms.map((a) => a.x), ys = S.atoms.map((a) => a.y);
  const span = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys), 1);
  const k = (size - 40) / span, ccx = (Math.max(...xs) + Math.min(...xs)) / 2, ccy = (Math.max(...ys) + Math.min(...ys)) / 2;
  const card = el("g", { class: "sketch", "data-sketch": state.level });
  // Leader from the nearest point of the card to the nearest point just off the bar (not through its label).
  const cxm = bx + size / 2, cym = by + size / 2;
  const lx = Math.min(Math.max(cxm, p.x0), p.x1);
  const ly = p.y + (cym > p.y ? 5 : -5), ex = lx;
  const cxp = Math.min(Math.max(ex, bx), bx + size), cyp = Math.min(Math.max(ly, by), by + size);
  card.append(el("line", { x1: cxp, y1: cyp, x2: ex, y2: ly, class: "sketch-leader" }));
  card.append(el("rect", { x: bx, y: by, width: size, height: size, rx: 8, class: "sketch-card" }));
  // Lobes are clipped to the card so large ones (SF6's ligand p lobes) do not spill over the diagram.
  const clip = el("clipPath", { id: "sketch-clip" });
  clip.append(el("rect", { x: bx, y: by, width: size, height: size, rx: 8 }));
  card.append(clip);
  const inner = el("g", { "clip-path": "url(#sketch-clip)" });
  drawOrbital(inner, S, (x) => bx + size / 2 + (x - ccx) * k, (y) => by + size / 2 - (y - ccy) * k, 0.42 * size / 108, false);
  card.append(inner);
  svg.append(card);
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
  wj.hidden = !(L && L.level !== null && view.schedule.s >= 1);
  if (!L) { $("ledger-hint").textContent = "Tap a molecular orbital (at any point of the construction) to see every contact that sets its energy."; return; }
  const when = view.schedule.s >= 1 ? "" : ` at this point of the construction (s = ${view.schedule.s.toFixed(2)})`;
  $("ledger-hint").textContent = L.fragment
    ? `${L.fragment_label}${L.multiplicity > 1 ? ` (one of ${L.multiplicity} degenerate components)` : ""}: ledger total ${fmt(L.total_eV, 6)} eV; orbital energy ${fmt(L.eigenvalue_eV, 6)} eV. Only contacts within its own set appear: it has not met the other side yet.`
    : `Level ${L.level + 1}${when}: ledger total ${fmt(L.total_eV, 6)} eV; eigenvalue ${fmt(L.eigenvalue_eV, 6)} eV.`;
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

// Photoelectron spectrum (Koopmans' picture): bands at -energy, height = electrons. Measured markers where curated.
function pes(view) {
  const svg = $("pes"); if (!svg) return;
  svg.replaceChildren();
  const P = view.pes, W = 360, H = 190, L = 34, R = 12, T = 18, B = 34;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const ref = P.reference;
  const xs = P.bands.map((b) => b.ie_eV).concat(ref ? [ref.first_band_adiabatic_eV, ref.second_band_near_eV, ref.he_i_eV] : []);
  const lo = Math.floor(Math.min(...xs) / 5) * 5 - 2, hi = Math.ceil(Math.max(...xs) / 5) * 5 + 2;
  const x = (e) => L + (e - lo) / (hi - lo) * (W - L - R);
  const maxE = Math.max(...P.bands.map((b) => b.electrons), 1);
  const y = (n) => H - B - n / maxE * (H - B - T - 14);
  svg.append(el("line", { x1: L, x2: W - R, y1: H - B, y2: H - B, class: "axis-line" }));
  for (let e = Math.ceil(lo / 5) * 5; e <= hi; e += 5) {
    svg.append(el("line", { x1: x(e), x2: x(e), y1: H - B, y2: H - B + 4, class: "axis-line" }));
    svg.append(el("text", { x: x(e), y: H - B + 16, "text-anchor": "middle", class: "axis-num" }, String(e)));
  }
  svg.append(el("text", { x: (L + W - R) / 2, y: H - 4, "text-anchor": "middle", class: "axis-title" }, "ionisation energy / eV (model: −orbital energy)"));
  if (ref) {
    const he = x(ref.he_i_eV);
    svg.append(el("line", { x1: he, x2: he, y1: T, y2: H - B, class: "pes-hei" }));
    for (const e of [ref.first_band_vertical_eV, ref.second_band_near_eV]) {
      const xe = x(e);
      svg.append(el("line", { x1: xe, x2: xe, y1: T, y2: H - B, class: "pes-ref" }));
    }
  }
  const leg = $("pes-legend");
  if (leg) leg.textContent = ref
    ? `Black bars: model bands. Orange dashed: measured band positions (${ref.first_band_vertical_eV.toFixed(2)} eV, NIST CCCBDB; near ${ref.second_band_near_eV} eV). Grey dotted: He I limit (${ref.he_i_eV.toFixed(2)} eV).`
    : "Black bars: model bands. No measured spectrum is curated for this case yet.";
  for (const b of P.bands) {
    const lv = view.diagram.centre[b.level];
    const g = el("g", { class: "pes-band", "data-ie": b.ie_eV, "data-electrons": b.electrons });
    g.append(el("rect", { x: x(b.ie_eV) - 4, y: y(b.electrons), width: 8, height: H - B - y(b.electrons), rx: 2 }));
    g.append(richText(labelParts(lv), { x: x(b.ie_eV), y: y(b.electrons) - 5, "text-anchor": "middle", class: "pes-label" }));
    g.append(el("title", {}, `${fmt(-b.ie_eV, 2)} eV level, ${b.electrons} electrons`));
    svg.append(g);
  }
  const ul = $("pes-notes"); ul.replaceChildren();
  for (const c of P.caveats) { const li = document.createElement("li"); li.textContent = c; ul.append(li); }
}

// ---- Release 2 dials: electron stepper (D6), degenerate-set dial (D4), lone-pair dial (D5), approximation slider (D2). ----
function dials(view) {
  const E = view.electrons;
  $("e-count").textContent = String(E.count);
  $("e-charge").textContent = E.charge === 0 ? "neutral" : `charge ${E.charge > 0 ? "+" : "−"}${Math.abs(E.charge)}`;
  $("e-minus").disabled = E.count <= 0; $("e-plus").disabled = E.count >= E.max;
  $("e-reset").hidden = E.count === E.default;
  const D = view.dial, box = $("dialbox");
  box.hidden = !D;
  if (D) {
    $("dial").value = D.angle_deg; $("dial-out").textContent = `${D.angle_deg.toFixed(0)}°`;
    const fmtW = (w) => D.atoms.map((a, i) => [a, w[i]]).filter(([, x]) => x > 0.005).map(([a, x]) => `${a} ${(100 * x).toFixed(0)} %`).join(", ");
    $("dial-note").textContent = `One of ${D.components} equivalent components, rotated by ${D.angle_deg.toFixed(0)}°. ` +
      `Its energy stays ${fmt(view.ledger.total_eV, 3)} eV. This component: ${fmtW(D.shown_weights)}. Averaged over the set (unchanged by the dial): ${fmtW(D.summed_weights)}.`;
  }
  const L = view.lone_pairs, lp = $("lpbox");
  lp.hidden = !L;
  if (L) {
    $("lpdial").value = L.angle_deg; $("lp-out").textContent = `${L.angle_deg.toFixed(0)}°`;
    // Side view along the molecular plane (water lies in yz, C2 along z): x across, -z up, so the hydrogens sit
    // below and the two lone pairs show as the familiar pair of lobes above oxygen.
    L.orbitals.forEach((o, i) => miniOrbital($(`lp${i + 1}`), o.structure, [[1, 0, 0], [0, 0, -1], [0, 1, 0]]));
    const [a, b] = L.orbitals.map((o) => o.energy_eV);
    $("lp-note").textContent = `Orbital energies ${fmt(a, 2)} and ${fmt(b, 2)} eV; their sum ${fmt(L.sum_eV, 4)} eV equals the canonical pair's ${fmt(L.canonical_sum_eV, 4)} eV, ` +
      `and the total electron density is unchanged (largest change ${L.density_change_max.toExponential(0)}). ` +
      (L.angle_deg >= 45 ? "At 45° the two are equivalent: mirror images through the molecular plane, the VSEPR-style lone pairs." : "At 0° they are the canonical MOs (one a1, one b1).");
  }
  const ap = $("approx");
  if (ap) {
    const c = caseInfo(view.case);
    const show = view.approx.available;
    $("approx-box").hidden = !show;
    if (show) ap.value = Math.round(100 * (view.approx.value ?? (c && c.variant === "lecture" ? 0 : 1)));
  }
}
// Re-project a structure spec onto chosen screen axes (u right, v up, w towards the viewer), from its 3D data.
function reproject(S, u, v, w) {
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  return { ...S,
    atoms: S.atoms.map((a) => ({ ...a, x: dot(a.xyz, u), y: dot(a.xyz, v), depth: dot(a.xyz, w) })),
    lobes: S.lobes.map((l) => (l.kind === "p" ? { ...l, dir2d: [dot(l.dir3d, u), dot(l.dir3d, v)] } : l)) };
}
function miniOrbital(svg, S, axes) {
  if (axes) S = reproject(S, ...axes);
  svg.replaceChildren();
  const size = 120; svg.setAttribute("viewBox", `0 0 ${size} ${size}`);
  const xs = S.atoms.map((a) => a.x), ys = S.atoms.map((a) => a.y);
  const span = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys), 1);
  const k = (size - 44) / span, cx = (Math.max(...xs) + Math.min(...xs)) / 2, cy = (Math.max(...ys) + Math.min(...ys)) / 2;
  drawOrbital(svg, S, (x) => size / 2 + (x - cx) * k, (y) => size / 2 - (y - cy) * k, 0.45, false);
}

// ---- D10 speed test: every case twice; the second (warm) pass is reported. Run on the reference phone. ----
async function speedTest() {
  const btn = $("speed-run"), out = $("speed-out");
  btn.disabled = true; out.hidden = false;
  const saved = { case: state.case, s: state.s, level: state.level, fragment: state.fragment };
  Object.assign(state, { s: 1, level: null, fragment: null, electrons: null, approx: null, dial: 0 });
  const rows = [];
  try {
    for (const pass of [0, 1]) {
      for (const c of state.cases) {
        out.textContent = `Pass ${pass + 1} of 2: ${c.name}…`;
        state.case = c.name;
        const t0 = performance.now();
        const v = await refresh();
        const total = performance.now() - t0;
        if (pass === 1) rows.push({ name: c.name, n: v.result.levels.reduce((a, l) => a + l.multiplicity, 0), compute: v.compute_ms, total });
      }
    }
  } finally {
    Object.assign(state, saved); $("case").value = state.case; await refresh();
    btn.disabled = false;
  }
  const med = (xs) => { const a = [...xs].sort((p, q) => p - q); return a[Math.floor(a.length / 2)]; };
  const worst = rows.reduce((a, r) => (r.total > a.total ? r : a));
  const info = state.info || {};
  const text = [
    `LGO Method Explorer speed test (D10), ${new Date().toISOString().slice(0, 16).replace("T", " ")} UTC`,
    `Device: ${navigator.userAgent}`,
    `Cores: ${navigator.hardwareConcurrency || "?"}; screen ${screen.width}x${screen.height} @${devicePixelRatio}x`,
    `First load (cold): ${info.load_ms ? (info.load_ms / 1000).toFixed(1) + " s" : "?"}`,
    `Warm engine recompute: median ${med(rows.map((r) => r.compute)).toFixed(1)} ms, max ${Math.max(...rows.map((r) => r.compute)).toFixed(1)} ms`,
    `Warm full update (engine + redraw): median ${med(rows.map((r) => r.total)).toFixed(1)} ms, max ${worst.total.toFixed(1)} ms (${worst.name}, ${worst.n} functions)`,
    `Budget (spec D10): under 100 ms warm. ${worst.total < 100 ? "Met on this device." : "NOT met on this device."}`,
  ].join("\n");
  out.textContent = text; $("speed-copy").hidden = false; out.dataset.done = "1";
  return rows;
}

function bonds(view) {
  const bo = view.result.bond_orders;
  const c = bo.counting, dm = bo.density_matrix;
  const vals = Object.values(dm.per_bond);
  const dmText = vals.length ? `${Math.min(...vals).toFixed(3)}${Math.max(...vals) - Math.min(...vals) > 1e-9 ? `–${Math.max(...vals).toFixed(3)}` : ""}` : "—";
  $("bonds").dataset.counting = c.per_bond === null ? "" : c.per_bond;
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
    `engine ${view.result.version} · CONVENTIONS v${(state.info && state.info.conventions && state.info.conventions[0]) || "?"}, ${(state.info && state.info.conventions && state.info.conventions[1]) || "status unknown"}` + (view.result.warnings.length ? ` · warnings: ${view.result.warnings.join(", ")}` : "");
}

let latest = 0;
async function refresh() {
  const mine = ++latest;
  const view = await call("view", { case: state.case, s: state.s, level: state.level, fragment: state.fragment, electrons: state.electrons,
    approx: state.approx, dial: state.dial, lone_pair_angle: state.case && state.case.startsWith("H2O") ? state.lp : null });
  if (mine !== latest) return state.view;  // a newer request superseded this one
  state.view = view;
  draw(view); levelTable(view); narrate(view); ledger(view); bonds(view); pes(view); dials(view); banner(view); about(view);
  structure(view); matching(view); slotBadges(view);
  if (proposal) proposalUpdate(view);
  if ($("t-3d").checked) view3d.render(view);
  document.body.dataset.computeMs = view.compute_ms.toFixed(1);
  document.body.dataset.dials = JSON.stringify({ electrons: view.electrons.count, approx: view.approx.value, dial: view.dial ? view.dial.angle_deg : null, lp: view.lone_pairs ? view.lone_pairs.angle_deg : null });
  document.body.dataset.rendered = `${view.case}|${view.schedule.s}|${view.ledger ? (view.ledger.fragment ?? view.ledger.level) : ""}`;
  return view;
}

function select(k) { glideToken++; state.level = k; state.fragment = null; state.dial = 0; refresh(); }
function selectFragment(id) { glideToken++; state.fragment = id; state.level = null; state.dial = 0; refresh(); }

// Glide: when jumping between construction steps (step buttons, "why does it jump"), the levels slide
// to their new energies through intermediate s values instead of snapping, so the eye can follow each
// level. Every frame is a real engine evaluation of H(s); nothing is interpolated on screen. Users who
// ask for reduced motion get the jump.
let glideToken = 0;
async function glide(target, level = null, ms = 420) {
  const mine = ++glideToken;
  const from = state.s;
  const reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  if (!reduce && Math.abs(target - from) > 1e-9) {
    state.level = null; state.fragment = null;
    const t0 = performance.now();
    for (;;) {
      const u = Math.min(1, (performance.now() - t0) / ms);
      if (u >= 1) break;
      const e = u < 0.5 ? 2 * u * u : 1 - 2 * (1 - u) * (1 - u);  // ease in-out
      state.s = from + (target - from) * e;
      $("s").value = state.s;
      await refresh();
      if (mine !== glideToken) return null;
      await new Promise((r) => requestAnimationFrame(r));
    }
  }
  if (mine !== glideToken) return null;
  state.s = target; state.level = level; $("s").value = target;
  return refresh();
}

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
  view3d.hot(on ? Number(row) : null);
}

// ---- Structure view (L2): atoms, orbital lobes coloured by phase, contact lines. ----
// Orbital drawing shared by the structure view and the in-diagram sketch: lobes coloured by phase
// (red positive, blue negative), sized by |coefficient|, with p lobes along their projected direction.
// The p functions on one atom add as a vector: c_x p_x + c_y p_y + c_z p_z is a single p orbital along
// (c_x, c_y, c_z), with size |c|. Drawing them as separate crossed dumbbells would show a shape the
// orbital does not have, so the drawings merge them per atom. Display only; the ledger keeps the basis.
function displayLobes(S) {
  const out = [], byAtom = new Map();
  for (const l of S.lobes) {
    if (l.kind !== "p") { out.push(l); continue; }
    const a = byAtom.get(l.atom) || { atom: l.atom, kind: "p", v3: [0, 0, 0], v2: [0, 0] };
    l.dir3d.forEach((d, i) => { a.v3[i] += l.coefficient * d; });
    l.dir2d.forEach((d, i) => { a.v2[i] += l.coefficient * d; });
    byAtom.set(l.atom, a);
  }
  for (const a of byAtom.values()) {
    const n = Math.hypot(...a.v3);
    if (n < 1e-6) continue;
    out.push({ atom: a.atom, kind: "p", coefficient: n, dir3d: a.v3.map((v) => v / n), dir2d: a.v2.map((v) => v / n) });
  }
  return out;
}

function drawOrbital(parent, S, X, Y, scale, labels) {
  const lobes = displayLobes(S);
  const maxc = Math.max(...lobes.map((l) => Math.abs(l.coefficient)), 1e-9);
  const order = S.atoms.slice().sort((a, b) => a.depth - b.depth);
  for (const a of order) {
    for (const l of lobes.filter((l) => l.atom === a.index)) {
      const r = (8 + 30 * Math.abs(l.coefficient) / maxc) * scale;
      const x = X(a.x), y = Y(a.y);
      if (l.kind === "s") {
        parent.append(el("circle", { cx: x, cy: y, r, class: l.coefficient > 0 ? "lobe-pos" : "lobe-neg" }));
      } else {
        const [dx, dy] = l.dir2d, n = Math.hypot(dx, dy);
        if (n < 1e-6) {  // p orbital pointing at the viewer: draw it face-on
          parent.append(el("circle", { cx: x, cy: y, r: r * 0.6, class: l.coefficient > 0 ? "lobe-pos" : "lobe-neg" }));
          continue;
        }
        const ang = Math.atan2(-dy, dx) * 180 / Math.PI, len = r * Math.max(0.35, n);
        for (const sgn of [1, -1]) {
          const positive = (sgn > 0) === (l.coefficient > 0);
          parent.append(el("ellipse", { cx: x + sgn * len * Math.cos(ang * Math.PI / 180), cy: y + sgn * len * Math.sin(ang * Math.PI / 180),
                                     rx: len, ry: r * 0.42, transform: `rotate(${ang} ${x + sgn * len * Math.cos(ang * Math.PI / 180)} ${y + sgn * len * Math.sin(ang * Math.PI / 180)})`,
                                     class: positive ? "lobe-pos" : "lobe-neg" }));
        }
      }
    }
    parent.append(el("circle", { cx: X(a.x), cy: Y(a.y), r: 9 * scale, class: "atom" }));
    if (labels) parent.append(el("text", { x: X(a.x), y: Y(a.y) + 4, "text-anchor": "middle" }, a.element));
  }
}

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
  drawOrbital(svg, S, X, Y, 1, true);
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
  let gl, prog, sphere, rot = [0, 0], lastView = null, lastCase = null, buf = {}, hotRow = null;
  const VS = `#version 300 es
  in vec3 p; uniform mat4 M; uniform mat4 P; uniform vec3 c; uniform vec3 s; uniform mat3 R; out vec3 n; out vec3 col;
  // R orients the ellipsoid: its first column is the lobe's long axis, in engine coordinates.
  void main(){ vec4 w = M * vec4(R * (p * s) + c, 1.0); n = normalize(mat3(M) * (R * (p / s))); col = vec3(0); gl_Position = P * w; }`;
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
    // Open in the same orientation as the 2D structure sketch (the engine's view_rotation: rows are
    // screen x, screen y and depth towards the viewer); dragging rotates from there. Reset per case.
    if (view.case !== lastCase) { lastCase = view.case; rot = [0, 0]; }
    const [ay, ax] = rot, cy = Math.cos(ay), sy = Math.sin(ay), cx = Math.cos(ax), sx = Math.sin(ax);
    const Ry = new Float32Array([cy, 0, -sy, 0, 0, 1, 0, 0, sy, 0, cy, 0, 0, 0, 0, 1]);
    const Rx = new Float32Array([1, 0, 0, 0, 0, cx, sx, 0, 0, -sx, cx, 0, 0, 0, 0, 1]);
    const V = S.view_rotation;
    const Up = new Float32Array([V[0][0], V[1][0], V[2][0], 0, V[0][1], V[1][1], V[2][1], 0, V[0][2], V[1][2], V[2][2], 0, 0, 0, 0, 1]);
    const M = matMul(Rx, matMul(Ry, Up));
    const a = canvas.width / canvas.height, sc = 1 / ext;
    const P = new Float32Array([sc / a, 0, 0, 0, 0, sc, 0, 0, 0, 0, -sc * 0.5, 0, 0, 0, 0, 1]);
    gl.uniformMatrix4fv(gl.getUniformLocation(prog, "M"), false, M);
    gl.uniformMatrix4fv(gl.getUniformLocation(prog, "P"), false, P);
    const I3 = new Float32Array([1, 0, 0, 0, 1, 0, 0, 0, 1]);
    // Orthonormal frame whose first axis is d (column-major, as WebGL expects).
    const frame = (d) => {
      const n = Math.hypot(...d) || 1, e1 = d.map((v) => v / n);
      const a = Math.abs(e1[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
      const cr = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
      let e2 = cr(e1, a); const m = Math.hypot(...e2); e2 = e2.map((v) => v / m);
      return new Float32Array([...e1, ...e2, ...cr(e1, e2)]);
    };
    const draw = (c, s, col, R = I3) => {
      gl.uniformMatrix3fv(gl.getUniformLocation(prog, "R"), false, R);
      gl.uniform3fv(gl.getUniformLocation(prog, "c"), c); gl.uniform3fv(gl.getUniformLocation(prog, "s"), s);
      gl.uniform4fv(gl.getUniformLocation(prog, "k"), col);
      gl.drawElements(gl.TRIANGLES, sphere.idx.length, gl.UNSIGNED_SHORT, 0);
    };
    const pos = rgb(css("--pos")), neg = rgb(css("--neg")), atom = rgb(css("--muted"));
    const lobes = displayLobes(S);
    const maxc = Math.max(...lobes.map((l) => Math.abs(l.coefficient)), 1e-9);
    for (const at of S.atoms) draw(at.xyz, [0.18, 0.18, 0.18], [...atom, 1]);
    // Ledger contacts as thin rods between atoms: dark = stabilising, light = destabilising; the hovered ledger row in the warn colour.
    const maxE = Math.max(...S.contacts.map((k) => Math.abs(k.energy_eV)), 1e-9);
    const fgc = rgb(css("--fg")), lite = rgb(css("--allowed")), warn = rgb(css("--warn"));
    for (const k of S.contacts) {
      const A = S.atoms[k.atoms[0]].xyz, B = S.atoms[k.atoms[1]].xyz;
      const d = B.map((v, i) => v - A[i]), len = Math.hypot(...d);
      if (len < 1e-6) continue;
      const isHot = hotRow === k.ledger_row, w = (0.025 + 0.05 * Math.abs(k.energy_eV) / maxE) * (isHot ? 1.8 : 1);
      draw(A.map((v, i) => v + d[i] / 2), [len / 2, w, w], [...(isHot ? warn : k.style === "solid" ? fgc : lite), 1], frame(d));
    }
    for (const l of lobes) {
      const at = S.atoms[l.atom].xyz, r = 0.12 + 0.42 * Math.abs(l.coefficient) / maxc;
      if (l.kind === "s") { draw(at, [r, r, r], [...(l.coefficient > 0 ? pos : neg), 1]); continue; }
      const R = frame(l.dir3d), d = R.slice(0, 3);
      for (const sgn of [1, -1]) {
        const positive = (sgn > 0) === (l.coefficient > 0);
        const c = at.map((v, i) => v + sgn * d[i] * r);
        // Prolate ellipsoid along the lobe direction, whatever that direction is (ligand p orbitals
        // in a symmetry-adapted combination need not lie along a Cartesian axis).
        draw(c, [r, r * 0.55, r * 0.55], [...(positive ? pos : neg), 1], R);
      }
    }
  }
  // Ledger hover: highlight that contact's rod (only when the 3D view is showing).
  function hot(row) { hotRow = row; if (lastView && !$("view3d").hidden) render(lastView); }
  return { render, hot };
})();


// ---- Design proposal preview (?design=proposal). Layout only; same engine data and the same draw code. ----
const proposal = new URLSearchParams(location.search).get("design") !== "classic";
function proposalSetup(cases) {
  document.body.classList.add("proposal");
  const main = document.querySelector("main");
  // Left rail: cases grouped by family, instead of a long select.
  const rail = document.createElement("nav"); rail.className = "rail"; rail.setAttribute("aria-label", "Cases");
  const fams = [...new Set(cases.map((c) => c.family))];
  for (const f of fams) {
    const h = document.createElement("h2"); h.textContent = f; rail.append(h);
    const list = document.createElement("div"); list.className = "list";
    // One button per entry: "H2O" stands for both its model and lecture presets, chosen with the toggle.
    const seen = new Set();
    for (const c of cases.filter((c) => c.family === f)) {
      if (seen.has(c.entry)) continue;
      seen.add(c.entry);
      const b = document.createElement("button"); b.className = "case"; b.type = "button"; b.dataset.entry = c.entry;
      b.innerHTML = `${c.entry}<small>${c.group}</small>`;
      b.addEventListener("click", () => pickEntry(c.entry));
      list.append(b);
    }
    rail.append(list);
  }
  main.prepend(rail);
  // Construction timeline above the diagram: named lecture steps plus play.
  const tl = document.createElement("div"); tl.className = "timeline";
  const play = $("watch"); const steps = document.createElement("div"); steps.className = "steps"; steps.id = "steps";
  const left = document.createElement("div"); left.append(play);
  const right = document.createElement("div"); right.append(steps, document.querySelector(".slider"));
  tl.append(left, right);
  document.querySelector(".diagram").prepend(variantSetup(), tl, $("narration"), $("badges"), $("rule-detail"), document.querySelector(".toggles"));
  // Inspector as tabs, one per h2 section.
  const insp = document.querySelector(".inspector");
  const kids = [...insp.children]; const tabs = document.createElement("div"); tabs.className = "tabs"; tabs.setAttribute("role", "tablist");
  const panels = []; let cur = null;
  for (const k of kids) {
    if (k.tagName === "H2") { cur = document.createElement("div"); cur.className = "tabpanel"; cur.setAttribute("role", "tabpanel"); panels.push([k.firstChild.textContent.trim(), cur]); }
    if (cur) cur.append(k);
  }
  insp.replaceChildren(tabs, ...panels.map((p) => p[1]));
  // WAI-ARIA tabs: arrow keys, Home and End move between tabs; only the selected tab is in the Tab order.
  tabs.setAttribute("aria-label", "Inspector sections");
  const choose = (i, focus) => panels.forEach(([, p], j) => {
    const t = tabs.children[j];
    p.hidden = j !== i; t.setAttribute("aria-selected", j === i); t.tabIndex = j === i ? 0 : -1;
    if (j === i && focus) t.focus();
  });
  panels.forEach(([name, panel], i) => {
    const b = document.createElement("button"); b.type = "button"; b.setAttribute("role", "tab"); b.textContent = name;
    b.id = `tab-${i}`; panel.id = `tabpanel-${i}`;
    b.setAttribute("aria-controls", panel.id); panel.setAttribute("aria-labelledby", b.id); panel.tabIndex = 0;
    b.addEventListener("click", () => choose(i, false));
    b.addEventListener("keydown", (e) => {
      const n = panels.length, k = { ArrowRight: i + 1, ArrowLeft: i - 1 + n, Home: 0, End: n - 1 }[e.key];
      if (k === undefined) return;
      e.preventDefault(); choose(k % n, true);
    });
    tabs.append(b);
  });
  choose(0, false);
}
// Model / lecture toggle. The remembered choice carries over between entries that have both.
function caseInfo(name) { return (state.cases || []).find((c) => c.name === name); }
function pickEntry(entry, variant) {
  const opts = state.cases.filter((c) => c.entry === entry);
  const want = variant || state.variant || "model";
  const c = opts.find((o) => o.variant === want) || opts[0];
  $("case").value = c.name; $("case").dispatchEvent(new Event("change"));
}
function variantSetup() {
  const box = document.createElement("div"); box.className = "variant"; box.id = "variant"; box.hidden = true;
  const lab = document.createElement("span"); lab.id = "variant-label"; lab.textContent = "Parameters";
  const seg = document.createElement("div"); seg.className = "seg"; seg.setAttribute("role", "radiogroup"); seg.setAttribute("aria-labelledby", lab.id);
  const names = { model: "Model", lecture: "Lecture" };
  const buttons = Object.keys(names).map((v, i, all) => {
    const b = document.createElement("button"); b.type = "button"; b.setAttribute("role", "radio"); b.dataset.variant = v; b.textContent = names[v];
    b.title = v === "model" ? "Every coupling from the parameter table" : "The simplifications used on the lecture slides";
    b.addEventListener("click", () => { state.variant = v; pickEntry(caseInfo(state.case).entry, v); });
    b.addEventListener("keydown", (e) => {
      const k = { ArrowRight: 1, ArrowDown: 1, ArrowLeft: -1, ArrowUp: -1 }[e.key];
      if (!k) return;
      e.preventDefault(); const n = all[(i + k + all.length) % all.length]; buttons.find((x) => x.dataset.variant === n).click();
    });
    seg.append(b); return b;
  });
  const note = document.createElement("span"); note.className = "note"; note.id = "variant-note";
  // D2: a continuous slider between the two presets, where they differ.
  const apb = document.createElement("label"); apb.id = "approx-box"; apb.className = "approx-box"; apb.hidden = true;
  apb.innerHTML = `<span>Lecture</span><input type="range" id="approx" min="0" max="100" step="1" aria-label="Lecture to model"><span>Model</span>`;
  apb.querySelector("input").addEventListener("input", (e) => {
    const v = Number(e.target.value) / 100, c = caseInfo(state.case);
    if (v <= 0) { state.approx = null; pickEntry(c.entry, "lecture"); return; }
    if (v >= 1) { state.approx = null; pickEntry(c.entry, "model"); return; }
    state.case = `${c.entry} model`; state.approx = v; state.level = null; state.fragment = null; refresh();
  });
  box.append(lab, seg, note, apb);
  return box;
}
function variantUpdate(view) {
  const box = $("variant"); if (!box) return;
  const c = caseInfo(view.case);
  box.hidden = !(c && c.variant);
  if (!c || !c.variant) return;
  for (const b of box.querySelectorAll("button")) {
    const on = b.dataset.variant === c.variant;
    b.setAttribute("aria-checked", on); b.tabIndex = on ? 0 : -1;
  }
  $("variant-note").textContent = c.variant_note || "";
  const focusIn = box.contains(document.activeElement);
  if (focusIn) box.querySelector('[aria-checked="true"]').focus();
}

function proposalUpdate(view) {
  const cur = caseInfo(view.case);
  for (const b of document.querySelectorAll(".rail button.case")) b.setAttribute("aria-current", !!cur && b.dataset.entry === cur.entry);
  variantUpdate(view);
  const steps = $("steps"); if (!steps) return;
  const sc = view.schedule, K = sc.stages.length;
  const names = { within_set: "Fix one, flip partner", combine: "Combine equal-energy sets", between_sets: "Put the central atom back" };
  const items = [["Separate atoms", 0], ...sc.stages.map((st, k) => [names[st.kind] || st.kind, (k + 1) / K])];
  steps.replaceChildren();
  for (const [label, s] of items) {
    const b = document.createElement("button"); b.type = "button"; b.className = "step"; b.textContent = label;
    b.setAttribute("aria-current", Math.abs(view.schedule.s - s) < 1e-9);
    b.addEventListener("click", () => glide(s));
    steps.append(b);
  }
}

(async () => {
  try {
    const info = await call("info");
    state.info = info;
    $("footer").textContent = `Pyodide · Python ${info.python} · NumPy ${info.numpy} · first load ${(info.load_ms / 1000).toFixed(1)} s`;
    const cases = await call("cases");
    state.cases = cases;
    const sel = $("case");
    for (const c of cases) sel.append(new Option(`${c.name} (${c.group})`, c.name));
    if (proposal) proposalSetup(cases);
    state.case = cases.find((c) => c.name === "SF6 model").name;
    if (await fromLink()) $("s").value = state.s;
    sel.value = state.case;
    sel.addEventListener("change", () => { glideToken++; state.case = sel.value; state.level = null; state.fragment = null; state.electrons = null; state.approx = null; state.dial = 0; state.lp = 0; state.s = 1; $("s").value = 1; refresh(); });
    $("s").addEventListener("input", () => { glideToken++; state.s = Number($("s").value); state.level = null; state.fragment = null; refresh(); });
    $("whyjump").addEventListener("click", () => {
      const o = state.view.origins[state.level];
      glide(o.s, o.level_index_there).then((v) => v && status(`Jumped to where this level appeared: ${o.label}`));
    });
    $("ex-svg").addEventListener("click", exportSVG);
    $("ex-png").addEventListener("click", exportPNG);
    $("ex-json").addEventListener("click", exportJSON);
    $("ex-link").addEventListener("click", shareLink);
    $("im-json").addEventListener("change", async (e) => { const f = e.target.files[0]; if (f) await importJSON(await f.text()); e.target.value = ""; });
    $("t-3d").addEventListener("change", () => { $("view3d").hidden = !$("t-3d").checked; $("structure").toggleAttribute("hidden", $("t-3d").checked); if ($("t-3d").checked) view3d.render(state.view); });
    window.lgo = { exportJSON, importJSON, shareLink, state };
    let rz; window.addEventListener("resize", () => { clearTimeout(rz); rz = setTimeout(() => state.view && draw(state.view), 150); });
    $("watch").addEventListener("click", watch);
    $("speed-run").addEventListener("click", speedTest);
    $("speed-copy").addEventListener("click", async () => {
      try { await navigator.clipboard.writeText($("speed-out").textContent); $("speed-copy").textContent = "Copied"; }
      catch { $("speed-copy").textContent = "Select the text above to copy"; }
    });
    const stepE = (d) => { const E = state.view.electrons; state.electrons = Math.min(Math.max(E.count + d, 0), E.max); if (state.electrons === E.default) state.electrons = null; refresh(); };
    $("e-minus").addEventListener("click", () => stepE(-1));
    $("e-plus").addEventListener("click", () => stepE(+1));
    $("e-reset").addEventListener("click", () => { state.electrons = null; refresh(); });
    $("dial").addEventListener("input", (e) => { state.dial = Number(e.target.value); refresh(); });
    $("lpdial").addEventListener("input", (e) => { state.lp = Number(e.target.value); refresh(); });
    $("t-symonly").addEventListener("change", () => draw(state.view));
    for (const id of ["t-mixing", "t-allowed", "t-irreps"]) $(id).addEventListener("change", () => draw(state.view));
    for (const id of ["case", "s", "watch"]) $(id).disabled = false;
    await refresh();
    document.body.dataset.ready = "1";
  } catch (e) {
    $("banner").textContent = `Engine failed to load: ${e.message}`;
    document.body.dataset.error = e.message;
  }
})();
