// Browser acceptance tests (gate C4, Release 1 subset) against the real Pyodide 0.27.7 runtime.
// Usage: node web/tests/web_test.js [screenshot-dir]. Needs network access to cdn.jsdelivr.net.
const { chromium } = require("playwright");
const { spawn } = require("child_process");
const fs = require("fs");
const path = require("path");
const ROOT = path.resolve(__dirname, "../..");
const REF = JSON.parse(fs.readFileSync(path.join(ROOT, "evidence/reference_spectra.json")));
// Gate C4 inputs: each fixture's basis (checks/crystod_cases.json) and the on-site table, for the s = 0 endpoint;
// counting bond orders from CONVENTIONS A.5.
const CASES = Object.fromEntries(JSON.parse(fs.readFileSync(path.join(ROOT, "checks/crystod_cases.json"))).map((c) => [c.case, c]));
const ONSITE = {};
for (const line of fs.readFileSync(path.join(ROOT, "params/onsite.csv"), "utf8").trim().split(/\r?\n/).slice(1)) {
  const [el, orb, e] = line.split(","); ONSITE[`${el}|${orb.slice(-1)}`] = Number(e);
}
const COUNTING = { "allyl cation pi": 0.5, "H3+": 1 / 3, BH3: 1, CH4: 1, "H2O lecture": 1, "SF6 model": 2 / 3, "SF6 lecture": 2 / 3,
  B2: 1, C2: 2, N2: 3, O2: 2, F2: 1, CO: 3, Al2: 1, Si2: 2, P2: 3, S2: 2, Cl2: 1 };
const shots = process.argv[2];
const results = [];
// The case list is a rail of buttons; the underlying <select> stays the single source of state.
const pick = (page, n) => page.$eval("#case", (s, n) => { s.value = n; s.dispatchEvent(new Event("change")); }, n);
const check = (name, ok, detail = "") => { results.push({ name, ok, detail }); console.log(`${ok ? "PASS" : "FAIL"} ${name}${detail ? " — " + detail : ""}`); };

(async () => {
  const port = 8000 + Math.floor(Math.random() * 1000);
  const server = spawn("python3", ["-m", "http.server", String(port), "--bind", "127.0.0.1"], { cwd: ROOT, stdio: "ignore" });
  await new Promise((r) => setTimeout(r, 800));
  const proxy = process.env.HTTPS_PROXY || process.env.https_proxy;
  const browser = await chromium.launch({ proxy: proxy ? { server: proxy, bypass: "127.0.0.1,localhost" } : undefined });
  try {
    const ctx = await browser.newContext({ ignoreHTTPSErrors: true, viewport: { width: 1200, height: 900 } });
    const page = await ctx.newPage();
    const errors = [];
    page.on("pageerror", (e) => errors.push(String(e)));
    const t0 = Date.now();
    await page.goto(`http://127.0.0.1:${port}/web/`);
    await page.waitForFunction(() => document.body.dataset.ready || document.body.dataset.error, null, { timeout: 180000 });
    const err = await page.evaluate(() => document.body.dataset.error);
    check("engine loads in Pyodide", !err, err || `cold load ${((Date.now() - t0) / 1000).toFixed(1)} s`);
    await page.waitForFunction(() => document.getElementById("loading").hidden, null, { timeout: 5000 });
    const loadSteps = await page.$$eval("#loading li", (l) => l.map((x) => x.className));
    check("loading screen shows the steps and goes away once the engine is ready", loadSteps.length === 5 && loadSteps.every((c) => c === "done"), loadSteps.join(","));
    const footer = await page.textContent("#footer");
    check("runtime is Python 3.12.7 + NumPy 2.0.2", /Python 3\.12\.7 · NumPy 2\.0\.2/.test(footer), footer);

    const names = await page.$$eval("#case option", (os) => os.map((o) => o.value));
    const rail = await page.$$eval(".rail button.case", (bs) => bs.map((b) => b.dataset.entry));
    const entries = [...new Set(names.map((n) => n.replace(/ (model|lecture)$/, "")))];
    check("case rail lists every entry once (model/lecture merged), grouped by family",
          JSON.stringify(rail.slice().sort()) === JSON.stringify(entries.sort()) && rail.length === entries.length && (await page.$$eval(".rail h2", (h) => h.length)) === 4, `${rail.length} buttons for ${names.length} fixtures`);
    // Toggle: H2O switches between its presets in place; cases without variants hide it.
    await page.$eval('.rail button.case[data-entry="H2O"]', (b) => b.click());
    await page.waitForFunction(() => document.body.dataset.rendered === "H2O model|1|");
    await page.$eval('#variant [data-variant="lecture"]', (b) => b.click());
    await page.waitForFunction(() => document.body.dataset.rendered === "H2O lecture|1|");
    const vnote = await page.textContent("#variant-note");
    await page.$eval('.rail button.case[data-entry="SF6"]', (b) => b.click());
    await page.waitForFunction(() => document.body.dataset.rendered === "SF6 lecture|1|");
    await page.$eval('.rail button.case[data-entry="N2"]', (b) => b.click());
    await page.waitForFunction(() => document.body.dataset.rendered === "N2|1|");
    const vHidden = await page.$eval("#variant", (v) => v.hidden);
    check("model/lecture toggle switches presets, is remembered, and hides for single-preset cases",
          vnote.includes("O 2s–H 1s coupling is switched off") && vHidden, vnote);
    await page.$eval('.rail button.case[data-entry="SF6"]', (b) => b.click());
    await page.$eval('#variant [data-variant="model"]', (b) => b.click());
    await page.waitForFunction(() => document.body.dataset.rendered === "SF6 model|1|");
    const fonts = await page.evaluate(async () => { await document.fonts.ready; return document.fonts.check('17px "STIX Two Text"'); });
    check("diagram labels use STIX Two Text", fonts);
    // Molecule packs (Release 3) have no committed reference spectrum: the pack gate checks them (F2.15).
    const fixtures = names.filter((n) => n in REF);
    check("all reference fixtures offered", JSON.stringify([...fixtures].sort()) === JSON.stringify(Object.keys(REF).sort()), `${fixtures.length} fixtures, ${names.length - fixtures.length} packs`);
    let worst = 0, timings = [];
    for (const n of fixtures) {
      await pick(page, n);
      await page.waitForFunction((n) => document.body.dataset.rendered === `${n}|1|`, n);
      const got = await page.$$eval("[data-mo]", (hs) => hs.map((h) => Number(h.dataset.energy)));
      const want = REF[n].map((x) => x.energy_eV);
      const err = got.length === want.length ? Math.max(...got.map((g, i) => Math.abs(g - want[i]))) : Infinity;
      worst = Math.max(worst, err);
      timings.push(Number(await page.evaluate(() => document.body.dataset.computeMs)));
    }
    check("browser energies equal committed reference spectra (1e-9 eV)", worst < 1e-9, `max error ${worst.toExponential(2)}`);

    // Gate C4: slider end point s = 0 (separate atoms) equals the on-site table, every case, with multiplicities;
    // the ledger identity (contacts + on-sites = eigenvalue) at s = 1 and mid-construction; counting bond orders (A.5).
    const setS = (v) => page.$eval("#s", (el, v) => { el.value = v; el.dispatchEvent(new Event("input")); }, v);
    let s0worst = 0, s0bad = [], ledgerWorst = 0, ledgerBad = [], boBad = [];
    for (const n of fixtures) {
      await pick(page, n);
      await page.waitForFunction((n) => document.body.dataset.rendered === `${n}|1|`, n);
      const bo = await page.$eval("#bonds", (b) => b.dataset.counting);
      if (n in COUNTING && Math.abs(Number(bo) - COUNTING[n]) > 1e-12) boBad.push(`${n}: ${bo}`);
      const top = await page.$$eval("[data-mo]", (hs) => hs.length - 1);
      await page.$eval(`[data-mo="${top}"]`, (h) => h.dispatchEvent(new Event("click")));
      await page.waitForFunction(([n, k]) => document.body.dataset.rendered === `${n}|1|${k}`, [n, top]);
      for (const [label, sel] of [["s=1", "#ledger table"]]) {
        const [t, e] = await page.$eval(sel, (t) => [Number(t.dataset.total), Number(t.dataset.eigen)]);
        ledgerWorst = Math.max(ledgerWorst, Math.abs(t - e)); if (!(Math.abs(t - e) < 1e-10)) ledgerBad.push(`${n} ${label}`);
      }
      await setS(0.5);
      await page.waitForFunction((n) => document.body.dataset.rendered === `${n}|0.5|`, n);
      await page.$eval('[data-slider-level="0"]', (h) => h.dispatchEvent(new Event("click")));
      await page.waitForFunction((n) => document.body.dataset.rendered === `${n}|0.5|0`, n);
      const [t2, e2] = await page.$eval("#ledger table", (t) => [Number(t.dataset.total), Number(t.dataset.eigen)]);
      ledgerWorst = Math.max(ledgerWorst, Math.abs(t2 - e2)); if (!(Math.abs(t2 - e2) < 1e-10)) ledgerBad.push(`${n} s=0.5`);
      await setS(0);
      await page.waitForFunction((n) => document.body.dataset.rendered === `${n}|0|`, n);
      const got = (await page.$$eval("[data-slider-level]", (hs) => hs.map((h) => [Number(h.dataset.energy), Number(h.dataset.multiplicity)])))
        .flatMap(([e, m]) => Array(m).fill(e)).sort((a, b) => a - b);
      const c = CASES[n];
      const want = c.basis.map((b) => ONSITE[`${c.geometry_angstrom[b.atom_index].element}|${b.type}`]).sort((a, b) => a - b);
      const err = got.length === want.length ? Math.max(...got.map((g, i) => Math.abs(g - want[i]))) : Infinity;
      s0worst = Math.max(s0worst, err); if (!(err < 1e-12)) s0bad.push(`${n} (${got.length} vs ${want.length} levels)`);
      await setS(1);
      await page.waitForFunction((n) => document.body.dataset.rendered === `${n}|1|`, n);
    }
    check("C4 slider: s = 0 end point equals the on-site table, every case (1e-12 eV)", s0bad.length === 0, s0bad.join("; ") || `max error ${s0worst.toExponential(2)}`);
    check("C4 ledger: contacts + on-sites = eigenvalue at s = 1 and s = 0.5, every case (1e-10 eV)", ledgerBad.length === 0, ledgerBad.join("; ") || `max error ${ledgerWorst.toExponential(2)}`);
    check("C4 counting bond orders equal CONVENTIONS A.5", boBad.length === 0, boBad.join("; ") || `${Object.keys(COUNTING).length} cases`);
    check("warm recompute under 100 ms (this machine, not the reference phone)", Math.max(...timings) < 100, `max ${Math.max(...timings).toFixed(1)} ms`);

    await pick(page, "SF6 model");
    await page.waitForFunction(() => document.body.dataset.rendered === "SF6 model|1|");
    await page.click('[data-mo="0"]');
    await page.waitForFunction(() => document.body.dataset.rendered === "SF6 model|1|0");
    if (shots) await page.screenshot({ path: path.join(shots, "inspector.png"), fullPage: true });
    await page.waitForSelector("#ledger table", { state: "attached" });
    const [tot, eig] = await page.$eval("#ledger table", (t) => [Number(t.dataset.total), Number(t.dataset.eigen)]);
    check("ledger total equals eigenvalue (1e-10 eV)", Math.abs(tot - eig) < 1e-10, `${tot} vs ${eig}`);
    if (shots) await page.screenshot({ path: path.join(shots, "sf6-ledger.png"), fullPage: true });

    const rows = await page.$$eval("#ledger tr[data-row]", (t) => t.filter((r) => !r.textContent.includes("on-site")).map((r) => r.dataset.row).sort());
    const lines = await page.$$eval("#structure line[data-row]", (l) => l.map((x) => x.dataset.row).sort());
    check("L2: one structure line per ledger contact row", JSON.stringify(rows) === JSON.stringify(lines) && rows.length > 0, `${rows.length} rows, ${lines.length} lines`);
    const badges = await page.$$eval(".badge-slot .badge", (b) => b.map((x) => x.dataset.rule));
    check("L3: rule badges on ledger, matching and bond order", ["R1", "R6", "R7", "R8"].every((r) => badges.includes(r)), badges.join(" "));
    const titles = await page.$$eval("#diagram title", (t) => t.length);
    check("F5.1: provenance tooltips on diagram items", titles > 10, `${titles} tooltips`);
    const zeros = await page.$$eval("#matching td.zero", (t) => t.length);
    check("P6: matching table shows zero tallies", zeros > 0, `${zeros} zero entries`);
    const eg = await page.$$eval("[data-mo]", (hs) => hs.findIndex((h) => Math.abs(Number(h.dataset.energy) + 14.850508798504222) < 1e-6));
    await page.click(`[data-mo="${eg}"]`);
    await page.waitForFunction((eg) => document.body.dataset.rendered === `SF6 model|1|${eg}`, eg);
    await page.$eval("#whyjump", (b) => b.click());
    await page.waitForFunction(() => document.body.dataset.rendered.startsWith("SF6 model|0.75|"));
    check("P3: why-jump takes eg to checkpoint 3 of 4", true, await page.textContent("#export-status"));
    const exported = await page.evaluate(() => window.lgo.exportJSON());
    await page.evaluate(() => { window.lgo.state.case = "BH3"; });
    await page.evaluate((d) => window.lgo.importJSON(JSON.stringify(d)), exported);
    await page.waitForFunction(() => document.body.dataset.rendered.startsWith("SF6 model|0.75|"));
    check("F5.4/F5.5: exported JSON re-imports to the same view", true);
    await page.evaluate(() => window.lgo.importJSON(JSON.stringify({ format: "lgo-explorer-export", config: { case: "XeF2" } })));
    const rej = await page.textContent("#export-status");
    check("F5.5: invalid file rejected with the failing field named", rej.includes("config.case"), rej);
    const link = await page.evaluate(() => window.lgo.shareLink());
    const p2 = await ctx.newPage();
    await p2.goto(link.replace(/^https?:\/\/[^/]+/, `http://127.0.0.1:${port}`));
    await p2.waitForFunction(() => document.body.dataset.ready, null, { timeout: 180000 });
    const r2 = await p2.evaluate(() => document.body.dataset.rendered);
    check("F5.6: share link reproduces the view in a new page", r2.startsWith("SF6 model|0.75|"), r2);
    await p2.close();
    // F3.7: a ligand group orbital is selectable; ledger, structure and 3D follow; hovering a ledger row highlights its 3D rod.
    await page.$eval('[data-frag="R2"]', (h) => h.dispatchEvent(new Event("click")));
    await page.waitForFunction(() => document.body.dataset.rendered === "SF6 model|0.75|R2" || document.body.dataset.rendered.endsWith("|R2"));
    const fragLedger = await page.$eval("#ledger table", (t) => Math.abs(Number(t.dataset.total) - Number(t.dataset.eigen)));
    const fragSel = await page.$eval('[data-frag="R2"]', (h) => h.getAttribute("aria-pressed"));
    check("F3.7: a ligand group orbital is selectable, with its own ledger", fragLedger < 1e-10 && fragSel === "true", `ledger error ${fragLedger.toExponential(1)}`);
    await page.check("#t-3d");
    await page.waitForTimeout(300);
    const painted = await page.$eval("#view3d", (c) => { const g = c.getContext("webgl2"); const px = new Uint8Array(4 * 50 * 50);
      g.readPixels(c.width / 2 - 25, c.height / 2 - 25, 50, 50, g.RGBA, g.UNSIGNED_BYTE, px); return new Set(Array.from(px)).size; });
    check("F3.7: 3D view renders the selected orbital", painted > 4, `${painted} distinct byte values`);
    const snap = () => page.$eval("#view3d", (c) => { const g = c.getContext("webgl2"); const px = new Uint8Array(4 * c.width * c.height);
      g.readPixels(0, 0, c.width, c.height, g.RGBA, g.UNSIGNED_BYTE, px); let h = 0; for (let i = 0; i < px.length; i += 97) h = (h * 31 + px[i]) >>> 0; return h; });
    const before3d = await snap();
    await page.evaluate(() => { const tr = [...document.querySelectorAll("#ledger tr[data-row]")].find((t) => !t.textContent.includes("on-site")); tr.dispatchEvent(new Event("mouseenter")); });
    await page.waitForTimeout(200);
    const after3d = await snap();
    await page.evaluate(() => { const tr = [...document.querySelectorAll("#ledger tr[data-row]")].find((t) => !t.textContent.includes("on-site")); tr.dispatchEvent(new Event("mouseleave")); });
    check("F3.7: hovering a ledger contact highlights it in the 3D view", before3d !== after3d);
    const both = await page.evaluate(() => [getComputedStyle(document.querySelector("#structure")).display, getComputedStyle(document.querySelector("#view3d")).display]);
    check("3D view replaces the 2D sketch while on", both[0] === "none" && both[1] !== "none", both.join(" / "));
    if (shots) await page.screenshot({ path: path.join(shots, "inspector-3d.png"), fullPage: true });
    await page.uncheck("#t-3d");
    await page.$eval("#s", (s) => { s.value = 1; s.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => document.body.dataset.rendered === "SF6 model|1|");

    await page.$eval("#s", (s) => { s.value = 0; s.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => document.querySelector("#narration").textContent.startsWith("Separate atoms"));
    check("slider at 0 shows separate atoms", true);
    await page.$eval("#s", (s) => { s.value = 0.125; s.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => /Stage 1 of 4/.test(document.querySelector("#narration").textContent));
    check("SF6 slider runs four stages", true);
    if (shots) await page.screenshot({ path: path.join(shots, "sf6-stage1.png"), fullPage: true });

    await pick(page, "N2");
    await page.waitForFunction(() => document.body.dataset.rendered === "N2|1|");
    const verdicts = await page.$$eval("#comparisons li", (ls) => ls.map((l) => l.dataset.agrees + ":" + l.textContent.split(":")[0]));
    check("N2 panel: ordering disagrees, spin agrees", JSON.stringify(verdicts) === JSON.stringify(["0:Orbital ordering", "1:Ground-state spin"]), verdicts.join(", "));
    await pick(page, "Si2");
    await page.waitForFunction(() => document.body.dataset.rendered === "Si2|1|");
    const siTitles = await page.$$eval("#diagram .col-title", (t) => t.map((x) => x.textContent));
    const siSpin = await page.$$eval("#comparisons li", (ls) => ls.map((l) => l.dataset.agrees + ":" + l.textContent.split(":")[0]));
    check("Si2: element symbol kept in column titles; spin compared, no ordering claim", /^Si \(/.test(siTitles[0]) && JSON.stringify(siSpin) === JSON.stringify(["1:Ground-state spin"]), `${siTitles.join(" | ")}; ${siSpin}`);
    await pick(page, "Al2");
    await page.waitForFunction(() => document.body.dataset.rendered === "Al2|1|");
    const alSpin = await page.$$eval("#comparisons li", (ls) => ls.map((l) => l.dataset.agrees + ":" + l.textContent.split(":")[0]));
    const alNote = await page.textContent("#notes");
    check("Al2: spin disagreement shown, with its reading note", JSON.stringify(alSpin) === JSON.stringify(["0:Ground-state spin"]) && alNote.includes("Spin disagrees"), `${alSpin}`);
    await pick(page, "CO");
    await page.waitForFunction(() => document.body.dataset.rendered === "CO|1|");
    const coNote = await page.textContent("#notes");
    check("CO reading note explains the starred HOMO", coNote.includes("carbon lone pair"));
    if (shots) await page.locator(".diagram").screenshot({ path: path.join(shots, "co-about.png") });

    await pick(page, "CH4");
    await page.waitForFunction(() => document.body.dataset.rendered === "CH4|1|");
    const pesBands = await page.$$eval("#pes .pes-band", (gs) => gs.map((g) => Number(g.dataset.electrons)).sort((a, b) => a - b));
    const pesLegend = await page.textContent("#pes-legend");
    const pesBadge = await page.$eval("#pes", (s) => !!s.closest(".tabpanel, .inspector").querySelector(".badge.approx"));
    check("CH4 photoelectron spectrum: two bands 3:1, badged approximate, measured positions cited",
          JSON.stringify(pesBands) === "[2,6]" && pesBadge && pesLegend.includes("NIST CCCBDB"), `${pesBands}`);
    const narr = await page.textContent("#narration");
    check("CH4 shows B5 skip note", narr.includes("No ligand–ligand coupling in this model."), narr.slice(0, 120));
    // Accessibility: named level buttons, a text table of levels, keyboard tabs.
    const a11y = await page.evaluate(() => ({
      hits: document.querySelectorAll("#diagram [data-mo]").length,
      named: [...document.querySelectorAll("#diagram [data-mo]")].filter((h) => h.getAttribute("role") === "button" && h.getAttribute("aria-label")).length,
      rows: document.querySelectorAll("#levels-sr tbody tr").length,
      homo: document.querySelector("#levels-sr").textContent.includes("(HOMO)") }));
    check("a11y: every level is a named button and listed in the text table", a11y.hits > 0 && a11y.named === a11y.hits && a11y.rows === a11y.hits && a11y.homo, JSON.stringify(a11y));
    await page.focus(".tabs [role=tab][aria-selected=true]");
    await page.keyboard.press("ArrowRight");
    const tabState = await page.evaluate(() => { const t = document.activeElement; const p = document.getElementById(t.getAttribute("aria-controls"));
      return { sel: t.getAttribute("aria-selected"), idx: t.id, shown: p && !p.hidden, others: [...document.querySelectorAll(".tabs [role=tab]")].filter((x) => x.tabIndex === 0).length }; });
    check("a11y: arrow keys move between inspector tabs (roving tabindex)", tabState.sel === "true" && tabState.idx === "tab-1" && tabState.shown && tabState.others === 1, JSON.stringify(tabState));
    await page.keyboard.press("Home");
    await page.focus("#diagram [data-mo]");
    await page.keyboard.press("Enter");
    await page.waitForFunction(() => document.body.dataset.rendered === "CH4|1|0");
    const kept = await page.evaluate(() => document.activeElement.dataset.mo);
    check("a11y: Enter on a focused level selects it and focus stays on it", kept === "0", `focused data-mo=${kept}`);
    await page.$eval("#s", (s) => { s.value = 1; s.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => document.body.dataset.rendered === "CH4|1|");

    // Glide: a step jump passes through intermediate s values; with reduced motion it does not.
    const glideSeen = async () => {
      await page.evaluate(() => { window.__seen = []; window.__obs = new MutationObserver(() => window.__seen.push(document.body.dataset.rendered));
        window.__obs.observe(document.body, { attributes: true, attributeFilter: ["data-rendered"] }); });
      await page.$eval("#steps .step", (b) => b.click());
      await page.waitForFunction(() => document.body.dataset.rendered === "CH4|0|");
      const seen = await page.evaluate(() => { window.__obs.disconnect(); return window.__seen; });
      await page.$eval("#s", (s) => { s.value = 1; s.dispatchEvent(new Event("input")); });
      await page.waitForFunction(() => document.body.dataset.rendered === "CH4|1|");
      return seen.filter((r) => !/\|(0|1)\|$/.test(r)).length;
    };
    const mid = await glideSeen();
    check("step jump glides through intermediate frames", mid >= 2, `${mid} intermediate frames`);
    await page.emulateMedia({ reducedMotion: "reduce" });
    const midR = await glideSeen();
    check("reduced motion: step jump does not animate", midR === 0, `${midR} intermediate frames`);
    await page.emulateMedia({ reducedMotion: null });
    await page.$eval("[data-mo]", (h) => h.dispatchEvent(new Event("click")));
    await page.waitForFunction(() => document.body.dataset.rendered === "CH4|1|0");
    const card = await page.evaluate(() => { const r = document.querySelector("#diagram .sketch rect").getBBox();
      const hits = [...document.querySelectorAll("#diagram text, #diagram line.level")].filter((e) => { const b = e.getBBox();
        return b.x < r.x + r.width && b.x + b.width > r.x && b.y < r.y + r.height && b.y + b.height > r.y; });
      return { lobes: document.querySelectorAll("#diagram .sketch ellipse, #diagram .sketch circle, #diagram .sketch path").length, hits: hits.length }; });
    check("selected level gets a sketch card clear of bars and labels", card.lobes > 0 && card.hits === 0, JSON.stringify(card));
    await page.click("#watch");
    await page.waitForFunction(() => document.body.dataset.watchDone === "1", null, { timeout: 60000 });
    const after = await page.$eval("#s", (s) => Number(s.value));
    check("Watch mode plays to the finished diagram", after === 1);

    // ---- Release 2, Milestone A: dials keep their invariants in the browser too. ----
    const dialsNow = () => page.evaluate(() => JSON.parse(document.body.dataset.dials));
    await pick(page, "H3+");
    await page.waitForFunction(() => document.body.dataset.rendered === "H3+|1|");
    await page.click("#e-plus");
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).electrons === 3);
    const d6 = [await page.textContent("#e-charge"), await page.textContent("#banner")];
    check("D6 electron stepper: H3+ plus one electron is neutral H3 with the Jahn–Teller warning", d6[0] === "neutral" && d6[1].includes("Jahn_Teller"), d6[0]);
    await pick(page, "H2O lecture");
    await page.waitForFunction(() => document.body.dataset.rendered === "H2O lecture|1|");
    const eLec = await page.$$eval("[data-mo]", (h) => Number(h[0].dataset.energy));
    await page.$eval("#approx", (x) => { x.value = 50; x.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).approx === 0.5);
    const eMid = await page.$$eval("[data-mo]", (h) => Number(h[0].dataset.energy));
    await page.$eval("#approx", (x) => { x.value = 100; x.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => document.body.dataset.rendered === "H2O model|1|" && JSON.parse(document.body.dataset.dials).approx === null);
    const eMod = await page.$$eval("[data-mo]", (h) => Number(h[0].dataset.energy));
    check("D2 approximation slider: lecture → model passes through intermediate values and lands on the model", eMod < eMid && eMid < eLec, `${eLec.toFixed(2)} → ${eMid.toFixed(2)} → ${eMod.toFixed(2)} eV`);
    await page.$eval("#lpdial", (x) => { x.value = 45; x.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).lp === 45);
    const lpNote = await page.textContent("#lp-note");
    const lpE = [...lpNote.matchAll(/([−+]\d+\.\d+) and ([−+]\d+\.\d+) eV/g)][0];
    check("D5 lone-pair dial: at 45° two equivalent lone pairs, sum and density unchanged", lpE && lpE[1] === lpE[2] && lpNote.includes("unchanged"), lpE ? `${lpE[1]} / ${lpE[2]} eV` : lpNote.slice(0, 80));
    await pick(page, "CH4");
    await page.waitForFunction(() => document.body.dataset.rendered === "CH4|1|");
    const t2 = await page.$$eval("[data-mo]", (h) => h.findIndex((x) => x.getAttribute("aria-label").includes("3-fold")));
    await page.$eval(`[data-mo="${t2}"]`, (h) => h.dispatchEvent(new Event("click")));
    await page.waitForFunction((k) => document.body.dataset.rendered === `CH4|1|${k}`, t2);
    const tot0 = await page.$eval("#ledger table", (t) => Number(t.dataset.total));
    await page.$eval("#dial", (x) => { x.value = 37; x.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).dial === 37);
    const tot37 = await page.$eval("#ledger table", (t) => Number(t.dataset.total));
    check("D4 degenerate-set dial: rotating a t2 component leaves its energy unchanged", Math.abs(tot37 - tot0) < 1e-10, `${(tot37 - tot0).toExponential(1)} eV`);
    const symBefore = await page.$$eval("#diagram .orb-label, #diagram line.mix, #diagram line.allowed", (e) => e.map((x) => x.textContent || x.getAttribute("class")).sort().join("|"));
    await page.check("#t-symonly");
    await page.waitForTimeout(200);
    const symAfter = await page.$$eval("#diagram .orb-label, #diagram line.mix, #diagram line.allowed", (e) => e.map((x) => x.textContent || x.getAttribute("class")).sort().join("|"));
    const noTicks = await page.$$eval("#diagram .axis-num", (t) => t.length);
    check("F3.6 symmetry-only view: same labels and connections, no energy scale", symBefore === symAfter && noTicks === 0, `${noTicks} ticks`);
    await page.uncheck("#t-symonly");
    await page.click("#speed-run");
    await page.waitForFunction(() => document.querySelector("#speed-out").dataset.done === "1", null, { timeout: 120000 });
    const speed = await page.textContent("#speed-out");
    check("D10 speed test runs every case and reports the warm budget", /Warm full update/.test(speed) && /Budget/.test(speed), speed.split("\n")[5]);

    // ---- Release 2, Milestone D: morph paths (D1) and the Walsh plot (F3.8). ----
    await pick(page, "H2O model");
    await page.waitForFunction(() => document.body.dataset.rendered === "H2O model|1|");
    await page.waitForFunction(() => document.querySelectorAll("#morph-path option").length > 1);
    await page.selectOption("#morph-path", "H2O model: bend");
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).morph !== null && document.querySelectorAll("#walsh .wl").length > 0);
    const atCase = await page.$$eval("[data-mo]", (h) => h.map((x) => Number(x.dataset.energy)));
    const caseErr = Math.max(...atCase.map((g, i) => Math.abs(g - REF["H2O model"][i].energy_eV)));
    const approxHidden = await page.$eval("#approx-box", (b) => b.hidden);
    await page.$eval("#morph-u", (x) => { x.value = 1000; x.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).morph === 1);
    const linear = await page.textContent("#morph-val");
    const walshLabels = await page.$$eval("#walsh .wlabel", (t) => t.length);
    check("D1/F3.8 water bend: opens at the H2O case point (equal to the fixture), slides to 180°, Walsh plot labelled by irrep",
          caseErr < 1e-9 && linear.includes("180.0°") && walshLabels === 6 && approxHidden, `case error ${caseErr.toExponential(1)}; ${linear}`);
    await pick(page, "H3+");
    await page.waitForFunction(() => document.body.dataset.rendered === "H3+|1|" && document.querySelectorAll("#morph-path option").length > 1);
    await page.selectOption("#morph-path", "H3+: assemble");
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).morph === 1);
    const h3end = await page.$$eval("[data-mo]", (h) => h.map((x) => Number(x.dataset.energy)));
    await page.$eval("#morph-u", (x) => { x.value = 500; x.dispatchEvent(new Event("input")); });
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).morph === 0.5);
    const chain = await page.$$eval("[data-mo]", (h) => h.length);
    check("D1 H2 + H → H3+: ends exactly at the H3+ case; at the midpoint (linear H3) the e′ pair is split",
          Math.max(...h3end.map((g, i) => Math.abs(g - REF["H3+"][i].energy_eV))) < 1e-9 && chain === 3, `${chain} levels mid-path`);
    await page.selectOption("#morph-path", "");
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).morph === null);

    // ---- Release 2, Milestone E: doubt toggles (P4) and alternative schedules (P5). ----
    const dl = () => page.evaluate(() => JSON.parse(document.body.dataset.dials));
    await pick(page, "H2O lecture");
    await page.waitForFunction(() => document.body.dataset.rendered === "H2O lecture|1|");
    for (const r of ["R2", "R3"]) {
      await page.click(`#doubt-chips button[data-rule="${r}"]`);
      await page.waitForFunction((r) => (JSON.parse(document.body.dataset.dials).doubt || []).includes(r), r);
    }
    const dsteps = await page.$$eval("#steps .doubt-step", (x) => x.map((y) => y.dataset.rule));
    await page.click("#watch");
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).open === "R3", null, { timeout: 20000 });
    const r3s = Number(await page.$eval("#s", (x) => x.value));
    const r3txt = await page.textContent("#doubtbox");
    await page.click("#doubt-continue");
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).open === "R2", null, { timeout: 20000 });
    const r2txt = await page.textContent("#doubtbox");
    await page.click("#doubt-continue");
    await page.waitForFunction(() => document.body.dataset.watchDone === "1" && JSON.parse(document.body.dataset.dials).open === null, null, { timeout: 20000 });
    check("P4 doubt toggles: doubting R2 and R3 adds exactly those checkpoints; Watch stops at each with the engine's check",
          JSON.stringify(dsteps) === JSON.stringify(["R3", "R2"]) && r3s === 0 && /B5/.test(r3txt) && /fails/.test(r2txt) && /O 2pz–LGO1/.test(r2txt),
          `${dsteps} · R3 at s=${r3s}`);
    for (const r of ["R2", "R3"]) await page.click(`#doubt-chips button[data-rule="${r}"]`);
    await page.waitForFunction(() => (JSON.parse(document.body.dataset.dials).doubt || []).length === 0);
    await pick(page, "SF6 model");
    await page.waitForFunction(() => document.body.dataset.rendered === "SF6 model|1|");
    const sopts = await page.$$eval("#sched option", (o) => o.map((x) => x.value));
    await page.selectOption("#sched", "pair:1");
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).schedule === "pair:1");
    const snote = await page.textContent("#sched-note");
    check("P5 alternative schedules: SF6 offers other pairings and orders; the endpoint is unchanged",
          sopts.length >= 4 && sopts[0] === "lecture" && /same endpoint/.test(snote) && /e-1[0-9]|0\.0e\+0/.test(snote), `${sopts.length} options · ${snote}`);
    await page.selectOption("#sched", "lecture");
    await page.waitForFunction(() => JSON.parse(document.body.dataset.dials).schedule === "lecture");

    // ---- Release 2, Milestone E: compare (D9) and linked highlighting (D7). ----
    await pick(page, "allyl cation pi");
    await page.waitForFunction(() => document.body.dataset.rendered === "allyl cation pi|1|");
    await page.selectOption("#compare", JSON.stringify(["allyl cation pi", 4]));
    await page.waitForFunction(() => document.body.dataset.compare === "allyl cation pi|4" && document.querySelectorAll("#diagram2 [data-peer-mo]").length === 3);
    const axes = await page.evaluate(() => ["#diagram", "#diagram2"].map((id) => [...document.querySelectorAll(`${id} .axis-num`)].map((t) => `${t.textContent}@${Math.round(+t.getAttribute("y"))}`).join(",")));
    const cap = await page.textContent("#cap-peer");
    await page.mouse.move(0, 0);
    await page.focus('#diagram [data-key="1A2"]');
    const peerHotN = await page.$$eval("#diagram2 g.peer-hot", (g) => g.length);
    check("D9 compare: allyl cation beside the anion on one energy axis; the same level highlights in both",
          axes[0] === axes[1] && axes[0].length > 0 && /anion/.test(cap) && peerHotN === 1, `${cap} · axis ${axes[0] === axes[1] ? "shared" : "differs"}`);
    await page.selectOption("#compare", "");
    await page.waitForFunction(() => document.body.dataset.compare === "");
    // CONVENTIONS v1.5: H2S vs H2O, the BH3 -> NH3 morph, the hydrides' lecture bond order.
    await pick(page, "H2S");
    await page.waitForFunction(() => document.body.dataset.rendered === "H2S|1|");
    await page.selectOption("#compare", JSON.stringify(["H2O model", null]));
    await page.waitForFunction(() => document.body.dataset.compare === "H2O model|" && document.querySelectorAll("#diagram2 [data-peer-mo]").length === 6);
    const h2sCap = await page.textContent("#cap-peer");
    await page.selectOption("#compare", "");
    await page.waitForFunction(() => document.body.dataset.compare === "");
    await pick(page, "BH3");
    await page.waitForFunction(() => document.body.dataset.rendered === "BH3|1|" && document.querySelectorAll("#morph-path option").length > 1);
    const bh3Paths = await page.$$eval("#morph-path option", (o) => o.map((x) => x.value));
    await pick(page, "HF lecture");
    await page.waitForFunction(() => document.body.dataset.rendered === "HF lecture|1|");
    const hfBO = await page.$eval("#bonds", (b) => b.dataset.counting);
    check("v1.5: H2S beside H2O on one axis; BH3 has the pyramidalisation morph; HF lecture bond order 1",
          /H2O/.test(h2sCap) && bh3Paths.includes("BH3: pyramidalise") && Number(hfBO) === 1, `${h2sCap} · ${bh3Paths.length - 1} BH3 paths · HF ${hfBO}`);
    await pick(page, "H2O model");
    await page.waitForFunction(() => document.body.dataset.rendered === "H2O model|1|");
    await page.mouse.move(0, 0);
    await page.focus('#diagram [data-mo="1"]');
    const lk = await page.evaluate(() => ({ frags: document.querySelectorAll("#diagram g.linked [data-frag]").length,
      lines: document.querySelectorAll("#diagram line.mix.linked").length, pes: document.querySelectorAll("#pes .pes-band.linked").length,
      match: document.querySelectorAll("#matching .linked").length }));
    await page.$eval('[data-mo="1"]', (h) => h.dispatchEvent(new Event("click")));
    await page.waitForFunction(() => document.body.dataset.rendered === "H2O model|1|1");
    await page.$eval('#structure circle.atom[data-atom="1"]', (a) => a.dispatchEvent(new MouseEvent("mouseover", { bubbles: true })));
    const rowsLit = await page.$$eval("#ledger tr.linked", (r) => r.length);
    check("D7 linked highlighting: a level lights its fragment orbitals, mixing lines, spectrum band and matching cells; an atom lights its ledger rows",
          lk.frags === 2 && lk.lines === 2 && lk.pes === 1 && lk.match > 0 && rowsLit > 0, JSON.stringify({ ...lk, rowsLit }));
    await page.mouse.move(0, 0);

    // ---- Release 3, Milestone A: molecule packs and the C1 fallback (F2.14). ----
    await pick(page, "HOF");
    await page.waitForFunction(() => document.body.dataset.rendered === "HOF|1|");
    const hof = await page.evaluate(() => ({ labels: [...document.querySelectorAll("#diagram [data-mo]")].length,
      notes: document.querySelector("#notes").textContent, rail: [...document.querySelectorAll(".rail button.case")].filter((b) => /pack/.test(b.textContent)).length }));
    check("Release 3 packs: HOF (Cs, not in the library) draws with the C1 fallback; packs are marked in the rail",
          hof.labels === 6 && /C1 fallback/.test(hof.notes) && hof.rail >= 3, JSON.stringify({ ...hof, notes: hof.notes.length }));

    // ---- Release 3, Milestone B: the case matcher. ----
    const findIt = async (t) => { await page.fill("#finder-input", t); await page.press("#finder-input", "Enter");
      await page.waitForFunction((t) => document.body.dataset.match && JSON.parse(document.body.dataset.match).input === t, t); return JSON.parse(await page.evaluate(() => document.body.dataset.match)); };
    const m1 = await findIt("C3H5+"), m2 = await findIt("H2Se");
    const nearTxt = await page.textContent("#finder-card .nearest");
    await page.click("#finder-card .actions .primary");
    await page.waitForFunction(() => document.body.dataset.rendered === "H2O model|1|");
    const m3 = await findIt("allyl anion");
    await page.click("#finder-card .actions .primary");
    await page.waitForFunction(() => document.body.dataset.rendered.startsWith("allyl cation pi|") && JSON.parse(document.body.dataset.dials).electrons === 4);
    check("Release 3 matcher: C3H5+ asks for a structure; H2Se is not built (Se) and opens its nearest case; the allyl anion opens with 4 electrons",
          m1.verdict === "choose" && m2.verdict === "not_built" && m2.code === "element_parameters" && /H2O model/.test(nearTxt) && m3.verdict === "drawn",
          `${m1.verdict} · ${m2.code} · ${nearTxt}`);
    await page.click("#finder-card .close");

    // ---- Release 2, Milestone C: node markers (L5), rule check (L4), phase painting (D3). ----
    await pick(page, "allyl cation pi");
    await page.waitForFunction(() => document.body.dataset.rendered === "allyl cation pi|1|");
    const nodeLabels = await page.$$eval("#diagram .nodes-tag", (t) => t.map((x) => x.textContent));
    await page.$eval('[data-mo="1"]', (h) => h.dispatchEvent(new Event("click")));
    await page.waitForFunction(() => document.body.dataset.rendered === "allyl cation pi|1|1");
    const psi2 = await page.evaluate(() => ({
      through: document.querySelectorAll("#structure .node-through").length,
      pos: document.querySelectorAll("#structure .lobe-pos").length, neg: document.querySelectorAll("#structure .lobe-neg").length,
      rows: [...document.querySelectorAll("#ledger tr.node-through")].map((r) => r.textContent).join(" | ") }));
    check("L5 allyl ψ2: one node through C2 linked to its zero-weight ledger rows; node counts beside the levels",
          psi2.through === 1 && psi2.rows.includes("node through C2") && JSON.stringify(nodeLabels) === JSON.stringify(["2 nodes", "1 node", "0 nodes"]), `${psi2.rows.slice(0, 60)}; ${nodeLabels}`);
    check("π lobes facing the viewer keep their sign (ψ2: one positive, one negative)", psi2.pos === 1 && psi2.neg === 1, JSON.stringify(psi2));
    await page.$eval('[data-mo="2"]', (h) => h.dispatchEvent(new Event("click")));
    await page.waitForFunction(() => document.body.dataset.rendered === "allyl cation pi|1|2");
    const psi3 = await page.$$eval("#structure .node-between", (m) => m.map((x) => x.dataset.row));
    await page.$eval("#structure .node-between", (m) => m.dispatchEvent(new Event("mouseenter")));
    const linked = await page.$eval(`#ledger tr[data-row="${psi3[0]}"]`, (r) => r.classList.contains("hot") && r.textContent.includes("out-of-phase"));
    check("L5 allyl ψ3: two nodes between atoms, each linked to its out-of-phase ledger row", psi3.length === 2 && linked, `${psi3.length} markers`);
    await page.$eval('[data-mo="0"]', (h) => h.dispatchEvent(new Event("click")));
    await page.waitForFunction(() => document.body.dataset.rendered === "allyl cation pi|1|0");
    const rc = await page.$$eval("#rulecheck .rc", (r) => r.map((x) => `${x.dataset.rule}:${x.querySelector(".badge").textContent}`));
    check("L4 allyl ψ1: R6 badged approximate and R9 exact, both shown", JSON.stringify(rc) === JSON.stringify(["R6:approximate", "R9:exact"]), rc.join(", "));
    await pick(page, "H3+");
    await page.waitForFunction(() => document.body.dataset.rendered === "H3+|1|");
    await page.$eval("#paint-toggle", (b) => b.click());
    await page.waitForFunction(() => !document.querySelector("#paintbox").hidden);
    const pureAll = await page.$eval("#paintbox", (b) => b.dataset.pure);
    const clickAtom = async (i) => { const n = await page.evaluate(() => document.body.dataset.dials); await page.$eval(`#structure circle.atom[data-atom="${i}"]`, (c) => c.dispatchEvent(new MouseEvent("click", { bubbles: true }))); await page.waitForFunction((n) => document.body.dataset.dials !== n, n); };
    await clickAtom(1); await clickAtom(2); await clickAtom(2);
    const ep = await page.$eval("#paintbox", (b) => [b.dataset.pure, b.textContent]);
    check("D3 phase painting: all-red H3+ is a true MO (a1′); painting (+, −, 0) gives 100 % e′", pureAll === "true" && ep[0] === "true" && ep[1].includes("100 % 1e"), ep[1].slice(-40));
    await page.$eval("#paint-toggle", (b) => b.click());

    // ---- Release 2, Milestone B: Practice and Exam modes. ----
    await pick(page, "N2");
    await page.waitForFunction(() => document.body.dataset.rendered === "N2|1|");
    await page.click('.modeswitch [data-mode="practice"]');
    await page.waitForSelector("#quiz .qstep h3");
    const qSteps = await page.$$eval("#quiz .qstep h3", (h) => h.length);
    const diagramHidden = await page.$eval("#diagram", (d) => d.getClientRects().length === 0);
    await page.$$eval("#quiz .qcheck", (b) => b[0].click());
    await page.waitForSelector('#quiz .qres[data-step="Order the atomic orbitals"]');
    const qfb = await page.textContent('#quiz .qres[data-step="Order the atomic orbitals"]');
    check("Practice mode: five steps for N2, diagram hidden, a step checked with its reason",
          qSteps === 5 && diagramHidden && /correct|Miessler/.test(qfb), `${qSteps} steps; ${qfb.slice(0, 60)}`);
    await page.click('.modeswitch [data-mode="exam"]');
    await page.$eval("#case", (s) => { s.value = "SF6 lecture"; s.dispatchEvent(new Event("change")); });
    await page.waitForFunction(() => document.querySelector("#quiz .qtext") && document.querySelector("#quiz .qtext").textContent.includes("SF6 lecture"));
    await page.$$eval("#quiz .qhead .seg button", (b) => b[1].click());
    await page.waitForSelector('#quiz input[data-field="bonds"]');
    for (const [k, v] of [["bonding_electrons", "8"], ["antibonding_electrons", "0"], ["bonds", "6"], ["bond_order", "2/3"]]) await page.fill(`#quiz input[data-field="${k}"]`, v);
    await page.click("#quiz .qsubmit");
    await page.waitForSelector("#qresult .qscore");
    const full = await page.$eval("#qresult", (r) => r.dataset.percent);
    await page.fill('#quiz input[data-field="bond_order"]', "1");
    await page.fill('#quiz input[data-field="bonding_electrons"]', "12");
    await page.click("#quiz .qsubmit");
    await page.waitForFunction(() => document.querySelector("#qresult").dataset.percent !== "100");
    const why = await page.textContent("#quiz .qres");
    check("Exam mode, Q2(d)-style: SF6 bond order 2/3 scores 100 %; '1' is explained by the nonbonding eg pair",
          full === "100" && why.includes("not a whole number") && why.includes("1eg"), `first ${full} %`);
    await page.click('.modeswitch [data-mode="explore"]');
    await page.waitForFunction(() => getComputedStyle(document.querySelector("#diagram")).display !== "none");
    check("Explore mode restores the diagram", true);

    // Banner shows the CONVENTIONS version read from the file, not a hard-coded one.
    const convVer = fs.readFileSync(path.join(ROOT, "CONVENTIONS.md"), "utf8").match(/\*\*Version: ([0-9.]+)/)[1];
    const banner = await page.textContent("#banner");
    check("banner shows the current CONVENTIONS version", banner.includes(`CONVENTIONS v${convVer},`), banner.split("·").slice(-1)[0].trim());
    // Layout: the full diagram (not the phone geometry) from 900 px up, and nothing clipped at any width.
    const layout = [];
    for (const w of [360, 1024, 1200, 1440]) {
      await page.setViewportSize({ width: w, height: 900 });
      await page.waitForTimeout(350);
      layout.push(await page.evaluate((w) => { const s = document.querySelector("#diagram"); const vb = s.viewBox.baseVal; const bb = s.getBBox();
        return { w, compact: s.classList.contains("compact"), clipped: bb.x < vb.x - 0.5 || bb.x + bb.width > vb.x + vb.width + 0.5, scroll: document.scrollingElement.scrollWidth > w }; }, w));
    }
    check("layout: full diagram from 1024 px up, phone geometry at 360 px, nothing clipped or scrolling sideways",
          layout.every((l) => l.compact === (l.w < 900) && !l.clipped && !l.scroll), JSON.stringify(layout.filter((l) => l.compact !== (l.w < 900) || l.clipped || l.scroll)));
    await page.setViewportSize({ width: 360, height: 800 });
    await page.waitForTimeout(300);
    const sw = await page.evaluate(() => document.scrollingElement.scrollWidth);
    check("no horizontal scroll at 360 px", sw <= 360, `scrollWidth ${sw}`);
    if (shots) await page.screenshot({ path: path.join(shots, "ch4-phone.png"), fullPage: true });
    check("no page errors", errors.length === 0, errors.join("; "));
  } finally {
    await browser.close();
    server.kill();
  }
  const failed = results.filter((r) => !r.ok).length;
  console.log(`${results.length - failed}/${results.length} browser checks passed`);
  process.exit(failed ? 1 : 0);
})().catch((e) => { console.error(e); process.exit(1); });
