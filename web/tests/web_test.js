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
    const footer = await page.textContent("#footer");
    check("runtime is Python 3.12.7 + NumPy 2.0.2", /Python 3\.12\.7 · NumPy 2\.0\.2/.test(footer), footer);

    const names = await page.$$eval("#case option", (os) => os.map((o) => o.value));
    const rail = await page.$$eval(".rail button.case", (bs) => bs.map((b) => b.dataset.entry));
    const entries = [...new Set(names.map((n) => n.replace(/ (model|lecture)$/, "")))];
    check("case rail lists every entry once (model/lecture merged), grouped by family",
          JSON.stringify(rail.slice().sort()) === JSON.stringify(entries.sort()) && rail.length === names.length - 2 && (await page.$$eval(".rail h2", (h) => h.length)) === 4, `${rail.length} buttons for ${names.length} fixtures`);
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
    check("all reference fixtures offered", JSON.stringify([...names].sort()) === JSON.stringify(Object.keys(REF).sort()), `${names.length} cases`);
    let worst = 0, timings = [];
    for (const n of names) {
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
    for (const n of names) {
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
    await page.check("#t-3d");
    await page.waitForTimeout(300);
    const painted = await page.$eval("#view3d", (c) => { const g = c.getContext("webgl2"); const px = new Uint8Array(4 * 50 * 50);
      g.readPixels(c.width / 2 - 25, c.height / 2 - 25, 50, 50, g.RGBA, g.UNSIGNED_BYTE, px); return new Set(Array.from(px)).size; });
    check("F3.7: 3D view renders the selected orbital", painted > 4, `${painted} distinct byte values`);
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
