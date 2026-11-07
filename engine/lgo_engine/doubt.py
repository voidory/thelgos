"""P4 doubt toggles (Release 2, Milestone E).

Each rule R1-R10 the current case uses has a place on the construction slider: the moment the lecture
applies it. A rule the user marks "not convinced" becomes a discrete checkpoint there, carrying the
engine's own check of the rule (exact numbers, never a restatement) and an alternative to try.

Where each rule sits, on the selected schedule (P5):
- R1, R5, R6, R9 when the ligand set is finished (end of the last fix-and-flip / combine stage); for a
  molecule that is one set (H2, H3+, allyl pi) at s = 1, where its combinations are the MOs.
- R7 at the start of the central-atom stage; R2 at the end of the last stage that couples orbitals of
  different energy; R3 where the omitted couplings would have switched on; R4 and R8 at s = 1.
A rule is "used" only if the schedule and the case give it something to do: R5 and R7 belong to the
lecture's symmetry route (not to "every coupling at once" or "atom by atom"), R2 needs a coupling
across an energy gap, R3 needs an atom pair left uncoupled, R6 and R9 a one-orbital-per-atom set.
R10 (sandwiches) is used by no Release 1 case.
"""
import dataclasses

import numpy as np

from .hamiltonian import build_hamiltonian
from .ledger import contact_ledger
from .schedule import hamiltonian_at
from .solve import cluster_levels

ORDER = ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10"]
LIGAND_KINDS = ("within_set", "combine")


def _f(x, n=2):
    return f"{x:.{n}f}".replace("-", "−")


def _positions(system, sched, H):
    """Slider positions of the stage boundaries each rule is tied to."""
    st = sched["stages"]
    K = len(st)
    lig = [k for k, s in enumerate(st) if s["kind"] in LIGAND_KINDS]
    cen = [k for k, s in enumerate(st) if s["kind"] == "between_sets"]
    gap = [k for k, s in enumerate(st) if np.any((abs(s["matrix"]) > 1e-12) &
                                                    (abs(np.subtract.outer(np.diag(H), np.diag(H))) > 1e-9))]
    pos = {"lig_end": (max(lig) + 1) / K if lig else 0.0,
           "cen_start": cen[0] / K if cen else None,
           "gap_end": (max(gap) + 1) / K if gap else None}
    return pos


def _levels_at(H):
    e, c = np.linalg.eigh(H)
    return e, c, cluster_levels(e, H)


def uses(system, sched, H, oset, two_col, matching):
    """The rules this case uses on this schedule, each with its slider position."""
    pos = _positions(system, sched, H)
    kinds = {s["kind"] for s in sched["stages"] + sched["skipped"]}
    symmetric = bool(kinds & set(LIGAND_KINDS)) or "between_sets" in kinds
    u = {"R1": 1.0, "R4": 1.0}
    if kinds & set(LIGAND_KINDS) and symmetric:
        u["R5"] = 1.0 if two_col else pos["lig_end"]
    if oset is not None and (two_col or (kinds & set(LIGAND_KINDS))):
        s6 = 1.0 if two_col else pos["lig_end"]
        u["R6"] = s6
        u["R9"] = s6
    if pos["cen_start"] is not None and system.central_atom is not None:
        u["R7"] = pos["cen_start"]
    if pos["gap_end"] is not None and matching is not None:
        u["R2"] = pos["gap_end"]
    n = len(system.atoms)
    edges = {tuple(sorted(e)) for e in system.edges}
    skipped_lig = [s for s in sched["skipped"] if s["kind"] in LIGAND_KINDS]
    if len(edges) < n * (n - 1) // 2:
        u["R3"] = pos["lig_end"] if skipped_lig else 1.0
    if system.electrons and any(True for _ in system.edges):
        u["R8"] = 1.0
    return u


def _missing_pairs(system):
    n = len(system.atoms)
    edges = {tuple(sorted(e)) for e in system.edges}
    return [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) not in edges]


def _name(system, a):
    return f"{system.atoms[a].element}{a + 1}"


def check(rule, ctx):
    """The engine's check of one rule for one case: claim, verdict, lines of evidence, alternative."""
    system, H, sched, s = ctx["system"], ctx["H"], ctx["sched"], ctx["uses"][rule]
    Hs = hamiltonian_at(sched, s)
    out = {"rule": rule, "s": s, "claim": ctx["rules"][rule][0], "lines": [], "alternative": None}
    if rule == "R4":
        n = len(system.basis)
        m = sum(lv["multiplicity"] for lv in ctx["result"]["levels"])
        out["verdict"] = "exact"
        out["lines"] = [f"{n} valence basis functions in, {m} molecular orbitals out "
                        f"({len(ctx['result']['levels'])} levels, counting each degenerate component)."]
    elif rule == "R1":
        e, c, lv = _levels_at(Hs)
        rows = []
        for ix in lv[:14]:
            L = contact_ledger(Hs, c[:, ix[0]])
            stab = sum(r["energy_eV"] for r in L["rows"] if r["kind"] == "contact" and r["energy_eV"] < 0)
            dest = sum(r["energy_eV"] for r in L["rows"] if r["kind"] == "contact" and r["energy_eV"] > 0)
            ons = sum(r["energy_eV"] for r in L["rows"] if r["kind"] == "onsite")
            rows.append([_f(float(e[ix[0]])), _f(ons), _f(stab), _f(dest),
                         "lowered" if stab + dest < -1e-6 else ("raised" if stab + dest > 1e-6 else "unchanged")])
        out["verdict"] = "exact"
        out["lines"] = ["Every same-colour contact (facing lobes of one sign) adds a negative energy and every opposite-colour "
                        "contact a positive one; a level sits below its starting energy exactly when the same-colour total wins. "
                        "At this point of the construction:"]
        out["table"] = {"head": ["Level (eV)", "Start (on-site)", "Same colour", "Opposite colour", "Net"], "rows": rows}
        if ctx["paint"]:
            out["alternative"] = {"kind": "paint", "label": "Paint your own colours and watch the energy"}
    elif rule == "R5":
        lgos = ctx["diagram"].get("right") or []
        if len(system.atoms) == 2 and system.central_atom is None:
            out["lines"] = ["Fix one atom and flip its partner: each orbital pairs with its copy on the other atom, in phase (g) "
                            "and out of phase (u). The pairs are 2s with 2s, 2pz with 2pz (σ), and 2px, 2py with their partners (π); "
                            "s and pz share σ symmetry, so those two pairs then mix."]
        elif ctx["two_col"] or not lgos:
            e, c, lv = _levels_at(Hs)
            names = [f"{_f(float(e[ix[0]]))} eV" + (f" (×{len(ix)})" if len(ix) > 1 else "") for ix in lv]
            out["lines"] = ["The symmetry combinations of the whole set are its molecular orbitals: " + ", ".join(names) + "."]
        else:
            out["lines"] = ["Fixing one atom and flipping its partner, then combining equal-energy sets, gives the same "
                            "orbitals as projecting onto each irreducible representation:",
                            ", ".join(f"{x['label']} ({x['irrep']}, {_f(x['energy_eV'])} eV" + (f", ×{x['multiplicity']})" if x['multiplicity'] > 1 else ")")
                                      for x in lgos) + "."]
        out["lines"].append("Exact. A degenerate set may come out rotated relative to the textbook's choice; any rotation is equally valid.")
        out["verdict"] = "exact"
        if len(ctx["options"]) > 1:
            out["alternative"] = {"kind": "schedule", "label": "Try another pairing or order (P5): the endpoint does not change"}
        elif ctx["paint"]:
            out["alternative"] = {"kind": "paint", "label": "Paint a combination that is not an LGO and see it decompose"}
    elif rule in ("R6", "R9"):
        rows, ok = [], True
        for label, (E, v), Vs in zip(ctx["set_names"], ctx["set_levels"], ctx["set_vecs"]):
            rc = ctx["rule_checks"](v)
            r = rc[rule]
            if r is None:
                continue
            if rule == "R6":
                # A degenerate level: the counting rule needs a ±1 pattern, which may be another member of the set.
                tries = [ctx["rule_checks"](Vs[:, k])["R6"] for k in range(Vs.shape[1])]
                good = next((t for t in tries if t["agrees"]), None)
                if good is not None and not r["agrees"]:
                    note = (f"agrees for the ±1 member of the degenerate set; the member shown elsewhere gives "
                            f"{_f(r['prediction_eV'])} eV: {r['reason']}")
                    r = good
                else:
                    note = "agrees" if r["agrees"] else (r["reason"] or "differs")
                rows.append([label, _f(r["prediction_eV"]), _f(r["exact_eV"]), note])
                ok &= r["agrees"]
            else:
                rows.append([label, "—" if r["nodes"] is None else str(r["nodes"]), _f(E), r["badge"]])
                ok &= r["agrees"] is not False
        if rule == "R6":
            out["lines"] = ["Counting rule E = ε + (2/N) Σ ±t against the exact energy of each combination:"]
            out["table"] = {"head": ["Combination", "R6 predicts (eV)", "Exact (eV)", ""], "rows": rows}
            out["verdict"] = "exact here" if ok else "approximate here"
        else:
            topo = ctx["oset"]["topology"]
            out["lines"] = [{"chain": "Open chain: more nodes means higher energy, exactly (oscillation theorem). It gives the order, not the gaps.",
                             "ring": "Ring: the node order is exact up to degenerate pairs. It gives the order, not the gaps.",
                             "none": "No contacts within the set, so nodes cost nothing and every combination has the same energy.",
                             }.get(topo, "A 3D set: node counting is only a heuristic for the order.")]
            out["table"] = {"head": ["Combination", "Nodes", "Exact energy (eV)", ""], "rows": rows}
            out["verdict"] = "exact here" if (ok and topo in ("chain", "ring")) else ("not applied" if topo == "none" else "heuristic")
        if ctx["paint"]:
            out["alternative"] = {"kind": "paint", "label": "Paint a pattern and compare its count with the exact energy"}
    elif rule == "R7":
        m = ctx["matching"]
        pairs, zero_same, cross = [], [], 0
        for i, r in enumerate(m["rows"]):
            for j, col in enumerate(m["cols"]):
                t = m["tally_eV"][i][j]
                if t != 0:
                    pairs.append(f"{r['label']}–{col['label']} ({_f(t)} eV)")
                    cross += r["irrep"] != col["irrep"]
        out["lines"] = ["Back at its own position, the central atom's orbitals meet each LGO with the signed sum of their lobe "
                        "contacts. Non-zero tallies: " + (", ".join(pairs) or "none") + ".",
                        "Every other pair sums to exactly zero: the contacts cancel by symmetry." +
                        (" All non-zero tallies join orbitals of the same irrep." if not cross else
                         f" {cross} non-zero tally joins different irreps (a symmetry error).")]
        out["verdict"] = "exact" if not cross else "violated"
        out["alternative"] = {"kind": "toggle", "target": "t-allowed", "label": "Show every symmetry-allowed line"}
    elif rule == "R2":
        m = ctx["matching"]
        left = ctx["diagram"].get("left") or []
        right = ctx["diagram"].get("right") or []
        el = {r["label"]: None for r in m["rows"]}
        def energy(lbl, side):
            for x in side:
                if lbl == x["label"] or lbl.startswith(x["label"]) or ("components" in x and lbl in x["components"]) or lbl in x["label"].split("/"):
                    return x["energy_eV"]
            return None
        rows, seen = [], set()
        for i, r in enumerate(m["rows"]):
            for j, col in enumerate(m["cols"]):
                t = m["tally_eV"][i][j]
                if t == 0:
                    continue
                ea, eb = energy(r["label"], left), energy(col["label"], right)
                if ea is None or eb is None:
                    continue
                d = eb - ea
                key = (*sorted((round(ea, 6), round(eb, 6))), round(abs(t), 6))
                if key in seen:
                    continue
                seen.add(key)
                if abs(d) < 1e-9:
                    rows.append([f"{r['label']}–{col['label']}", _f(abs(d)), _f(abs(t)), "equal energies", _f(abs(t)), "0.50"])
                    continue
                ex = np.sqrt(d * d / 4 + t * t) - abs(d) / 2
                w_ex = 0.5 * (1 - abs(d) / np.sqrt(d * d + 4 * t * t))
                if abs(t / d) > 0.5:
                    est_s = est_w = "fails"
                else:
                    est_s, est_w = _f(t * t / abs(d)), f"{(t / d) ** 2:.2f}"
                rows.append([f"{r['label']}–{col['label']}", _f(abs(d)), _f(abs(t)), _f(abs(t / d)), f"{est_s} / {_f(ex)}",
                             f"{est_w} / {w_ex:.2f}"])
        out["lines"] = ["The gap rule is a yes/no version of mixing that never switches off. For each coupled pair, taken two at a "
                        "time: gap, coupling, their ratio, the shift (perturbative t²/Δ against the two-level exact value) and "
                        "how much of the partner mixes in (t²/Δ² against exact). Where |t| is more than half the gap the "
                        "perturbative estimate fails outright: the pair is strongly mixed, however large the gap looks."]
        out["table"] = {"head": ["Pair", "Gap (eV)", "|t| (eV)", "|t|/gap", "Shift: estimate / exact (eV)", "Admixture: estimate / exact"],
                        "rows": rows}
        res = ctx["result"]["levels"]
        mix = []
        names = ctx["mo_names"]
        for i, f in enumerate(system.basis):
            ws = [(k, lv["composition_ao"][i] * lv["multiplicity"]) for k, lv in enumerate(res)]
            big = [(k, w) for k, w in ws if w > 0.05]
            if len(big) > 1 and (system.central_atom is None and f.atom == 0 or f.atom == system.central_atom):
                mix.append(f"{ctx['ao_name'](f)} spreads over " + ", ".join(f"{names[k]} ({w:.0%})" for k, w in big))
        if mix:
            out["lines"].append("In the full calculation: " + "; ".join(mix[:6]) + ".")
        out["verdict"] = "approximation"
        if ctx["approx"]:
            out["alternative"] = {"kind": "approx", "label": "Slide from the lecture's on/off choice to actual mixing (D2)"}
        elif ctx["sweep"]:
            out["alternative"] = {"kind": "morph", "target": ctx["sweep"], "label": "Sweep the central energy through resonance"}
    elif rule == "R3":
        miss = _missing_pairs(system)
        pos = system.positions()
        try:
            full = dataclasses.replace(system, edges=sorted({tuple(sorted(e)) for e in system.edges} | set(miss)))
            Hf = build_hamiltonian(full, ctx["params"])
            Hb = build_hamiltonian(system, ctx["params"])
            Hw = H + (Hf - Hb)
        except Exception:   # no tabulated coupling for some pair
            Hw = None
        groups_ = {}
        for a, b in miss:
            d = float(np.linalg.norm(pos[a] - pos[b]))
            big = 0.0
            if Hw is not None:
                ia = [i for i, f in enumerate(system.basis) if f.atom == a]
                jb = [j for j, f in enumerate(system.basis) if f.atom == b]
                big = max((abs(Hw[i, j]) for i in ia for j in jb), default=0.0)
            key = ("–".join(sorted((system.atoms[a].element, system.atoms[b].element))), round(d, 2), round(big, 2))
            groups_.setdefault(key, []).append(f"{_name(system, a)}–{_name(system, b)}")
        bits = []
        for (pair, d, big), names in groups_.items():
            who = names[0] if len(names) == 1 else f"{len(names)} {pair} pairs ({', '.join(names[:4])}{' …' if len(names) > 4 else ''})"
            bits.append(f"{who} at {d:.2f} Å" + (f", |coupling| up to {_f(big)} eV each if switched on" if Hw is not None else ""))
        why = " The lecture leaves ligand–ligand couplings out (decision B5)." if any(
            st["kind"] in LIGAND_KINDS for st in sched["skipped"]) else ""
        out["lines"] = ["Left out as too far apart: " + "; ".join(bits) + "." + why]
        if Hw is not None:
            e0 = np.linalg.eigvalsh(H)
            e1 = np.linalg.eigvalsh(Hw)
            from .result import run
            r1 = run(system, ctx["params"], ctx["groups"], H=Hw)
            lab0 = [lv["irrep"] for lv in ctx["result"]["levels"]]
            lab1 = [lv["irrep"] for lv in r1["levels"]]
            same = lab0 == lab1
            out["lines"].append(f"Switching them on at their distance-law values moves the levels by up to {_f(float(np.max(abs(e1 - e0))))} eV; "
                                + ("the order of the labelled levels does not change." if same else
                                   "the order of the labelled levels changes: " + ", ".join(lab1) + " instead of " + ", ".join(lab0) + "."))
            out["verdict"] = "qualitative: unchanged order" if same else "qualitative: the order changes"
        else:
            out["verdict"] = "qualitative"
        if ctx["ring"]:
            out["alternative"] = {"kind": "morph", "target": ctx["ring"], "label": "Switch the missing contact on gradually"}
    elif rule == "R8":
        bo = ctx["result"]["bond_orders"]
        c = bo["counting"]
        dm = bo["density_matrix"]
        vals = list(dm["per_bond"].values())
        out["lines"] = [f"Counting: ({c['bonding_electrons']} bonding − {c['antibonding_electrons']} antibonding) / 2 over {c['bonds']} "
                        f"equivalent bond{'s' if c['bonds'] != 1 else ''} = {_f(c['per_bond'], 3) if c['per_bond'] is not None else '—'} per bond.",
                        f"From the density matrix ({dm['definition']}): " + (", ".join(_f(v, 3) for v in vals[:6]) if vals else "—") + " per bond.",
                        "Both are conventions, not observables; they need not agree, and the counting one ignores how bonding each orbital is."]
        out["verdict"] = "convention"
        out["alternative"] = {"kind": "electrons", "label": "Add or remove electrons and watch the bond order"}
    else:
        out["verdict"] = "not used"
    return out


def doubt_payload(rules_table, ctx, doubted):
    used = ctx["uses"]
    cps = [check(r, {**ctx, "rules": rules_table}) for r in ORDER if r in used and r in (doubted or [])]
    cps.sort(key=lambda x: (x["s"], ORDER.index(x["rule"])))
    return {"used": [r for r in ORDER if r in used], "positions": {r: used[r] for r in ORDER if r in used}, "checkpoints": cps}
