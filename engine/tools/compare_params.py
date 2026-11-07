"""Compare every case's qualitative picture under two parameter sets.

Usage: python engine/tools/compare_params.py OLD_ROOT [NEW_ROOT]

Each root must contain params/ (groups/ and evidence/ are read from this repository). To compare
against an earlier commit: `git archive <commit> params | tar -x -C /some/dir`, then pass /some/dir.

Reported as a change ("flip"):
- the energy order of the finished levels (by irrep);
- electrons per level, HOMO/LUMO, and open-shell warnings;
- bonding/antibonding character (the star);
- comparisons with experiment and the counting bond order;
- construction stages run or skipped;
- the energy order at every construction checkpoint;
- for each MO, which fragment orbital dominates it, with both atoms of a diatomic summed (so a
  50/50 left/right tie is not reported).

Energies themselves are not reported; they always move when parameters change.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from load import load_params, load_groups, load_extras, ROOT  # noqa: E402
from lgo_engine import app_api  # noqa: E402


def summary(root):
    P, G, X = load_params(root), load_groups(), load_extras()
    out = {}
    for c in app_api.case_list(P):
        v = app_api.view(P, G, c["name"], extras=X)
        r, d = v["result"], v["diagram"]
        lv = r["levels"]
        lab = lambda x: x.get("linear_label") or x["irrep"]  # noqa: E731
        steps = []
        for cp in v["schedule"]["checkpoints"]:
            w = app_api.view(P, G, c["name"], s=cp["s"])
            sl = w["diagram"]["slider_levels"] if cp["s"] < 1 else w["diagram"]["centre"]
            steps.append([x["multiplicity"] for x in sorted(sl, key=lambda x: x["energy_eV"])])
        frag = {x["id"]: (x.get("label") or x.get("irrep") or x["id"]) for x in d["left"] + d["right"]}
        agg = {}
        for ln in d["lines"]:
            if ln["layer"] != "mixing":
                continue
            m = ln["to"] if str(ln["to"]).startswith("M") else ln["from"]
            f = ln["from"] if m == ln["to"] else ln["to"]
            agg.setdefault(m, {}).setdefault(frag.get(f, f), 0.0)
            agg[m][frag.get(f, f)] += ln["weight"]
        out[c["name"]] = {
            "level order": [lab(x) for x in lv],
            "electrons": [x["electrons"] for x in lv],
            "HOMO": lab(lv[r["frontier"]["HOMO"]]) if r["frontier"]["HOMO"] is not None else None,
            "LUMO": lab(lv[r["frontier"]["LUMO"]]) if r["frontier"]["LUMO"] is not None else None,
            "warnings": sorted(r["warnings"]),
            "character": [x.get("character") for x in d["centre"]],
            "comparisons": [(x["category"], x["agrees"]) for x in v["comparisons"]],
            "counting bond order": r["bond_orders"]["counting"]["per_bond"],
            "stages": ([s["kind"] for s in v["schedule"]["stages"]], [s["kind"] for s in v["schedule"]["skipped"]]),
            "checkpoint order (multiplicities, lowest first)": steps,
            "dominant fragment": {m: max(w.items(), key=lambda t: t[1]) for m, w in sorted(agg.items())},
        }
    return out


def main(old, new):
    A, B = summary(Path(old)), summary(Path(new))
    flips = 0
    for name in B:
        a, b = A.get(name), B[name]
        if a is None:
            print(f"{name}: only in the new set")
            continue
        changed = []
        for key in b:
            if key == "dominant fragment":
                diff = {m: (a[key].get(m), b[key][m]) for m in b[key]
                        if a[key].get(m, (None,))[0] != b[key][m][0]}
                if diff:
                    changed.append((key, "; ".join(f"{m}: {x[0] if x else None} {x[1] if x else 0:.0%} -> {y[0]} {y[1]:.0%}"
                                                   for m, (x, y) in diff.items())))
            elif a[key] != b[key]:
                changed.append((key, f"{a[key]} -> {b[key]}"))
        print(f"{name}: " + ("no qualitative change" if not changed else ""))
        for key, text in changed:
            flips += 1
            print(f"    {key}: {text}")
    print(f"\n{flips} change(s) across {len(B)} cases")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else ROOT)
