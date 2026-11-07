"""Entry points the browser worker calls. Pure functions returning JSON-ready dicts."""
import numpy as np

from . import cases
from .hamiltonian import build_hamiltonian
from .result import run
from .schedule import build_schedule, hamiltonian_at, checkpoints, lambdas
from .symmetry import Group
from .diagram import diagram_spec
from .ledger import contact_ledger
from .solve import cluster_levels
from . import fidelity
from . import views
from .symmetry import projectors


def _family(s):
    if len(s.atoms) == 2 and s.central_atom is None:
        return "Diatomics"
    if s.central_atom is not None:
        return "Central atom with ligands"
    if s.name.startswith(("six H", "four H")):
        return "Ligand sets on their own"
    return "Rings and chains"


VARIANTS = ("model", "lecture")


def _entry(name):
    """'H2O lecture' -> ('H2O', 'lecture'); cases without parameter variants -> (name, None)."""
    base, _, last = name.rpartition(" ")
    return (base, last) if last in VARIANTS and base else (name, None)


def _variant_note(params, preset, other):
    """What the lecture preset changes relative to the model one, read from the presets themselves."""
    a, b = params.presets[preset], params.presets[other]
    keys = ("coupling_overrides", "edge_policy", "geometry_overrides", "hydrogen_scale")
    if all(a.get(k) == b.get(k) for k in keys):
        return "In this release the two presets have identical parameters, so the diagrams are the same."
    bits = []
    for o in (a if a["provenance_tag"] == "lecture" else b).get("coupling_overrides", []):
        na = params.orbital_name.get((o["element_a"], o["orbital_a"]), o["orbital_a"])
        nb = params.orbital_name.get((o["element_b"], o["orbital_b"]), o["orbital_b"])
        what = "switched off" if o["value_eV"] == 0 else f"set to {o['value_eV']:g} eV"
        bits.append(f"the {o['element_a']} {na}–{o['element_b']} {nb} coupling is {what}")
    for k in ("edge_policy", "geometry_overrides", "hydrogen_scale"):
        if a.get(k) != b.get(k):
            bits.append(f"different {k.replace('_', ' ')}")
    return "In the lecture version " + "; ".join(bits) + "."


def case_list(params):
    out = []
    systems = cases.release1(params)
    presets = {s.name: s.preset for s in systems}
    for s in systems:
        entry, variant = _entry(s.name)
        row = {"name": s.name, "preset": s.preset, "provenance_tag": params.presets[s.preset]["provenance_tag"],
               "group": s.group, "electrons": s.electrons, "family": _family(s), "entry": entry, "variant": variant}
        if variant:
            other = f"{entry} {'lecture' if variant == 'model' else 'model'}"
            if other in presets:
                row["variant_note"] = _variant_note(params, s.preset, presets[other])
        out.append(row)
    return out


def _system(params, name):
    for s in cases.release1(params):
        if s.name == name:
            return s
    raise ValueError(f"unknown case {name!r}")


def view(params, groups, name, s=1.0, level=None, conventions_sha256=None, extras=None, fragment=None):
    """Everything one screen needs: result, diagram, schedule, slider state and ledger."""
    system = _system(params, name)
    result = run(system, params, groups, conventions_sha256)
    H = np.array(result["hamiltonian_eV"])
    sched = build_schedule(system, Group(groups[system.group]) if system.group != "C1" else None, H)
    Hs = hamiltonian_at(sched, s)
    diagram = diagram_spec(system, H, result, groups, H_now=Hs, params=params)
    frag_vectors = diagram.pop("_fragment_vectors", {})
    K = len(sched["stages"])
    lam = lambdas(s, K)
    active = next((k for k, x in enumerate(lam) if 0 < x < 1), None)
    ledger = structure = None
    if level is not None:
        # At s < 1 the ledger and structure describe level `level` of H(s) (P3: "that moment's ledger").
        Huse = Hs if s < 1 else H
        e_s, c_s = np.linalg.eigh(Huse)
        lv_s = cluster_levels(e_s, Huse)
        ix = lv_s[min(level, len(lv_s) - 1)]
        c = c_s[:, ix[0]]
        k = int(np.flatnonzero(abs(c) > 1e-8)[0]) if np.any(abs(c) > 1e-8) else 0
        c = c if c[k] >= 0 else -c
        ledger = contact_ledger(Huse, c)
        ledger["level"] = level
        ledger["eigenvalue_eV"] = float(e_s[ix[0]])
        ledger["multiplicity"] = len(ix)
        ledger["basis_labels"] = [f.label for f in system.basis]
        structure = views.structure_spec(system, Huse, c, ledger)
    elif fragment is not None:
        # A selected fragment orbital (an LGO, or a symmetry-adapted central-atom/atom orbital): its ledger and
        # structure, from the full Hamiltonian. Only one component of a degenerate set is shown.
        if fragment not in frag_vectors:
            raise ValueError(f"unknown fragment orbital {fragment!r}")
        c = np.array(frag_vectors[fragment])[:, 0]
        c = c / np.linalg.norm(c)
        k = int(np.flatnonzero(abs(c) > 1e-8)[0])
        c = c if c[k] >= 0 else -c
        info = next(x for side in ("left", "right") for x in diagram.get(side, []) if x["id"] == fragment)
        ledger = contact_ledger(H, c)
        ledger.update({"level": None, "fragment": fragment, "fragment_label": info["label"], "eigenvalue_eV": info["energy_eV"],
                       "multiplicity": info["multiplicity"], "basis_labels": [f.label for f in system.basis]})
        structure = views.structure_spec(system, H, c, ledger)
    origin = views.origins(sched, result["levels"])
    for st in sched["stages"]:
        st["rules"] = views.STAGE_RULES.get(st["kind"], [])
    fr = result["frontier"]
    notes = []
    if extras:
        n = extras["notes"]
        notes = list(n.get("all", []))
        if any(a.element == "H" for a in system.atoms):
            notes += n.get("has_hydrogen", [])
        notes += n.get("cases", {}).get(name, [])
    comparisons = fidelity.compare(system, result, extras["curated"]) if extras else []
    return {
        "case": name,
        "notes": notes,
        "comparisons": comparisons,
        "preset": {"id": system.preset, "provenance_tag": params.presets[system.preset]["provenance_tag"]},
        "result": {k: result[k] for k in ("levels", "frontier", "bond_orders", "warnings", "checks", "version", "conventions_sha256")},
        "diagram": diagram,
        "schedule": {"stages": [{"kind": st["kind"], "edges": st["edges"], "narration": st["narration"], "rules": st.get("rules", [])} for st in sched["stages"]],
                     "skipped": [{"kind": st["kind"], "edges": st["edges"], "narration": st["narration"], "note": st["note"]} for st in sched["skipped"]],
                     "checkpoints": checkpoints(sched), "s": s, "lambdas": lam, "active_stage": active},
        "ledger": ledger,
        "structure": structure,
        "origins": origin,
        "rules": views.rules_payload(),
        "explanation": views.explanation(system, result, diagram.get("matching"), (fr["HOMO"], fr["LUMO"])),
        "provenance": views.provenance(system, params),
        "pes": views.pes_spec(system, result, extras["curated"] if extras else None),
    }


def export(params, groups, name, s=1.0, level=None, conventions_sha256=None):
    """F5.4: config, full Result, software version and conventions hash, all as JSON."""
    from . import __version__
    system = _system(params, name)
    return {"format": "lgo-explorer-export", "format_version": 1,
            "config": {"case": name, "preset": system.preset, "s": s, "level": level},
            "software": {"engine": __version__}, "conventions_sha256": conventions_sha256,
            "result": run(system, params, groups, conventions_sha256)}


def validate_import(params, data, conventions_sha256=None):
    """F5.5: reject malformed files with the failing field named; report version or hash drift."""
    if not isinstance(data, dict) or data.get("format") != "lgo-explorer-export":
        raise ValueError("format: not an LGO Explorer export")
    cfg = data.get("config")
    if not isinstance(cfg, dict):
        raise ValueError("config: missing")
    names = {c["name"] for c in case_list(params)}
    if cfg.get("case") not in names:
        raise ValueError(f"config.case: unknown case {cfg.get('case')!r}")
    s = cfg.get("s", 1.0)
    if not isinstance(s, (int, float)) or not 0 <= s <= 1:
        raise ValueError("config.s: must be a number between 0 and 1")
    lvl = cfg.get("level")
    if lvl is not None and (not isinstance(lvl, int) or lvl < 0):
        raise ValueError("config.level: must be a non-negative integer or null")
    notes = []
    if conventions_sha256 and data.get("conventions_sha256") != conventions_sha256:
        notes.append("This file was made under a different CONVENTIONS version; the diagram is recomputed with the current one.")
    return {"config": {"case": cfg["case"], "s": float(s), "level": lvl}, "notes": notes}
