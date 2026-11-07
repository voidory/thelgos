"""Entry points the browser worker calls. Pure functions returning JSON-ready dicts."""
import dataclasses

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


def _atom_weights(system, V):
    """Per-atom weight of one orbital (V a vector) or averaged over a degenerate set (V a matrix)."""
    V = np.atleast_2d(np.asarray(V).T).T
    w = np.zeros(len(system.atoms))
    for i, f in enumerate(system.basis):
        w[f.atom] += float(np.sum(V[i, :] ** 2)) / V.shape[1]
    return [float(x) for x in w]


def _dialled(system, V, angle_deg):
    """D4: one member of a degenerate set, rotated by angle_deg between its first two components."""
    t = np.deg2rad(angle_deg or 0.0)
    c = V[:, 0] if V.shape[1] < 2 else np.cos(t) * V[:, 0] + np.sin(t) * V[:, 1]
    c = c / np.linalg.norm(c)
    k = int(np.flatnonzero(abs(c) > 1e-8)[0]) if np.any(abs(c) > 1e-8) else 0
    c = c if c[k] >= 0 else -c
    info = None
    if V.shape[1] > 1:
        info = {"angle_deg": float(angle_deg or 0.0), "components": int(V.shape[1]),
                "shown_weights": _atom_weights(system, c), "summed_weights": _atom_weights(system, V),
                "atoms": [a.element + str(i + 1) for i, a in enumerate(system.atoms)]}
    return c, info


def approx_pair(params, name):
    """D2: the (lecture, model) fixture names for an entry whose two presets differ, else None."""
    entry, variant = _entry(name)
    if not variant:
        return None
    names = {s.name: s.preset for s in cases.release1(params)}
    lec, mod = f"{entry} lecture", f"{entry} model"
    if lec not in names or mod not in names:
        return None
    a, b = params.presets[names[lec]], params.presets[names[mod]]
    keys = ("coupling_overrides", "edge_policy", "geometry_overrides", "hydrogen_scale")
    return None if all(a.get(k) == b.get(k) for k in keys) else (lec, mod)


def view(params, groups, name, s=1.0, level=None, conventions_sha256=None, extras=None, fragment=None,
         electrons=None, approx=None, dial=0.0, lone_pair_angle=None):
    """Everything one screen needs: result, diagram, schedule, slider state and ledger.

    Release 2 dials: electrons (D6, electron stepper), approx (D2, 0 = lecture preset .. 1 = model
    preset, linear in the couplings that differ), dial (D4, rotation within a selected degenerate set)
    and lone_pair_angle (D5, water's nonbonding pair)."""
    system = _system(params, name)
    default_electrons = system.electrons
    if electrons is not None:
        electrons = int(electrons)
        if not 0 <= electrons <= 2 * len(system.basis):
            raise ValueError(f"electrons: must be 0..{2 * len(system.basis)}")
        system = dataclasses.replace(system, electrons=electrons)
    Happ = None
    if approx is not None:
        pair = approx_pair(params, name)
        if pair is None:
            raise ValueError(f"approx: {name!r} has no lecture/model pair with different parameters")
        approx = min(max(float(approx), 0.0), 1.0)
        Hl, Hm = (build_hamiltonian(_system(params, n), params) for n in pair)
        Happ = Hl + approx * (Hm - Hl)
    result = run(system, params, groups, conventions_sha256, H=Happ)
    H = np.array(result["hamiltonian_eV"])
    sched = build_schedule(system, Group(groups[system.group]) if system.group != "C1" else None, H)
    Hs = hamiltonian_at(sched, s)
    diagram = diagram_spec(system, H, result, groups, H_now=Hs, params=params)
    frag_vectors = diagram.pop("_fragment_vectors", {})
    K = len(sched["stages"])
    lam = lambdas(s, K)
    active = next((k for k, x in enumerate(lam) if 0 < x < 1), None)
    ledger = structure = dial_info = None
    if level is not None:
        # At s < 1 the ledger and structure describe level `level` of H(s) (P3: "that moment's ledger").
        Huse = Hs if s < 1 else H
        e_s, c_s = np.linalg.eigh(Huse)
        lv_s = cluster_levels(e_s, Huse)
        ix = lv_s[min(level, len(lv_s) - 1)]
        c, dial_info = _dialled(system, c_s[:, ix], dial)
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
        c, dial_info = _dialled(system, np.array(frag_vectors[fragment]), dial)
        info = next(x for side in ("left", "right") for x in diagram.get(side, []) if x["id"] == fragment)
        # The fragment's own Hamiltonian: for a central-atom or diatomic layout, the full H restricted to the
        # fragment (cross-fragment terms vanish for a vector on one side); for the two-column layout the
        # "Atomic orbitals" are separate atoms, degenerate only under the on-site part.
        Hfrag = H if "left" in diagram["columns"] else np.diag(np.diag(H))
        ledger = contact_ledger(Hfrag, c)
        ledger.update({"level": None, "fragment": fragment, "fragment_label": info["label"], "eigenvalue_eV": info["energy_eV"],
                       "multiplicity": info["multiplicity"], "basis_labels": [f.label for f in system.basis]})
        structure = views.structure_spec(system, Hfrag, c, ledger)
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
        "electrons": {"count": system.electrons, "default": default_electrons, "neutral": views.neutral_electrons(system),
                      "charge": views.neutral_electrons(system) - system.electrons,
                      "charge_label": views.charge_label(views.neutral_electrons(system) - system.electrons),
                      "max": 2 * len(system.basis)},
        "approx": {"value": approx, "available": approx_pair(params, name) is not None},
        "dial": dial_info,
        "lone_pairs": views.lone_pairs(system, H, result, lone_pair_angle) if (lone_pair_angle is not None and name.startswith("H2O")) else None,
    }


def _quiz_context(params, groups, name):
    from . import exam
    system = _system(params, name)
    result = run(system, params, groups)
    H = np.array(result["hamiltonian_eV"])
    diagram = diagram_spec(system, H, result, groups, H_now=None, params=params)
    frags = diagram.pop("_fragment_vectors", {})
    ctx = {"result": {"levels": result["levels"], "bond_orders": result["bond_orders"]}, "diagram": diagram,
           "_eigenvectors": result["eigenvectors"]}
    return exam, system, ctx, H, frags


def quiz_cases(params):
    """Cases offered in Practice and Exam modes (Release 2, Milestone B)."""
    from .exam import QUIZ_EXCLUDE
    return [c["name"] for c in case_list(params) if c["name"] not in QUIZ_EXCLUDE]


def exam_task(params, groups, name, kind="build"):
    exam, system, ctx, H, frags = _quiz_context(params, groups, name)
    return exam.task(system, params, ctx, H, frags, kind)


def exam_check(params, groups, name, answers, kind="build"):
    exam, system, ctx, H, frags = _quiz_context(params, groups, name)
    return exam.check(system, params, ctx, H, frags, answers, kind)


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
