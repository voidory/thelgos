"""Entry points the browser worker calls. Pure functions returning JSON-ready dicts."""
import numpy as np

from . import cases
from .hamiltonian import build_hamiltonian
from .result import run
from .schedule import build_schedule, hamiltonian_at, checkpoints, lambdas
from .symmetry import Group
from .diagram import diagram_spec
from .ledger import contact_ledger


def case_list(params):
    return [{"name": s.name, "preset": s.preset, "provenance_tag": params.presets[s.preset]["provenance_tag"],
             "group": s.group, "electrons": s.electrons} for s in cases.release1(params)]


def _system(params, name):
    for s in cases.release1(params):
        if s.name == name:
            return s
    raise ValueError(f"unknown case {name!r}")


def view(params, groups, name, s=1.0, level=None, conventions_sha256=None):
    """Everything one screen needs: result, diagram, schedule, slider state and ledger."""
    system = _system(params, name)
    result = run(system, params, groups, conventions_sha256)
    H = np.array(result["hamiltonian_eV"])
    sched = build_schedule(system, Group(groups[system.group]) if system.group != "C1" else None, H)
    Hs = hamiltonian_at(sched, s)
    diagram = diagram_spec(system, H, result, groups, H_now=Hs)
    K = len(sched["stages"])
    lam = lambdas(s, K)
    active = next((k for k, x in enumerate(lam) if 0 < x < 1), None)
    ledger = None
    if level is not None:
        lv = result["levels"][level]
        c = np.array(result["eigenvectors"])[:, lv["indices"][0]]
        ledger = contact_ledger(H, c)
        ledger["level"] = level
        ledger["eigenvalue_eV"] = lv["energy_eV"]
        ledger["basis_labels"] = [f.label for f in system.basis]
    return {
        "case": name,
        "preset": {"id": system.preset, "provenance_tag": params.presets[system.preset]["provenance_tag"]},
        "result": {k: result[k] for k in ("levels", "frontier", "bond_orders", "warnings", "checks", "version", "conventions_sha256")},
        "diagram": diagram,
        "schedule": {"stages": [{"kind": st["kind"], "edges": st["edges"], "narration": st["narration"]} for st in sched["stages"]],
                     "skipped": [{"kind": st["kind"], "edges": st["edges"], "narration": st["narration"], "note": st["note"]} for st in sched["skipped"]],
                     "checkpoints": checkpoints(sched), "s": s, "lambdas": lam, "active_stage": active},
        "ledger": ledger,
    }
