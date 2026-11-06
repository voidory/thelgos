"""Run a system end to end and return a JSON-serialisable Result (F2.9)."""
from . import __name__ as _pkg
from .hamiltonian import build_hamiltonian
from .solve import solve, cluster_levels
from .symmetry import Group, label_levels
from .analysis import analyse


def run(system, params, groups, conventions_sha256=None):
    from . import __version__
    H = build_hamiltonian(system, params)
    e, c, checks = solve(H)
    levels = cluster_levels(e, H)
    labels, linear, comm = (None, None, 0.0)
    if system.group != "C1":
        labels, linear, comm = label_levels(system, Group(groups[system.group]), H, c, levels)
    out, frontier, bo = analyse(system, H, e, c, levels)
    for k, lv in enumerate(out):
        lv["irrep"] = labels[k] if labels else "—"
        if linear:
            lv["linear_label"] = linear[k]
        lv["indices"] = levels[k]
    return {
        "engine": _pkg, "version": __version__, "conventions_sha256": conventions_sha256,
        "system": {"name": system.name, "group": system.group, "preset": system.preset,
                   "electrons": system.electrons, "counting_convention": system.counting_convention,
                   "atoms": [{"element": a.element, "position": list(a.position)} for a in system.atoms],
                   "basis": [{"atom": f.atom, "kind": f.kind, "direction": None if f.direction is None else list(f.direction),
                              "label": f.label} for f in system.basis],
                   "edges": [list(x) for x in system.edges]},
        "hamiltonian_eV": H.tolist(),
        "energies_eV": e.tolist(),
        "eigenvectors": c.tolist(),
        "levels": out,
        "frontier": frontier,
        "bond_orders": bo,
        "checks": {**checks, "commutator_max_error_eV": comm},
        "warnings": sorted({w for lv in out for w in lv["warnings"]}),
    }
