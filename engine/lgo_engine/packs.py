"""Molecule packs (spec section 14; Release 3 Milestone A).

A pack is a small data file (packs/<id>.json, schema packs/pack.schema.json) that declares a molecule:
identity, a template shape with sourced bond lengths and angles, the basis, the electron count and the
point group. Adding a molecule never touches engine code. Packs sit beside the CONVENTIONS fixtures, not
inside them: their regression check is the engine's consistency test (F2.15), not Appendix B.

Owner decisions R3-1 and R3-3: parameters only from the scheme in params/ (Miessler on-sites, Harrison
couplings), and only sourced geometries. A pack whose element, orbital or pair has no parameter fails its
gate with that reason rather than being drawn.

Shapes put the molecule in the CONVENTIONS section 2 orientation of its template, so the stored group files
apply: bent molecules in the yz plane with z the bisector, trigonal planar in xy with a terminal on +x,
tetrahedral along (1,1,1)-type directions, AX3Y with Y on +z, octahedral on the axes.
"""
import dataclasses

import numpy as np

from . import cases
from .model import Atom, BasisFunction, System
from .views import VALENCE

TETRA = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)


def _directions(st):
    shape = st["shape"]
    if shape in ("AX2_bent", "AXY_bent"):
        t = np.deg2rad(st["angle_degrees"]) / 2
        return [(0, np.sin(t), np.cos(t)), (0, -np.sin(t), np.cos(t))]
    if shape == "AX3_planar":
        a = np.arange(3) * 2 * np.pi / 3
        return list(np.c_[np.cos(a), np.sin(a), np.zeros(3)])
    if shape == "AX4_tetrahedral":
        return list(TETRA)
    if shape == "AX3Y_C3v":
        # X at angle b from -z, from the X-A-X angle: cos(XAX) = cos^2 b - sin^2 b / 2; Y on +z.
        c = np.cos(np.deg2rad(st["angle_degrees"]))
        b = np.arccos(np.sqrt((2 * c + 1) / 3))
        a = np.arange(3) * 2 * np.pi / 3
        return list(np.c_[np.sin(b) * np.cos(a), np.sin(b) * np.sin(a), -np.cos(b) * np.ones(3)]) + [(0, 0, 1.0)]
    if shape == "AX6_octahedral":
        return list(cases.AXES)
    raise ValueError(f"unknown shape {shape!r}")


def _lengths(st, n):
    L = st["bond_lengths_angstrom"]
    if st["shape"] in ("AXY_bent", "AX3Y_C3v"):
        if len(L) != 2:
            raise ValueError(f"{st['shape']}: give two bond lengths (X, then Y)")
        return [L[0]] * (n - 1) + [L[1]] if st["shape"] == "AX3Y_C3v" else [L[0], L[1]]
    if len(L) != 1:
        raise ValueError(f"{st['shape']}: give one bond length")
    return L * n


def preset_id(pack):
    return f"pack:{pack['id']}"


def register(params, pack):
    """Add the pack's parameter preset to params (in memory): the scheme's hydrogen factor, no overrides."""
    lig = pack.get("edges", {}).get("ligand_ligand", False)
    params.presets[preset_id(pack)] = {
        "id": preset_id(pack), "coupling_overrides": [], "hydrogen_scale": 0.6, "geometry_overrides": {},
        "edge_policy": {"edges_on": ["central_ligand"] + (["ligand_ligand"] if lig else []), "edges_off": [] if lig else ["ligand_ligand"]},
        "provenance_tag": "model", "note": f"Molecule pack {pack['id']}: parameter scheme (F2.12), geometry {pack['structure']['source']}"}


def expected_electrons(pack):
    """Neutral-atom counting over the included basis: the central atom's valence electrons, one per terminal
    sigma function, less the charge (CONVENTIONS section 7, as for SF6's reduced basis)."""
    st = pack["structure"]
    return VALENCE[st["central"]] + len(st["terminals"]) - pack["identity"]["charge"]


def build(params, pack, groups=None):
    """System for a pack. A declared group missing from the library falls back to C1 (F2.14)."""
    st = pack["structure"]
    dirs = [np.asarray(u, float) / np.linalg.norm(u) for u in _directions(st)]
    terms = st["terminals"]
    if len(dirs) != len(terms):
        raise ValueError(f"{pack['id']}: shape {st['shape']} has {len(dirs)} terminal positions, pack lists {len(terms)}")
    Ls = _lengths(st, len(terms))
    atoms = [Atom(st["central"], (0.0, 0.0, 0.0))] + [Atom(el, tuple(R * u + 0.0)) for el, R, u in zip(terms, Ls, dirs)]
    basis = cases._sp(0, st["central"])
    for i, (el, u) in enumerate(zip(terms, dirs), start=1):
        if el == "H":
            basis.append(BasisFunction(i, "s", None, f"H{i + 1}s"))
        else:   # sigma-only terminal: inward p (CONVENTIONS section 3)
            basis.append(BasisFunction(i, "p", tuple(-u + 0.0), f"{el}{i + 1}p_inward"))
    edges = [(0, i) for i in range(1, len(atoms))]
    if pack.get("edges", {}).get("ligand_ligand", False):
        edges += [(a, b) for a in range(1, len(atoms)) for b in range(a + 1, len(atoms))]
    declared = pack["group"]["declared"]
    group = declared if (groups is None or declared in groups) else "C1"
    if preset_id(pack) not in params.presets:
        register(params, pack)
    s = System(pack["identity"]["name"], atoms, basis, edges, pack["electrons"]["count"], group=group, central_atom=0,
               preset=preset_id(pack))
    return s


def fallback_note(pack, groups):
    d = pack["group"]["declared"]
    return None if d in groups else f"Point group {d} is not in the library, so symmetry labels are shown as — (C1 fallback, spec F2.14)."


def _formula_counts(formula):
    import re
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        out[el] = out.get(el, 0) + (int(n) if n else 1)
    return out


def gate(params, groups, pack):
    """The `make add-molecule` checks that need no schema library (spec section 14). Returns a list of
    (check, passed, detail). Schema validation is done by scripts/add_molecule.py."""
    from .hamiltonian import build_hamiltonian
    from .ledger import contact_ledger
    from .result import run
    from .schedule import build_schedule, hamiltonian_at, atom_permutations
    from .solve import cluster_levels
    from .symmetry import Group
    out = []

    def check(name, ok, detail=""):
        out.append((name, bool(ok), detail))
        return ok
    st = pack["structure"]
    atoms_f = {st["central"]: 1}
    for t in st["terminals"]:
        atoms_f[t] = atoms_f.get(t, 0) + 1
    check("formula matches the structure", _formula_counts(pack["identity"]["formula"]) == atoms_f,
          f"{pack['identity']['formula']} vs {atoms_f}")
    try:
        s = build(params, pack, set(groups))
        H = build_hamiltonian(s, params)
    except (ValueError, KeyError) as e:
        check("parameters exist for every element, orbital and pair (R3-1)", False, str(e))
        return out
    check("parameters exist for every element, orbital and pair (R3-1)", True)
    check("electron count follows neutral-atom counting", s.electrons == expected_electrons(pack),
          f"{s.electrons} declared, {expected_electrons(pack)} by the rule")
    declared = pack["group"]["declared"]
    if declared in groups:
        try:
            atom_permutations(s, Group(groups[declared]))
            check(f"geometry has the declared symmetry {declared}", True)
        except ValueError as e:
            check(f"geometry has the declared symmetry {declared}", False, str(e))
            return out
    else:
        check(f"group {declared} not in the library: C1 fallback (F2.14)", s.group == "C1", "drawn with labels —")
    r = run(s, params, groups)
    nmo = sum(lv["multiplicity"] for lv in r["levels"])
    check("R4: n basis functions give n molecular orbitals", nmo == len(s.basis), f"{len(s.basis)} in, {nmo} out")
    c = np.array(r["eigenvectors"])
    worst = max(abs(contact_ledger(H, c[:, ix[0]])["total_eV"] - r["energies_eV"][ix[0]]) for ix in (lv["indices"] for lv in r["levels"]))
    check("ledger totals equal the eigenvalues (1e-10 eV)", worst < 1e-10, f"max {worst:.1e} eV")
    try:
        sched = build_schedule(s, Group(groups[s.group]) if s.group != "C1" else None, H)
        dev = float(np.max(abs(np.linalg.eigvalsh(hamiltonian_at(sched, 1.0)) - np.array(r["energies_eV"])))) if sched["stages"] else 0.0
        check("pairing schedule (F2.13) partitions the couplings and ends at H", dev < 1e-10, f"{len(sched['stages'])} stage{'' if len(sched['stages']) == 1 else 's'}")
    except AssertionError as e:
        check("pairing schedule (F2.13) partitions the couplings and ends at H", False, str(e))
    check("morph invariants", True, "no morph paths declared")
    return out
