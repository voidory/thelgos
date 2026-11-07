"""Release 1 geometry presets (F1.1), built from unit directions per CONVENTIONS section 2.

Each builder returns a System. Bond lengths and angles come from the preset's
geometry_overrides, so params/presets.json stays the single source for them.
"""
import re

import numpy as np
from .model import Atom, BasisFunction, System

X, Y, Z = np.eye(3)


def _sp(atom, prefix=""):
    return [BasisFunction(atom, "s", None, f"{prefix}s")] + [
        BasisFunction(atom, "p", tuple(v), f"{prefix}p{a}") for v, a in zip(np.eye(3), "xyz")]


def _star(name, centre, ligand, dirs, R, preset, group, electrons, edges_hh=False):
    atoms = [Atom(centre, (0.0, 0.0, 0.0))] + [Atom(ligand, tuple(R * np.asarray(u, float))) for u in dirs]
    basis = _sp(0, centre) + [BasisFunction(i, "s", None, f"{ligand}{i}s") for i in range(1, len(atoms))]
    edges = [(0, i) for i in range(1, len(atoms))]
    return System(name, atoms, basis, edges, electrons, group=group, central_atom=0, preset=preset)


def water(params, preset="H2O_model"):
    g = params.presets[preset]["geometry_overrides"]
    t = np.deg2rad(g["angle_degrees"]) / 2
    dirs = [(0, np.sin(t), np.cos(t)), (0, -np.sin(t), np.cos(t))]
    return _star("H2O " + preset.split("_")[1], "O", "H", dirs, g["bond_length_angstrom"], preset, "C2v", 8)


def methane(params, preset="CH4_model"):
    dirs = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
    return _star("CH4", "C", "H", dirs, params.presets[preset]["geometry_overrides"]["bond_length_angstrom"], preset, "Td", 8)


def borane(params, preset="BH3_model"):
    a = np.arange(3) * 2 * np.pi / 3
    dirs = np.c_[np.cos(a), np.sin(a), np.zeros(3)]
    return _star("BH3", "B", "H", dirs, params.presets[preset]["geometry_overrides"]["bond_length_angstrom"], preset, "D3h", 6)


AXES = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


def sf6(params, preset="SF6_model"):
    R = params.presets[preset]["geometry_overrides"]["bond_length_angstrom"]
    atoms = [Atom("S", (0.0, 0.0, 0.0))] + [Atom("F", tuple(R * u)) for u in AXES]
    basis = _sp(0, "S") + [BasisFunction(i + 1, "p", tuple(-u + 0.0), f"F{i + 1}p_inward") for i, u in enumerate(AXES)]
    edges = [(a, b) for a in range(7) for b in range(a + 1, 7)]  # central-F and all F-F on (section 9.2)
    return System("SF6 " + preset.split("_")[1], atoms, basis, edges, 12, group="Oh", central_atom=0, preset=preset)


def ring_h(params, n_electrons, name, preset="H3_fixture"):
    R = params.presets[preset]["geometry_overrides"]["bond_length_angstrom"]
    a = np.arange(3) * 2 * np.pi / 3
    atoms = [Atom("H", (R / np.sqrt(3) * np.cos(t), R / np.sqrt(3) * np.sin(t), 0.0)) for t in a]
    return System(name, atoms, [BasisFunction(i, "s", None, f"H{i + 1}s") for i in range(3)],
                  [(0, 1), (0, 2), (1, 2)], n_electrons, group="D3h", preset=preset)


def h2(params, preset="H2_fixture"):
    R = params.presets[preset]["geometry_overrides"]["bond_length_angstrom"]
    atoms = [Atom("H", (0.0, 0.0, -R / 2)), Atom("H", (0.0, 0.0, R / 2))]
    return System("H2", atoms, [BasisFunction(i, "s", None, f"H{i + 1}s") for i in range(2)], [(0, 1)], 2,
                  group="D2h", linear=True, preset=preset)


def allyl(params, n_electrons=2, preset="allyl_pi"):
    g = params.presets[preset]["geometry_overrides"]
    R, t = g["bond_length_angstrom"], np.deg2rad(g["angle_degrees"]) / 2
    atoms = [Atom("C", (0.0, R * np.sin(t), R * np.cos(t))), Atom("C", (0.0, 0.0, 0.0)),
             Atom("C", (0.0, -R * np.sin(t), R * np.cos(t)))]
    names = {2: "allyl cation pi", 3: "allyl radical pi", 4: "allyl anion pi"}
    return System(names[n_electrons], atoms, [BasisFunction(i, "p", (1.0, 0.0, 0.0), f"C{i + 1}px") for i in range(3)],
                  [(0, 1), (1, 2)], n_electrons, group="C2v", preset=preset)


def six_h(params, preset="sixH_lecture"):
    R = params.presets[preset]["geometry_overrides"]["radius_angstrom"]
    atoms = [Atom("H", tuple(R * u)) for u in AXES]
    return System("six H Oh", atoms, [BasisFunction(i, "s", None, f"H{i + 1}s") for i in range(6)],
                  [(a, b) for a in range(6) for b in range(a + 1, 6)], 6, group="Oh", preset=preset)


def four_h(params, preset="fourH_diagnostic"):
    R = params.presets[preset]["geometry_overrides"]["radius_angstrom"]
    dirs = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
    atoms = [Atom("H", tuple(R * u)) for u in dirs]
    return System("four H Td", atoms, [BasisFunction(i, "s", None, f"H{i + 1}s") for i in range(4)],
                  [(a, b) for a in range(4) for b in range(a + 1, 4)], 4, group="Td", preset=preset)


DIATOMIC_ELECTRONS = {"B2": 6, "C2": 8, "N2": 10, "O2": 12, "F2": 14, "CO": 10, "Al2": 6, "Si2": 8, "P2": 10, "S2": 12, "Cl2": 14}


def diatomic(params, name):
    els = re.findall(r"[A-Z][a-z]?", name)  # "Si2" -> Si, Si; "CO" -> C, O
    a, b = (els[0], els[0]) if name.endswith("2") else (els[0], els[1])
    preset = f"{name}_model"
    R = params.presets[preset]["geometry_overrides"]["bond_length_angstrom"]
    atoms = [Atom(a, (0.0, 0.0, -R / 2)), Atom(b, (0.0, 0.0, R / 2))]
    return System(name, atoms, _sp(0, a + "1") + _sp(1, b + "2"), [(0, 1)], DIATOMIC_ELECTRONS[name],
                  group="D2h" if a == b else "C2v", linear=True, preset=preset)


# ---- CONVENTIONS v1.5 (owner sign-off of the draft, 2026-11-07) ----

def ammonia(params, preset="NH3_model"):
    """N=1 at the origin, C3 along z, lone pair towards +z; H2-4 at azimuth 0, 120, 240 degrees, each at angle b
    from -z with cos(HNH) = cos^2 b - sin^2 b / 2. H2 lies in the xz mirror plane."""
    g = params.presets[preset]["geometry_overrides"]
    c = np.cos(np.deg2rad(g["angle_degrees"]))
    b = np.arccos(np.sqrt((2 * c + 1) / 3))
    a = np.arange(3) * 2 * np.pi / 3
    dirs = np.c_[np.sin(b) * np.cos(a), np.sin(b) * np.sin(a), -np.cos(b) * np.ones(3)]
    return _star("NH3", "N", "H", dirs, g["bond_length_angstrom"], preset, "C3v", 8)


def hydrogen_sulfide(params, preset="H2S_model"):
    g = params.presets[preset]["geometry_overrides"]
    t = np.deg2rad(g["angle_degrees"]) / 2
    return _star("H2S", "S", "H", [(0, np.sin(t), np.cos(t)), (0, -np.sin(t), np.cos(t))], g["bond_length_angstrom"], preset, "C2v", 8)


def silane(params, preset="SiH4_model"):
    dirs = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
    return _star("SiH4", "Si", "H", dirs, params.presets[preset]["geometry_overrides"]["bond_length_angstrom"], preset, "Td", 8)


HYDRIDES = {"NH": ("N", 6), "OH": ("O", 7), "HF": ("F", 8)}


def hydride(params, name, variant="model"):
    """Diatomic hydride: heavy atom 1 at the origin, H 2 on +z (C2v finite subgroup of Cinf-v, as for CO)."""
    el, ne = HYDRIDES[name]
    preset = f"{name}_{variant}"
    R = params.presets[preset]["geometry_overrides"]["bond_length_angstrom"]
    atoms = [Atom(el, (0.0, 0.0, 0.0)), Atom("H", (0.0, 0.0, R))]
    return System(f"{name} {variant}", atoms, _sp(0, el + "1") + [BasisFunction(1, "s", None, "H2s")], [(0, 1)], ne,
                  group="C2v", linear=True, preset=preset)


def release1(params):
    """All fixtures (Release 1 and CONVENTIONS v1.5), in the order of checks/release1_energies.csv."""
    return [water(params, "H2O_model"), water(params, "H2O_lecture"), methane(params), borane(params),
            ammonia(params), hydrogen_sulfide(params), silane(params)] + [
        hydride(params, n, v) for n in HYDRIDES for v in ("model", "lecture")] + [
            sf6(params, "SF6_model"), sf6(params, "SF6_lecture"), h2(params),
            ring_h(params, 2, "H3+"), ring_h(params, 3, "neutral H3"), allyl(params), six_h(params), four_h(params)] + [
        diatomic(params, n) for n in ["B2", "C2", "N2", "O2", "F2", "CO", "Al2", "Si2", "P2", "S2", "Cl2"]]
