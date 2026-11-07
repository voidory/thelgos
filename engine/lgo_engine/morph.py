"""Morph paths between cases (D1) and Walsh-type correlation plots (F3.8), Release 2 Milestone D.

A path maps u in [0, 1] to a System and its Hamiltonian. Each path passes through a worked case at a
stated u, and there the frame reproduces the case exactly (tested). Couplings and geometry come from the
same parameter tables as every case. Some paths change geometry (water bending, CH4 flattening); others
switch contacts on (closing the allyl chain, assembling H3+, inserting sulfur), which in this model is
what "bringing atoms together" means.

Correlation plot (F3.8): levels are tracked by irrep wherever the path keeps a symmetry. The Hamiltonian
is block-diagonalised by the group's projectors, so levels of different irreps may cross and levels of
the same irrep never do (they are followed by rank within their block). Where the path has no symmetry,
levels are followed by maximum overlap between successive eigenvectors.
"""
import dataclasses

import numpy as np

from . import cases
from .hamiltonian import build_hamiltonian
from .symmetry import Group, projectors

TETRA = float(np.degrees(np.arccos(1 / np.sqrt(3))))   # 54.7356 deg: polar angle of a CH4 bond


def _water(params, preset, angle):
    g = params.presets[preset]["geometry_overrides"]
    t = np.deg2rad(angle) / 2
    dirs = [(0, np.sin(t), np.cos(t)), (0, -np.sin(t), np.cos(t))]
    return cases._star("H2O " + preset.split("_")[1], "O", "H", dirs, g["bond_length_angstrom"], preset, "C2v", 8)


def _ch4(params, beta):
    b = np.deg2rad(beta)
    dirs = [(np.sin(b) * np.cos(np.deg2rad(phi)), np.sin(b) * np.sin(np.deg2rad(phi)), z * np.cos(b))
            for phi, z in ((45, 1), (-45, -1), (135, -1), (-135, 1))]
    s = cases._star("CH4", "C", "H", dirs, params.presets["CH4_model"]["geometry_overrides"]["bond_length_angstrom"], "CH4_model", "D2d", 8)
    return s


def _ah3(params, centre, preset, beta, electrons):
    """AH3 with each A-H bond at angle beta (degrees) from -z; 90 is planar. Azimuths 0, 120, 240 from +x, as in
    BH3 and NH3 (CONVENTIONS section 2), so both cases sit on the same path family."""
    b = np.deg2rad(beta)
    a = np.arange(3) * 2 * np.pi / 3
    dirs = np.c_[np.sin(b) * np.cos(a), np.sin(b) * np.sin(a), -np.cos(b) * np.ones(3)]
    R = params.presets[preset]["geometry_overrides"]["bond_length_angstrom"]
    return cases._star("BH3" if centre == "B" else "NH3", centre, "H", dirs, R, preset, "C3v", electrons)


def _nh3_beta(params):
    c = np.cos(np.deg2rad(params.presets["NH3_model"]["geometry_overrides"]["angle_degrees"]))
    return float(np.degrees(np.arccos(np.sqrt((2 * c + 1) / 3))))


def _mask(system, a_set, b_set):
    """Basis-index mask for couplings between atom sets a_set and b_set."""
    N = len(system.basis)
    M = np.zeros((N, N))
    for i, f in enumerate(system.basis):
        for j, g in enumerate(system.basis):
            if (f.atom in a_set and g.atom in b_set) or (f.atom in b_set and g.atom in a_set):
                M[i, j] = 1.0
    return M


def _central_type_mask(system, kind):
    return np.diag([1.0 if f.atom == system.central_atom and f.kind == kind else 0.0 for f in system.basis])


# Each path: title, the case it reproduces and where, the parameter shown, and a builder u -> (system, H).
def _paths(params):
    p = {}
    for water in ("H2O model", "H2O lecture"):
        preset = "H2O_model" if "model" in water else "H2O_lecture"
        th0 = params.presets[preset]["geometry_overrides"]["angle_degrees"]

        def build(u, preset=preset):
            sys_ = _water(params, preset, 90 + 90 * u)
            return sys_, build_hamiltonian(sys_, params)
        p[f"{water}: bend"] = {"title": "Bend water from 90° to 180° (Walsh diagram)", "case": water, "case_u": (th0 - 90) / 90,
                               "param": ("H–O–H angle", lambda u: f"{90 + 90 * u:.1f}°"), "build": build, "group": "C2v",
                               "note": "Water is bent because the occupied 2a1 level (3a1 in textbooks that count the O 1s core) drops steeply as "
                                       "the molecule bends, since O 2pz can mix in only when it is bent, while 1b2 rises only a little. "
                                       "At 180° 2a1 and 1b1 become the degenerate πu pair."}

    def ch4(u):
        sys_ = _ch4(params, TETRA + (90 - TETRA) * u)
        return sys_, build_hamiltonian(sys_, params)
    p["CH4: flatten"] = {"title": "Flatten CH4 from tetrahedral to square planar", "case": "CH4", "case_u": 0.0,
                         "param": ("H–C–z angle", lambda u: f"{TETRA + (90 - TETRA) * u:.1f}°"), "build": ch4, "group": "D2d",
                         "note": "Labelled in D2d, the symmetry kept along the whole path. One t2 component (b2 in D2d) pairs C 2pz with a "
                                 "hydrogen combination; as the hydrogens reach the plane that pairing vanishes, leaving a nonbonding C 2pz and a "
                                 "nonbonding hydrogen combination. In this model (no H–H coupling, decision B5) the hydrogen one lies lower, at "
                                 "the H 1s energy, so it is occupied and C 2pz is the empty LUMO. Textbooks make C 2pz the occupied lone pair: "
                                 "H–H repulsion, absent here, pushes that hydrogen combination up."}

    # BH3 -> NH3 pyramidalisation (spec D1; needs C3v, CONVENTIONS v1.5). Both directions, each with its own electrons.
    if "NH3_model" in params.presets:
        b0 = _nh3_beta(params)

        def pyr(u, b0=b0):
            sys_ = cases.borane(params) if u == 0 else _ah3(params, "B", "BH3_model", 90 - (90 - b0) * u, 6)
            return sys_, build_hamiltonian(sys_, params)

        def flat(u, b0=b0):
            sys_ = cases.ammonia(params) if u == 0 else _ah3(params, "N", "NH3_model", b0 + (90 - b0) * u, 8)
            return sys_, build_hamiltonian(sys_, params)
        ang = lambda beta: float(np.degrees(np.arccos(np.clip(np.cos(np.deg2rad(beta)) ** 2 - 0.5 * np.sin(np.deg2rad(beta)) ** 2, -1, 1))))
        p["BH3: pyramidalise"] = {"title": "BH3 → NH3 shape: bend BH3 into a pyramid (Walsh diagram)", "case": "BH3", "case_u": 0.0,
                                  "param": ("H–B–H angle", lambda u, b0=b0: f"{ang(90 - (90 - b0) * u):.1f}°"), "build": pyr, "group": "C3v",
                                  "note": "Labelled in C3v, the symmetry kept along the path. Planar BH3's empty pz (a2″) becomes a1 as the "
                                          "molecule bends and mixes with the a1 bonding set: the Lewis-acid orbital drops but stays empty. "
                                          "With six electrons nothing occupied is stabilised by bending, so BH3 stays planar."}
        p["NH3: flatten"] = {"title": "Flatten NH3 (the inversion path, Walsh diagram)", "case": "NH3", "case_u": 0.0,
                             "param": ("H–N–H angle", lambda u, b0=b0: f"{ang(b0 + (90 - b0) * u):.1f}°"), "build": flat, "group": "C3v",
                             "note": "As NH3 flattens, the lone pair 3a1 loses its N 2s and H character and rises to a pure, nonbonding "
                                     "N 2pz at 120° (a2″ in D3h). With eight electrons that level is occupied, so bending lowers the total: "
                                     "NH3 is pyramidal, BH3 (six electrons, same orbitals) planar."}

    allyl = cases.allyl(params)
    Ha = build_hamiltonian(allyl, params)
    t = Ha[0, 1]
    E13 = np.zeros_like(Ha); E13[0, 2] = E13[2, 0] = t

    def close(u):
        return allyl, Ha + u * E13
    p["allyl cation pi: close the ring"] = {"title": "Close the allyl chain into a ring", "case": "allyl cation pi", "case_u": 0.0,
                                            "param": ("C1–C3 contact", lambda u: f"{100 * u:.0f} % of C1–C2"), "build": close, "group": "C2v",
                                            "note": "Closing the chain adds one contact: ε + √2 t, ε, ε − √2 t become ε + 2t, ε − t, ε − t "
                                                    "(Tutorial A becomes Tutorial B). The one-node and two-node chain orbitals become the ring's degenerate pair."}

    h3 = cases.ring_h(params, 2, "H3+")
    Hh = build_hamiltonian(h3, params)
    D = np.diag(np.diag(Hh))
    C = {k: np.zeros_like(Hh) for k in ("12", "23", "13")}
    for (i, j), k in (((0, 1), "12"), ((1, 2), "23"), ((0, 2), "13")):
        C[k][i, j] = C[k][j, i] = Hh[i, j]

    def assemble(u):
        return h3, D + C["12"] + min(1.0, 2 * u) * C["23"] + max(0.0, 2 * u - 1) * C["13"]
    p["H3+: assemble"] = {"title": "H2 + H → linear H3 → triangular H3+", "case": "H3+", "case_u": 1.0,
                          "param": ("contacts", lambda u: "H2 and a separate H" if u == 0 else ("H1–H2–H3 chain" if u <= 0.5 else "closing the triangle")),
                          "build": assemble, "group": "C1",
                          "note": "First H2–H3 switches on (a linear chain, at the midpoint), then H1–H3 (the triangle). The e′ pair is "
                                  "degenerate only when the triangle is complete. No symmetry is kept along the way, so levels are followed by overlap."}

    sf6 = cases.sf6(params)
    Hs = build_hamiltonian(sf6, params)
    lig = set(range(1, len(sf6.atoms)))
    Mc = _mask(sf6, {0}, lig)

    def insert(u):
        return sf6, Hs * (1 - Mc) + u * Hs * Mc
    p["SF6 model: insert sulfur"] = {"title": "Six-fluorine set → SF6: switch on the sulfur", "case": "SF6 model", "case_u": 1.0,
                                     "param": ("S–F coupling", lambda u: f"{100 * u:.0f} %"), "build": insert, "group": "Oh",
                                     "note": "LGOs with a sulfur partner of the same symmetry (a1g, t1u) split into bonding and antibonding "
                                             "pairs; the eg set has no partner and stays put."}

    for case, builder in (("H2O model", lambda: cases.water(params, "H2O_model")), ("CH4", lambda: cases.methane(params)),
                          ("BH3", lambda: cases.borane(params)), ("SF6 model", lambda: cases.sf6(params))):
        for kind in ("s", "p"):
            sys0 = builder()
            H0 = build_hamiltonian(sys0, params)
            M = _central_type_mask(sys0, kind)
            el = sys0.atoms[sys0.central_atom].element
            name = params.orbital_name.get((el, kind), kind)

            def sweep(u, sys0=sys0, H0=H0, M=M):
                return sys0, H0 + (-12 + 24 * u) * M
            p[f"{case}: sweep {el} {name}"] = {
                "title": f"Sweep the {el} {name} energy through the ligand levels", "case": case, "case_u": 0.5,
                "param": (f"{el} {name} shift", lambda u: f"{-12 + 24 * u:+.1f} eV".replace("-", "−")), "build": sweep, "group": sys0.group,
                "note": "Rule R2 as a continuum: mixing is strongest near resonance with an LGO of the same symmetry, and the "
                        "bonding level's character swaps from mostly central to mostly ligand across the avoided crossing."}
    return p


_CACHE = {}


def paths(params):
    key = id(params)
    if key not in _CACHE:
        _CACHE[key] = _paths(params)
    return _CACHE[key]


def path_list(params, case=None):
    return [{"id": k, "title": v["title"], "case": v["case"], "case_u": v["case_u"], "param": v["param"][0], "note": v["note"]}
            for k, v in paths(params).items() if case is None or v["case"] == case]


def frame(params, path_id, u):
    P = paths(params)[path_id]
    u = min(max(float(u), 0.0), 1.0)
    system, H = P["build"](u)
    if P["group"] == "C1" and abs(u - P["case_u"]) > 1e-12:
        system = dataclasses.replace(system, group="C1")   # off the case point the path keeps no symmetry
    return P, system, H, u


def walsh(params, groups, path_id, samples=61):
    """F3.8: level energies along the path, tracked by irrep (or by overlap when the path has no symmetry)."""
    P = paths(params)[path_id]
    us = np.linspace(0.0, 1.0, samples)
    lines = {}
    if P["group"] != "C1":
        G = Group(groups[P["group"]])
        for k, u in enumerate(us):
            system, H = P["build"](u)
            proj, _ = projectors(system, G)
            dims = {lab: d for lab, d, _ in G.irreps}
            for lab, Pl in proj.items():
                w, v = np.linalg.eigh((Pl + Pl.T) / 2)
                B = v[:, w > 0.5]
                if B.shape[1] == 0:
                    continue
                e = np.sort(np.linalg.eigvalsh(B.T @ H @ B))
                d = dims[lab]
                for r in range(len(e) // d):
                    lines.setdefault((lab, r), [None] * samples)[k] = float(np.mean(e[r * d:(r + 1) * d]))
        out = [{"irrep": lab, "rank": r, "multiplicity": dims_of(G, lab), "energies": es} for (lab, r), es in lines.items()]
    else:
        prev = None
        tracks = None
        for k, u in enumerate(us):
            system, H = P["build"](u)
            e, c = np.linalg.eigh(H)
            if prev is None:
                tracks = [[None] * samples for _ in range(len(e))]
                order = list(range(len(e)))
            else:
                ov = np.abs(prev.T @ c)
                order = [None] * len(e)
                taken = set()
                for flat in np.argsort(-ov, axis=None):
                    i, j = divmod(int(flat), len(e))
                    if order[i] is None and j not in taken:
                        order[i] = j; taken.add(j)
            for i, j in enumerate(order):
                tracks[i][k] = float(e[j])
            prev = c[:, order]
        out = [{"irrep": None, "rank": i, "multiplicity": 1, "energies": es} for i, es in enumerate(tracks)]
    out.sort(key=lambda L: L["energies"][0])
    system, _ = P["build"](0.0)
    return {"path": path_id, "u": [float(x) for x in us], "lines": out, "group": P["group"], "case_u": P["case_u"],
            "param": P["param"][0], "labels": [P["param"][1](float(x)) for x in (0.0, 0.5, 1.0)], "electrons": system.electrons}


def dims_of(G, lab):
    return next(d for l, d, _ in G.irreps if l == lab)
