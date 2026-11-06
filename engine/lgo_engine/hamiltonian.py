"""Hamiltonian assembly (F2.1, F2.2) with s/p Slater-Koster factors.

Signs follow the authoritative +z table in CONVENTIONS section 1: for unit vector n from a to b,
  <s_a|H|s_b> = V_sss
  <s_a|H|p_b(d)> = +(d.n) V_sps        <p_a(d)|H|s_b> = -(d.n) V_sps
  <p_a(d1)|H|p_b(d2)> = (d1.n)(d2.n)(V_pps - V_ppp) + (d1.d2) V_ppp
p functions are arbitrary real directions d . (px, py, pz), which covers Cartesian p, inward
sigma donors and local pi frames with one formula. Each a<b block is computed once and its
transpose fills the reverse block.
"""
import numpy as np

CHANNELS = {("s", "s"): ("ss_sigma",), ("s", "p"): ("sp_sigma",), ("p", "s"): ("sp_sigma",),
            ("p", "p"): ("pp_sigma", "pp_pi")}


def element(fa, fb, n, V):
    """Matrix element between basis functions fa (on atom a) and fb (on atom b)."""
    if fa.kind == "s" and fb.kind == "s":
        return V["ss_sigma"]
    if fa.kind == "s":
        return float(np.dot(fb.direction, n)) * V["sp_sigma"]
    if fb.kind == "s":
        return -float(np.dot(fa.direction, n)) * V["sp_sigma"]
    d1, d2 = np.asarray(fa.direction, float), np.asarray(fb.direction, float)
    return float((d1 @ n) * (d2 @ n) * (V["pp_sigma"] - V["pp_pi"]) + (d1 @ d2) * V["pp_pi"])


def build_hamiltonian(system, params):
    system.validate()
    pos = system.positions()
    N = len(system.basis)
    H = np.zeros((N, N))
    els = [a.element for a in system.atoms]
    for i, f in enumerate(system.basis):
        H[i, i] = params.onsite_energy(els[f.atom], f.kind)
    by_atom = {}
    for i, f in enumerate(system.basis):
        by_atom.setdefault(f.atom, []).append(i)
    for a, b in system.edges:
        if a > b:
            a, b = b, a
        r = pos[b] - pos[a]
        d = float(np.linalg.norm(r))
        n = r / d
        needed = {ch for i in by_atom.get(a, []) for j in by_atom.get(b, [])
                  for ch in CHANNELS[(system.basis[i].kind, system.basis[j].kind)]}
        V = {ch: params.coupling(els[a], els[b], ch, d, system.preset) for ch in needed}
        for i in by_atom.get(a, []):
            for j in by_atom.get(b, []):
                fa, fb = system.basis[i], system.basis[j]
                ov = params.override(system.preset, els[a], fa.kind, els[b], fb.kind)
                H[i, j] = H[j, i] = ov if ov is not None else element(fa, fb, n, V)
    if np.max(abs(H - H.T)) > 1e-12:
        raise AssertionError("Hamiltonian not symmetric")
    return H
