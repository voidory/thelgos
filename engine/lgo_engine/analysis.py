"""Composition, bonding character, filling, frontier orbitals, bond orders (F2.6-F2.8, F2.10, F2.11)."""
import numpy as np

NONBONDING_TOL_EV = 1e-6


def warn_shell(d, n):
    """CONVENTIONS section 7: any degenerate level with 0 < n < 2d warns."""
    if d > 1 and 0 < n < 2 * d:
        return ["occupation_not_unique", "Jahn_Teller_warning", "spin_state_not_determined"]
    return []


def bonding_scope(system):
    """Atom pairs whose contacts define bonding character: central-ligand pairs when there is a
    central atom, otherwise every enabled edge."""
    if system.central_atom is not None:
        return [e for e in system.edges if system.central_atom in e]
    return list(system.edges)


def analyse(system, H, e, c, levels):
    basis = system.basis
    N = len(basis)
    scope = {tuple(sorted(p)) for p in bonding_scope(system)}
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)
             if tuple(sorted((basis[i].atom, basis[j].atom))) in scope]
    left = system.electrons
    out = []
    for ix in levels:
        V = c[:, ix]
        d = len(ix)
        n = min(left, 2 * d)
        left -= n
        bond = sum(2 * float(V[i] @ V[j]) * H[i, j] for i, j in pairs) / d
        character = "bonding" if bond < -NONBONDING_TOL_EV else ("antibonding" if bond > NONBONDING_TOL_EV else "nonbonding")
        comp_ao = (V ** 2).sum(axis=1) / d
        comp_atom = {}
        for i, f in enumerate(basis):
            comp_atom[f.atom] = comp_atom.get(f.atom, 0.0) + float(comp_ao[i])
        out.append({"energy_eV": float(e[ix[0]]), "multiplicity": d, "electrons": n,
                    "warnings": warn_shell(d, n), "bond_energy_eV": float(bond), "character": character,
                    "composition_ao": [float(x) for x in comp_ao],
                    "composition_atom": {str(k): v for k, v in sorted(comp_atom.items())}})
    if left:
        raise AssertionError("electron count not conserved")
    occ = [k for k, lv in enumerate(out) if lv["electrons"]]
    vac = [k for k, lv in enumerate(out) if lv["electrons"] < 2 * lv["multiplicity"]]
    frontier = {"HOMO": occ[-1] if occ else None, "LUMO": vac[0] if vac else None}
    for key in ("HOMO", "LUMO"):
        k = frontier[key]
        if k is not None:
            frontier[key + "_flags"] = [x for x, on in [("degenerate", out[k]["multiplicity"] > 1),
                                                         ("nonbonding", out[k]["character"] == "nonbonding")] if on]
    return out, frontier, bond_orders(system, c, levels, out, scope)


def bond_orders(system, c, levels, out, scope):
    """(i) Counting: (bonding - antibonding electrons)/2 shared over the bonding-scope edges, which
    are symmetry-equivalent in every Release 1 case. (ii) Density matrix: Coulson p_ab for one
    orbital per atom, Wiberg W_AB otherwise. The two are reported side by side, never forced to agree."""
    nb = sum(lv["electrons"] for lv in out if lv["character"] == "bonding")
    na = sum(lv["electrons"] for lv in out if lv["character"] == "antibonding")
    edges = sorted(scope)
    counting = (nb - na) / 2 / len(edges) if edges else None
    occ = np.zeros(c.shape[1])
    for ix, lv in zip(levels, out):
        occ[ix] = lv["electrons"] / len(ix)
    P = (c * occ) @ c.T
    one_per_atom = all(sum(1 for f in system.basis if f.atom == a) == 1 for a in range(len(system.atoms)))
    dm = {}
    for a, b in edges:
        ia = [i for i, f in enumerate(system.basis) if f.atom == a]
        jb = [j for j, f in enumerate(system.basis) if f.atom == b]
        dm[f"{a}-{b}"] = float(P[ia[0], jb[0]]) if one_per_atom else float(sum(P[i, j] ** 2 for i in ia for j in jb))
    return {"counting": {"definition": "(bonding - antibonding electrons)/2 over equivalent bonds; nonbonding count zero",
                         "per_bond": counting, "bonding_electrons": nb, "antibonding_electrons": na, "bonds": len(edges)},
            "density_matrix": {"definition": "Coulson p_ab" if one_per_atom else "Wiberg W_AB", "per_bond": dm}}
