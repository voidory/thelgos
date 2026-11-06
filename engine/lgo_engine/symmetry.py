"""Generic symmetry analysis from group data (F2.4, F2.5).

Groups are passed in as the dicts stored in groups/*.json (operations, irreps with characters).
The basis representation is built from atom permutations and the p transformation; levels are
labelled by projecting each degenerate eigenspace. For linear molecules run in a finite subgroup
a pi pair spans two irreps; such sets are labelled by their exact decomposition and mapped to
Dinfh/Cinfv names (K14 in CHECKER-CHANGES.md).
"""
import numpy as np

LINEAR_LABELS = {"D2h": {"Ag": "sigma_g", "B1u": "sigma_u", "B2u+B3u": "pi_u", "B2g+B3g": "pi_g"},
                 "C2v": {"A1": "sigma", "B1+B2": "pi"}}  # s/p basis only; delta needs d functions


class Group:
    def __init__(self, data):
        self.name = data["group"]
        self.ops = [np.array(op["matrix"], float) for op in data["operations"]]
        self.irreps = [(r["label"], int(r["dimension"]), np.array(r["characters"], float)) for r in data["irreps"]]


def basis_representation(system, group):
    pos = system.positions()
    els = [a.element for a in system.atoms]
    N = len(system.basis)
    mats = []
    for R in group.ops:
        D = np.zeros((N, N))
        for j, f in enumerate(system.basis):
            img = R @ pos[f.atom]
            dst = [b for b in range(len(pos)) if els[b] == els[f.atom] and np.linalg.norm(pos[b] - img) < 1e-9]
            if len(dst) != 1:
                raise ValueError(f"{system.name}: geometry is not invariant under the {group.name} operations")
            for i, g in enumerate(system.basis):
                if g.atom != dst[0] or g.kind != f.kind:
                    continue
                D[i, j] = 1.0 if f.kind == "s" else float(np.asarray(g.direction) @ R @ np.asarray(f.direction))
        if np.max(abs(D.T @ D - np.eye(N))) > 1e-10:
            raise ValueError(f"{system.name}: basis is not closed under {group.name} (check p directions)")
        mats.append(D)
    return np.array(mats)


def projectors(system, group):
    D = basis_representation(system, group)
    order = len(group.ops)
    return {label: dim / order * np.einsum("g,gij->ij", chars, D) for label, dim, chars in group.irreps}, D


def label_levels(system, group, H, c, levels):
    P, D = projectors(system, group)
    comm = float(max(np.max(abs(H @ g - g @ H)) for g in D))
    if comm > 1e-10:
        raise AssertionError(f"{system.name}: H does not commute with {group.name} (max {comm})")
    labels = []
    for ix in levels:
        V = c[:, ix]
        n = {l: float(np.trace(V.T @ p @ V)) for l, p in P.items()}
        best = max(n, key=n.get)
        if n[best] / len(ix) > 0.999:
            labels.append(best)
            continue
        if not system.linear:
            raise AssertionError(f"{system.name}: level is not a clean irrep: {n}")
        if any(abs(x - round(x)) > 1e-3 for x in n.values()) or abs(sum(n.values()) - len(ix)) > 1e-6:
            raise AssertionError(f"{system.name}: non-integer irrep decomposition {n}")
        labels.append("+".join(sorted(l for l, x in n.items() if x > 0.5)))
    linear = [LINEAR_LABELS.get(group.name, {}).get(l) for l in labels] if system.linear else None
    return labels, linear, comm
