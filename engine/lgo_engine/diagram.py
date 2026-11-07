"""DiagramSpec (spec section 10, F3.1-F3.4): everything the renderer draws, computed here.

Columns: central-atom orbitals (left), MOs (centre), ligand group orbitals (right). Systems with
no central atom use two columns: atomic orbitals and MOs. Two line layers (F3.2), never mixed:
  - "allowed": dashed, a symmetry-allowed interaction (same irrep on both ends);
  - "mixing": solid, thickness proportional to the composition weight, hidden below 5%.
Renderers draw this spec and compute nothing.
"""
import numpy as np

from .solve import cluster_levels
from .symmetry import Group, projectors

MIX_MIN = 0.05
AXIS_LABEL = "Energy (model units, eV scale)"


def _fragment_levels(H, idx, P, N):
    """Symmetry-adapted fragment orbitals: within each irrep, diagonalise the fragment block on the
    range of that irrep's projector. This keeps accidental degeneracies between irreps apart (e.g. the
    H 1s combinations of CH4, all at eps_H when H-H coupling is off) and labels every level."""
    Hf = H[np.ix_(idx, idx)]
    out = []
    if P is None:
        e, c = np.linalg.eigh(Hf)
        for ix in cluster_levels(e, Hf):
            out.append((float(e[ix[0]]), len(ix), "—", c[:, ix]))
        return sorted(out, key=lambda x: x[0])
    for label, Pl in P.items():
        Pf = Pl[np.ix_(idx, idx)]
        w, v = np.linalg.eigh((Pf + Pf.T) / 2)
        B = v[:, w > 0.5]
        if B.shape[1] == 0:
            continue
        e, c = np.linalg.eigh(B.T @ Hf @ B)
        vec = B @ c
        for ix in cluster_levels(e, B.T @ Hf @ B):
            out.append((float(e[ix[0]]), len(ix), label, vec[:, ix]))
    return sorted(out, key=lambda x: (x[0], x[2]))


def diagram_spec(system, H, result, groups, H_now=None, params=None):
    """H is the full Hamiltonian (for the fragment columns); H_now the slider Hamiltonian, if any.
    Diatomics use the textbook layout: atom 1 orbitals | MOs | atom 2 orbitals."""
    N = len(system.basis)
    P = None
    if system.group != "C1":
        P, _ = projectors(system, Group(groups[system.group]))
    centre = system.central_atom
    left_idx = [i for i, f in enumerate(system.basis) if f.atom == centre] if centre is not None else []
    right_idx = [i for i, f in enumerate(system.basis) if f.atom != centre] if centre is not None else list(range(N))
    cols = {}
    diatomic = centre is None and len(system.atoms) == 2
    if diatomic:
        # Each side lists one atom's valence orbitals (s, then the p set), at its on-site energy.
        left_idx = [i for i, f in enumerate(system.basis) if f.atom == 0]
        right_idx = [i for i, f in enumerate(system.basis) if f.atom == 1]
        for side, idx, tag in (("left", left_idx, "L"), ("right", right_idx, "R")):
            el = system.atoms[system.basis[idx[0]].atom].element
            rows = []
            for kind in ("s", "p"):
                sub = [k for k, i in enumerate(idx) if system.basis[i].kind == kind]
                if not sub:
                    continue
                v = np.zeros((len(idx), len(sub)))
                v[sub, range(len(sub))] = 1.0
                name = params.orbital_name.get((el, kind), kind) if params is not None else kind
                rows.append({"energy_eV": float(H[idx[sub[0]], idx[sub[0]]]), "multiplicity": len(sub), "irrep": "—",
                             "label": f"{el} {name}", "_v": v})
            cols[side] = [{"id": f"{tag}{k}", **r} for k, r in enumerate(rows)]
    elif centre is not None:
        el = system.atoms[centre].element

        def ao_name(i):
            f = system.basis[left_idx[i]]
            name = params.orbital_name.get((el, f.kind), f.kind) if params is not None else f.kind
            if f.kind == "p":
                name += "xyz"[int(np.argmax(np.abs(f.direction)))]
            return f"{el} {name}"
        cols["left"] = []
        for k, (e, m, l, v) in enumerate(_fragment_levels(H, left_idx, P, N)):
            comps = [ao_name(int(np.argmax(abs(v[:, j])))) for j in range(m)]
            label = comps[0] if m == 1 else (comps[0][:-1] if all(c[:-1] == comps[0][:-1] for c in comps) else "/".join(comps))
            cols["left"].append({"id": f"L{k}", "energy_eV": e, "multiplicity": m, "irrep": l, "label": label, "components": comps, "_v": v})
        cols["right"] = [{"id": f"R{k}", "energy_eV": e, "multiplicity": m, "irrep": l, "label": f"LGO{k + 1}", "_v": v}
                         for k, (e, m, l, v) in enumerate(_fragment_levels(H, right_idx, P, N))]
    else:
        # No central atom: the second column is the isolated atomic orbitals, grouped by energy.
        Hd = np.diag(np.diag(H))
        e, c = np.linalg.eigh(Hd)
        cols["right"] = [{"id": f"R{k}", "energy_eV": float(e[ix[0]]), "multiplicity": len(ix), "irrep": "—",
                          "label": "/".join(sorted({system.basis[int(np.argmax(abs(c[:, i])))].label.rstrip("0123456789") for i in ix})),
                          "_v": c[:, ix]} for k, ix in enumerate(cluster_levels(e, Hd))]
    mos = []
    count = {}  # valence-only numbering within each symmetry, ascending energy (CONVENTIONS §6)
    for k, lv_ in enumerate(result["levels"]):
        key = lv_.get("linear_label") or lv_["irrep"]
        count[key] = count.get(key, 0) + 1
        mos.append({"id": f"M{k}", "number": count[key], "energy_eV": lv_["energy_eV"], "multiplicity": lv_["multiplicity"],
                    "electrons": lv_["electrons"], "irrep": lv_["irrep"], "linear_label": lv_.get("linear_label"),
                    "character": lv_["character"], "warnings": lv_["warnings"]})
    C = np.array(result["eigenvectors"])
    lines = []
    for side, idx in (("left", left_idx), ("right", right_idx)):
        for frag in cols.get(side, []):
            Fv = np.zeros((N, frag["multiplicity"]))
            Fv[idx, :] = frag["_v"]
            for k, lv_ in enumerate(result["levels"]):
                V = C[:, lv_["indices"]]
                w = float(np.sum((Fv.T @ V) ** 2) / lv_["multiplicity"])
                if frag["irrep"] != "—" and frag["irrep"] == lv_["irrep"]:
                    lines.append({"layer": "allowed", "from": frag["id"], "to": f"M{k}", "style": "dashed"})
                if w >= MIX_MIN:
                    lines.append({"layer": "mixing", "from": frag["id"], "to": f"M{k}", "style": "solid", "weight": w})
    for side, idx in (("left", left_idx), ("right", right_idx)):
        for frag in cols.get(side, []):
            full = np.zeros((N, frag["multiplicity"]))
            full[idx, :] = frag["_v"]
            frag["_v_full"] = full
    from .views import matching_spec
    matching = matching_spec(system, H, cols.get("left", []), cols.get("right", []), P) if (centre is not None or diatomic) else None
    vectors = {}
    for side in cols:
        for frag in cols[side]:
            frag.pop("_v")
            vectors[frag["id"]] = frag.pop("_v_full")
    if diatomic:
        titles = {"left": f"{system.atoms[0].element} (atom 1)", "centre": "Molecular orbitals", "right": f"{system.atoms[1].element} (atom 2)"}
    elif centre is not None:
        titles = {"left": f"{system.atoms[centre].element} orbitals", "centre": "Molecular orbitals", "right": "Ligand group orbitals"}
    else:
        titles = {"right": "Atomic orbitals", "centre": "Molecular orbitals"}
    spec = {"columns": (["left", "centre", "right"] if (centre is not None or diatomic) else ["right", "centre"]),
            "column_titles": titles,
            "left": cols.get("left", []), "right": cols["right"], "centre": mos, "lines": lines,
            "axis_label": AXIS_LABEL, "electrons": system.electrons, "matching": matching}
    spec["_fragment_vectors"] = vectors  # internal: app_api uses it for a selected fragment orbital, then drops it
    if H_now is not None:
        e = np.linalg.eigvalsh(H_now)
        spec["slider_levels"] = [{"energy_eV": float(e[ix[0]]), "multiplicity": len(ix)} for ix in cluster_levels(e, H_now)]
    return spec
