"""Engine-side data for the Release 1 inspector views. Renderers only draw these.

- structure_spec (L2): atoms projected to 2D, orbital lobes coloured by phase, one contact line
  per ledger contact row (solid stabilising, dashed destabilising).
- matching_spec (P6): coupling tallies between central-atom orbitals and LGO components.
- origins (P3): for every MO level, the slider checkpoint at which it appeared.
- RULES and stage badges (L3): the instructor's rules R1-R10 and the badge each earns.
- explanation (F5.2): deterministic text built from the Result.
- provenance (F5.1): source tags for levels, labels and lines.
"""
import numpy as np

from .ledger import contact_ledger

from .solve import cluster_levels
from .schedule import hamiltonian_at

RULES = {
    "R1": ("Same colour → bonding", "Sign of each contact's contribution", "Exact, given a fixed phase convention"),
    "R2": ("Similar energy → interact", "Energy shift ~ coupling²/gap; admixed weight ~ (coupling/gap)²", "Approximation: never strictly on/off"),
    "R3": ("Close enough → interact", "Coupling magnitude falls with distance", "Qualitative: it is a parameter"),
    "R4": ("n atomic orbitals → n MOs", "Basis dimension", "Exact"),
    "R5": ("Fix one, flip the partner, combine equal-energy sets", "Stepwise symmetry adaptation", "Exact; may differ from textbook LGOs by a rotation within a degenerate set"),
    "R6": ("Count in-phase vs out-of-phase contacts", "E = ε + (2t/N)(n_in − n_out)", "Exact for ±1 coefficients and equal contacts; approximate otherwise"),
    "R7": ("Put the central atom back at its original position", "Couplings depend on direction through angular factors", "Exact"),
    "R8": ("Bond order = (bonding − antibonding)/2", "Counting definition", "Convention, not an observable"),
    "R9": ("More nodes → higher energy", "Sign changes in the coefficients", "Exact for open chains; up to degenerate pairs on rings; heuristic in 3D"),
    "R10": ("Sandwich matching by nodal planes", "|m| matching through the axis", "Exact for an ideal sandwich"),
}
STAGE_RULES = {"within_set": ["R5"], "combine": ["R5"], "between_sets": ["R7", "R2"]}


def rules_payload():
    return {k: {"rule": a, "model": b, "badge": c} for k, (a, b, c) in RULES.items()}


def _projection(pos):
    """2D view: planar molecules face the viewer, linear ones lie horizontal, 3D ones get a fixed
    oblique view. Returns the 3x3 rotation whose first two rows are screen x and screen y."""
    c = pos - pos.mean(axis=0)
    if len(pos) == 1:
        return np.eye(3)
    rank = np.linalg.matrix_rank(c, tol=1e-6)
    if rank == 1:
        axis = c[np.argmax(np.linalg.norm(c, axis=1))]
        u = axis / np.linalg.norm(axis)
        up = np.array([0, 0, 1.0]) if abs(u[2]) < 0.9 else np.array([1.0, 0, 0])
        v = up - (up @ u) * u
        v /= np.linalg.norm(v)
        return np.array([u, v, np.cross(u, v)])
    if rank == 2:
        _, _, vt = np.linalg.svd(c)
        n = vt[2]
        # Keep the principal symmetry axis vertical where the plane contains z (water, allyl).
        up = np.array([0, 0, 1.0]) if abs(n[2]) < 0.9 else np.array([0, 1.0, 0])
        v = up - (up @ n) * n
        v /= np.linalg.norm(v)
        return np.array([np.cross(v, n), v, n])
    a, b = np.deg2rad(28), np.deg2rad(-18)
    Ry = np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    Rx = np.array([[1, 0, 0], [0, np.cos(b), -np.sin(b)], [0, np.sin(b), np.cos(b)]])
    M = Rx @ Ry
    return np.array([M[0], M[2], M[1]])  # screen y = original z (up)


def structure_spec(system, H, c_vec, ledger):
    pos = system.positions()
    R = _projection(pos)
    p2 = (pos - pos.mean(axis=0)) @ R.T
    centred = pos - pos.mean(axis=0)
    atoms = [{"index": i, "element": a.element, "x": float(p2[i, 0]), "y": float(p2[i, 1]), "depth": float(p2[i, 2]),
              "xyz": [float(v) for v in centred[i]]} for i, a in enumerate(system.atoms)]
    lobes = []
    for i, f in enumerate(system.basis):
        c = float(c_vec[i])
        if abs(c) < 1e-6:
            continue
        entry = {"basis": i, "atom": f.atom, "kind": f.kind, "coefficient": c, "label": f.label}
        if f.kind == "p":
            d = np.asarray(f.direction, float) @ R.T
            entry["dir2d"] = [float(d[0]), float(d[1])]
            entry["dir3d"] = [float(v) for v in f.direction]
            entry["phase_plus_lobe"] = "positive" if c > 0 else "negative"
        else:
            entry["phase"] = "positive" if c > 0 else "negative"
        lobes.append(entry)
    contacts = []
    for k, row in enumerate(ledger["rows"]):
        if row["kind"] != "contact" or abs(row["energy_eV"]) < 1e-9:  # same cut as the ledger table
            continue
        a, b = system.basis[row["mu"]].atom, system.basis[row["nu"]].atom
        contacts.append({"ledger_row": k, "atoms": [a, b], "mu": row["mu"], "nu": row["nu"], "energy_eV": row["energy_eV"],
                         "style": "solid" if row["energy_eV"] < 0 else "dashed", "phase": row["phase"]})
    return {"atoms": atoms, "lobes": lobes, "contacts": contacts, "view_rotation": R.tolist(),
            "line_rule": "solid = stabilising, dashed = destabilising; width ∝ |energy|; red/blue reserved for orbital phase"}


def _align(Mblock):
    """Rotate degenerate LGO components so the coupling block reads like the lecture's one-to-one
    matching (orthogonal Procrustes). Returns the rotation W (columns = new components)."""
    U, _, Vt = np.linalg.svd(Mblock)
    k = min(Mblock.shape)
    return Vt.T[:, :k] @ U[:, :k].T if Mblock.shape[0] == Mblock.shape[1] else Vt.T


def matching_spec(system, H, frag_left, frag_right, P):
    """P6: tallies <central orbital|H|LGO component>, each the signed sum of its per-lobe contacts.
    Zero means no interaction. frag_* are the diagram's fragment levels with their vectors."""
    if not frag_left or not frag_right:
        return None
    N = len(system.basis)
    rows, cols, vals = [], [], []
    left = [(lv, k) for lv in frag_left for k in range(lv["multiplicity"])]
    for lv in frag_right:
        V = lv["_v_full"]
        same = [(l, k) for (l, k) in left if l["irrep"] == lv["irrep"] and l["irrep"] != "—"]
        if same and lv["multiplicity"] > 1 and len(same) == lv["multiplicity"]:
            A = np.array([l["_v_full"][:, k] for l, k in same])
            V = V @ _align(A @ H @ V)
        suffix = "xyz" if (lv["irrep"] == "—" and lv["multiplicity"] == 3) else "abc"
        for k in range(lv["multiplicity"]):
            cols.append({"label": lv["label"] + ("" if lv["multiplicity"] == 1 else suffix[k]), "irrep": lv["irrep"], "v": V[:, k]})
    for l, k in left:
        parts = l["label"].split("/")
        if "components" in l:
            name = l["components"][k]
        elif len(parts) == l["multiplicity"] and l["multiplicity"] > 1:
            name = parts[k]
        elif l["multiplicity"] == 3:
            name = f"{l['label']}{'xyz'[k]}"
        else:
            name = l["label"]
        rows.append({"label": name, "irrep": l["irrep"], "v": l["_v_full"][:, k]})
    table = [[float(r["v"] @ H @ c["v"]) for c in cols] for r in rows]
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            if abs(table[i][j]) < 1e-9:
                table[i][j] = 0.0
    return {"rows": [{"label": r["label"], "irrep": r["irrep"]} for r in rows],
            "cols": [{"label": c["label"], "irrep": c["irrep"]} for c in cols],
            "tally_eV": table,
            "rule": "Each entry is ⟨central orbital|H|LGO⟩, the signed sum of its per-lobe contacts; 0 means no interaction (P6)."}


def origins(schedule, final_levels, tol=1e-7):
    """P3: for each final level, the first checkpoint from which a level at its energy and with at
    least its multiplicity is present at every later checkpoint."""
    K = len(schedule["stages"])
    plateaus = []
    for k in range(K + 1):
        Hk = hamiltonian_at(schedule, k / K if K else 1.0)
        e = np.linalg.eigvalsh(Hk)
        plateaus.append([(float(e[ix[0]]), len(ix)) for ix in cluster_levels(e, Hk)])
    out = []
    for lv in final_levels:
        first = K
        for k in range(K, -1, -1):
            if any(abs(E - lv["energy_eV"]) < tol and m >= lv["multiplicity"] for E, m in plateaus[k]):
                first = k
            else:
                break
        label = "Present from the start (separate atoms)" if first == 0 else schedule["stages"][first - 1]["narration"]
        index = next(i for i, (E, m) in enumerate(plateaus[first]) if abs(E - lv["energy_eV"]) < tol and m >= lv["multiplicity"])
        out.append({"s": first / K if K else 1.0, "checkpoint": first, "label": label, "level_index_there": index})
    return out


def _name(level):
    ll = level.get("linear_label")
    if not ll:
        return level["irrep"]
    base, _, par = ll.partition("_")
    return {"sigma": "σ", "pi": "π"}.get(base, base) + par


def explanation(system, result, matching, frontier_levels):
    """F5.2: same config, same text. Built only from Result values."""
    lv = result["levels"]
    n = len(system.basis)
    name = system.name
    parts = [f"{name} ({system.group}, {system.preset} preset): {n} valence basis functions give {n} molecular orbitals (R4)."]
    if matching and system.central_atom is None:
        pairs = [f"{r['label']}–{c['label']}" for i, r in enumerate(matching["rows"]) for j, c in enumerate(matching["cols"])
                 if matching["tally_eV"][i][j] != 0]
        parts.append("The two atoms' orbitals couple only where symmetry allows: " + ", ".join(pairs) +
                     ". s and pz mix into σ levels; px and py form the π pairs.")
    elif matching:
        by_irrep, alone = {}, set()
        for j, col in enumerate(matching["cols"]):
            hits = [matching["rows"][i]["label"] for i in range(len(matching["rows"])) if matching["tally_eV"][i][j] != 0]
            if hits:
                by_irrep.setdefault(col["irrep"], set()).update(hits)
            else:
                alone.add(col["irrep"])
        if by_irrep:
            parts.append("Matching by symmetry: " + "; ".join(f"{irr} ligand combinations pair with {', '.join(sorted(h))}"
                                                            for irr, h in by_irrep.items()) + ".")
        alone -= set(by_irrep)
        if alone:
            parts.append(f"Ligand combinations of symmetry {', '.join(sorted(alone))} have no partner on the central atom and stay nonbonding with respect to it.")
    filled = sum(x["electrons"] for x in lv)
    parts.append(f"{filled} electrons fill the levels from the bottom (neutral-atom counting).")
    h, l = frontier_levels
    if h is not None:
        parts.append(f"HOMO: {_name(lv[h])} ({lv[h]['character']}).")
    if l is not None:
        parts.append(f"LUMO: {_name(lv[l])} ({lv[l]['character']}).")
    bo = result["bond_orders"]["counting"]
    if bo["per_bond"] is not None:
        parts.append(f"Counting bond order (R8, a convention): {bo['per_bond']:.3f} per bond over {bo['bonds']} bond{'s' if bo['bonds'] != 1 else ''}.")
    if result["warnings"]:
        parts.append("A degenerate level is partly filled, so the one-electron picture does not fix the occupation or spin.")
    return " ".join(parts)


VALENCE = {"H": 1, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7}


def neutral_electrons(system):
    """Electrons of the neutral system under the lecture's neutral-atom counting: an atom with its full
    valence basis gives its valence electrons; an atom represented only by p functions (a pi-only
    carbon, an inward-pointing ligand p) gives one electron per function."""
    n = 0
    for i, a in enumerate(system.atoms):
        kinds = [f.kind for f in system.basis if f.atom == i]
        n += VALENCE[a.element] if "s" in kinds else len(kinds)
    return n


def charge_label(q):
    return "neutral" if q == 0 else (f"{'+' if q > 0 else '−'}{abs(q) if abs(q) > 1 else ''}")


def lone_pairs(system, H, result, angle_deg):
    """D5: rotate water's two nonbonding orbitals (the highest a1 and b1 occupied levels) into a pair of
    equivalent lone pairs. Individual orbital energies change with the angle; their sum and the total
    density do not."""
    lv = result["levels"]
    occ = [k for k, x in enumerate(lv) if x["electrons"] > 0 and x["multiplicity"] == 1]
    pick = [k for k in occ if lv[k]["irrep"] in ("A1", "B1")][-2:]
    if len(pick) != 2:
        return None
    C = np.array(result["eigenvectors"])
    a, b = (C[:, lv[k]["indices"][0]] for k in pick)
    t = np.deg2rad(angle_deg)
    p1, p2 = np.cos(t) * a + np.sin(t) * b, -np.sin(t) * a + np.cos(t) * b
    rho0 = np.outer(a, a) + np.outer(b, b)
    rho = np.outer(p1, p1) + np.outer(p2, p2)
    out = []
    for v in (p1, p2):
        k = int(np.flatnonzero(abs(v) > 1e-8)[0])
        v = v if v[k] >= 0 else -v
        L = contact_ledger(H, v)
        out.append({"energy_eV": float(v @ H @ v), "structure": structure_spec(system, H, v, L)})
    return {"angle_deg": float(angle_deg), "from_levels": pick, "orbitals": out,
            "sum_eV": out[0]["energy_eV"] + out[1]["energy_eV"],
            "canonical_sum_eV": float(lv[pick[0]]["energy_eV"] + lv[pick[1]]["energy_eV"]),
            "density_change_max": float(np.max(np.abs(rho - rho0)))}


def pes_spec(system, result, curated=None):
    """Photoelectron spectrum in Koopmans' picture: one band per occupied level at -energy, height = electrons.

    Model orbital energies, not predicted ionisation energies; heights count electrons, not
    photoionisation cross-sections. Measured reference values, where curated, are returned alongside
    so the view can show the gap honestly (spec Release 1, CH4).
    """
    bands = [{"level": k, "ie_eV": -lv["energy_eV"], "electrons": lv["electrons"], "multiplicity": lv["multiplicity"]}
             for k, lv in enumerate(result["levels"]) if lv["electrons"] > 0]
    ref = ((curated or {}).get("photoelectron") or {}).get(system.name)
    caveats = [
        "Approximate. Band height counts electrons, so a triply degenerate level gives a band three times as tall. "
        "Measured intensities follow photoionisation cross-sections, so real band heights differ.",
        "Band positions are this model's orbital energies with the sign flipped (Koopmans' picture). "
        "They are not predicted ionisation energies.",
    ]
    if any(lv["warnings"] for lv in result["levels"] if lv["electrons"] > 0):
        caveats.append("A degenerate level is partly filled, so the one-electron picture does not fix which band appears first.")
    if ref:
        first = min(bands, key=lambda b: b["ie_eV"])
        caveats.append(
            f"Measured, the first band starts at {ref['first_band_adiabatic_eV']:.2f} eV and peaks at {ref['first_band_vertical_eV']:.2f} eV "
            f"(NIST CCCBDB); this model puts it at {first['ie_eV']:.1f} eV. This simple model binds bonding levels too strongly, "
            "so read the number of bands and their 3:1 heights, not the positions. "
            f"The second band lies near {ref['second_band_near_eV']:g} eV, beyond a He I lamp ({ref['he_i_eV']:.2f} eV), "
            "so a He I spectrum shows only the first band.")
    return {"bands": bands, "reference": ref, "caveats": caveats}


def provenance(system, params):
    preset = params.presets[system.preset]
    return {
        "level_energy": f"Computed: eigenvalue of H built from on-site energies (params/onsite.csv, Miessler table) and couplings (params/pairs.csv, Harrison law), preset {system.preset} ({preset['provenance_tag']}).",
        "level_label": f"Symmetry: projection onto the irreps of {system.group} (groups/{system.group}.json).",
        "fragment": "Symmetry: symmetry-adapted combination of the fragment's orbitals, diagonalised within each irrep.",
        "mixing_line": "Computed: weight of the fragment orbital in the MO (composition); shown when ≥ 5 % (convention).",
        "allowed_line": "Symmetry: same irrep on both ends.",
        "electrons": "Convention: neutral-atom electron count, filled from the bottom.",
        "geometry": f"Geometry: CONVENTIONS §2 and preset {system.preset}.",
    }
