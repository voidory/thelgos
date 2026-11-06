"""Construction-slider schedule (spec section 4, P1-P2; F2.13 pairing generator).

H(s) = H0 + sum_k lambda_k(s) H_k, where H0 holds the on-site energies and each H_k the couplings
switched on at stage k. Stages come from the lecture procedure:

1. Split the atoms into sets of symmetry-equivalent atoms (orbits of the point group).
2. Within each set, "fix one and flip its partner": pair atoms swapped by one symmetry operation,
   then pair the pairs, recursively. Each round of pairing is one stage. Its couplings are the
   edges between atoms that first become part of the same group in that round.
3. A last stage switches on every coupling between different sets ("put the central atom back").

Every enabled edge lands in exactly one stage, so the endpoint s = 1 is the full Hamiltonian
whatever the order. A stage whose couplings are all exactly zero is skipped and reported with the
note "No ligand-ligand coupling in this model" (DECISIONS.md B5).
"""
import numpy as np

SKIP_NOTE = "No ligand–ligand coupling in this model."
SKIP_NOTE_NO_CENTRE = "No coupling between these atoms in this model."


def atom_permutations(system, group):
    pos = system.positions()
    els = [a.element for a in system.atoms]
    perms = []
    for R in group.ops:
        img = (R @ pos.T).T
        perm = []
        for i in range(len(pos)):
            hit = [j for j in range(len(pos)) if els[j] == els[i] and np.linalg.norm(pos[j] - img[i]) < 1e-9]
            if len(hit) != 1:
                raise ValueError(f"{system.name}: geometry not invariant under {group.name}")
            perm.append(hit[0])
        perms.append((R, perm))
    return perms


def orbits(system, perms):
    seen, out = set(), []
    for i in range(len(system.atoms)):
        if i in seen:
            continue
        orb = sorted({p[i] for _, p in perms})
        seen.update(orb)
        out.append(orb)
    return out


def _preference(R):
    """Deterministic order for choosing the pairing operation: inversion first (lecture SF6 pairs
    each axis), then C2 about z, x, y (lecture CH4 pairs top and bottom), then any other involution."""
    I = np.eye(3)
    named = [-I, np.diag([-1, -1, 1]), np.diag([1, -1, -1]), np.diag([-1, 1, -1])]
    for k, M in enumerate(named):
        if np.allclose(R, M):
            return k
    return len(named)


def first_pairs(orbit, perms):
    """Pair atoms of one orbit with the involution that forms the most 2-cycles inside it."""
    if len(orbit) == 1:
        return [[orbit[0]]], None
    best = None
    for R, p in perms:
        if all(p[p[i]] == i for i in orbit) and any(p[i] != i for i in orbit):
            cycles = sorted({tuple(sorted((i, p[i]))) for i in orbit})
            key = (-sum(1 for c in cycles if c[0] != c[1]), _preference(R))
            if best is None or key < best[0]:
                best = (key, cycles, R)
    if best is None:
        return [[i] for i in orbit], None
    groups = [list(c) if c[0] != c[1] else [c[0]] for c in best[1]]
    return sorted(groups), best[2]


def build_schedule(system, group, H):
    """Return the stages, each with its edges, coupling matrix and narration."""
    edges = {tuple(sorted(e)) for e in system.edges}
    atom_of = [f.atom for f in system.basis]
    perms = atom_permutations(system, group) if group is not None else [(np.eye(3), list(range(len(system.atoms))))]
    orbs = orbits(system, perms)
    raw = []
    for orb in orbs:
        if len(orb) == 1:
            continue
        groups, op = first_pairs(orb, perms)
        rounds = [[g for g in groups if len(g) == 2]]
        narr = [f"Fix one, flip its partner: pair atoms {', '.join('-'.join(str(a + 1) for a in g) for g in groups if len(g) == 2)}."]
        while len(groups) > 1:
            merged, new_edges = [], []
            for k in range(0, len(groups) - 1, 2):
                merged.append(groups[k] + groups[k + 1])
                new_edges.append((groups[k], groups[k + 1]))
            if len(groups) % 2:
                merged.append(groups[-1])
            rounds.append(new_edges)
            narr.append("Combine equal-energy sets: " + "; ".join(
                f"({' '.join(str(a + 1) for a in g1)}) with ({' '.join(str(a + 1) for a in g2)})" for g1, g2 in new_edges) + ".")
            groups = merged
        for r, (spec, text) in enumerate(zip(rounds, narr)):
            if r == 0:
                stage_edges = {tuple(sorted(g)) for g in spec}
            else:
                stage_edges = {tuple(sorted((a, b))) for g1, g2 in spec for a in g1 for b in g2}
            raw.append({"kind": "within_set" if r == 0 else "combine", "set": orb, "edges": sorted(stage_edges & edges), "narration": text})
    inter = sorted(e for e in edges if not any(e[0] in o and e[1] in o for o in orbs))
    if inter:
        raw.append({"kind": "between_sets", "set": None, "edges": inter,
                    "narration": "Put the central atom back and switch on couplings between sets." if system.central_atom is not None
                    else "Switch on couplings between the sets of equivalent atoms."})
    stages, skipped = [], []
    for st in raw:
        mask = np.zeros_like(H, dtype=bool)
        for a, b in st["edges"]:
            ia = [i for i, x in enumerate(atom_of) if x == a]
            ib = [i for i, x in enumerate(atom_of) if x == b]
            mask[np.ix_(ia, ib)] = mask[np.ix_(ib, ia)] = True
        Hk = np.where(mask, H, 0.0)
        if not np.any(Hk):
            skipped.append({**st, "note": SKIP_NOTE if system.central_atom is not None else SKIP_NOTE_NO_CENTRE})
            continue
        stages.append({**st, "matrix": Hk})
    covered = sorted(e for st in stages + skipped for e in st["edges"])
    if covered != sorted(edges):
        raise AssertionError(f"{system.name}: stages do not partition the edges ({covered} vs {sorted(edges)})")
    H0 = np.diag(np.diag(H))
    if np.max(abs(H0 + sum((st["matrix"] for st in stages), np.zeros_like(H)) - H)) > 1e-12:
        raise AssertionError(f"{system.name}: stage matrices do not sum to H")
    return {"H0": H0, "stages": stages, "skipped": skipped}


def lambdas(s, K):
    """lambda_k(s): stage k ramps linearly from 0 to 1 over s in [k/K, (k+1)/K]."""
    return [float(np.clip(s * K - k, 0.0, 1.0)) for k in range(K)]


def hamiltonian_at(schedule, s):
    K = len(schedule["stages"])
    H = schedule["H0"].copy()
    for lam, st in zip(lambdas(s, K) if K else [], schedule["stages"]):
        H = H + lam * st["matrix"]
    return H


def checkpoints(schedule):
    """Slider ticks at stage boundaries (P2): each lands on a plateau, where every ramp is 0 or 1."""
    K = len(schedule["stages"])
    ticks = [{"s": 0.0, "label": "Separate atoms"}]
    for k, st in enumerate(schedule["stages"]):
        ticks.append({"s": (k + 1) / K, "label": st["narration"]})
    return ticks
