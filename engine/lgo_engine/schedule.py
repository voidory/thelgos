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


_PERM_CACHE = {}
_OPTIONS_CACHE = {}


def atom_permutations(system, group):
    pos = system.positions()
    key = (group.name if hasattr(group, "name") else id(group), pos.round(12).tobytes(), tuple(a.element for a in system.atoms))
    if key in _PERM_CACHE:
        return _PERM_CACHE[key]
    _PERM_CACHE[key] = _atom_permutations(system, group, pos)
    return _PERM_CACHE[key]


def _atom_permutations(system, group, pos):
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


_NAMED = [-np.eye(3), np.diag([-1.0, -1, 1]), np.diag([1.0, -1, -1]), np.diag([-1.0, 1, -1])]


def _preference(R):
    """Deterministic order for choosing the pairing operation: inversion first (lecture SF6 pairs
    each axis), then C2 about z, x, y (lecture CH4 pairs top and bottom), then any other involution."""
    for k, M in enumerate(_NAMED):
        if float(np.max(np.abs(R - M))) < 1e-8:
            return k
    return len(_NAMED)


def _involution_pairs(orbit, p):
    """2-cycles of permutation p inside orbit, or None if p is not an involution on it that pairs anything."""
    if not all(p[p[i]] == i for i in orbit) or not any(p[i] != i for i in orbit):
        return None
    return sorted({tuple(sorted((i, p[i]))) for i in orbit})


_CAND_CACHE = {}


def _candidates(orbit, perms):
    """Involutions that pair the most atoms of this orbit, best (lecture-preferred) first: [(perm, cycles, R)]."""
    key = (tuple(tuple(p) for _, p in perms), tuple(orbit))
    if key not in _CAND_CACHE:
        found = []
        for R, p in perms:
            cycles = _involution_pairs(orbit, p)
            if cycles is not None:
                found.append((-sum(1 for c in cycles if c[0] != c[1]), _preference(R), p, cycles, R))
        found.sort(key=lambda x: (x[0], x[1]))
        best = found[0][0] if found else None
        _CAND_CACHE[key] = [(p, cycles, R) for n, _, p, cycles, R in found if n == best]
    return _CAND_CACHE[key]


def first_pairs(orbit, perms, forced=None):
    """Pair atoms of one orbit with the involution that forms the most 2-cycles inside it. With
    `forced` (a permutation), use that involution instead wherever it forms as many 2-cycles (P5)."""
    if len(orbit) == 1:
        return [[orbit[0]]], None
    cand = _candidates(orbit, perms)
    if not cand:
        return [[i] for i in orbit], None
    p, cycles, R = cand[0]
    if forced is not None:
        fc = _involution_pairs(orbit, forced)
        if fc is not None and sum(1 for c in fc if c[0] != c[1]) == sum(1 for c in cycles if c[0] != c[1]):
            cycles, R = fc, None
    groups = [list(c) if c[0] != c[1] else [c[0]] for c in cycles]
    return sorted(groups), R


VARIANT_TITLES = {
    "lecture": "Lecture order: pair, combine, then the central atom",
    "central_first": "Central atom first, then pair and combine the ligands",
    "one_step": "Every coupling at once",
    "atom_by_atom": "Atom by atom, in numbering order",
}
NAMES = {"within_set": "Fix one, flip partner", "combine": "Combine equal-energy sets", "between_sets": "Put the central atom back",
         "all_at_once": "Switch everything on", "add_atom": "Add an atom"}


def _symmetry_stages(system, perms, orbs, edges, forced=None):
    raw = []
    for orb in orbs:
        if len(orb) == 1:
            continue
        groups, op = first_pairs(orb, perms, forced)
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
    tail = []
    if inter:
        tail.append({"kind": "between_sets", "set": None, "edges": inter,
                     "narration": "Put the central atom back and switch on couplings between sets." if system.central_atom is not None
                     else "Switch on couplings between the sets of equivalent atoms."})
    return raw, tail


def _pairing_ops(system, perms, orbs):
    """P5: the distinct pairings (one per involution) that pair as many atoms as the lecture's choice in every set."""
    out, seen = [], set()
    # Only involutions that are a best pairing for some set can change the schedule.
    pool, ids = [], set()
    for orb in orbs:
        if len(orb) > 1:
            for p, _, _ in _candidates(orb, perms):
                if tuple(p) not in ids:
                    ids.add(tuple(p)); pool.append((None, p))
    for R, p in pool:
        sig = []
        for orb in orbs:
            if len(orb) == 1:
                continue
            g, _ = first_pairs(orb, perms, p)
            sig.append(tuple(map(tuple, g)))
        sig = tuple(sig)
        if not sig:
            continue
        # Pairings related by a symmetry operation give the same plateaus: keep one per class.
        cls = min(tuple(tuple(sorted(tuple(sorted(q[a] for a in g)) for g in orb_sig)) for orb_sig in sig) for _, q in perms)
        if cls not in seen:
            seen.add(cls)
            out.append((p, sig))
    return out


def build_schedule(system, group, H, variant="lecture"):
    """Return the stages, each with its edges, coupling matrix and narration.

    variant (P5): "lecture" (default), "pair:k" (the k-th alternative pairing), "central_first",
    "one_step" or "atom_by_atom". Every variant partitions the same couplings, so H(1) = H."""
    edges = {tuple(sorted(e)) for e in system.edges}
    atom_of = [f.atom for f in system.basis]
    perms = atom_permutations(system, group) if group is not None else [(np.eye(3), list(range(len(system.atoms))))]
    orbs = orbits(system, perms)
    if variant == "lecture":
        lig, tail = _symmetry_stages(system, perms, orbs, edges)
        raw = lig + tail
    elif variant.startswith("pair:"):
        ops = _pairing_ops(system, perms, orbs)
        k = int(variant.split(":")[1])
        if not 0 <= k < len(ops):
            raise ValueError(f"schedule: no pairing {k} for {system.name}")
        lig, tail = _symmetry_stages(system, perms, orbs, edges, forced=ops[k][0])
        raw = lig + tail
    elif variant == "central_first":
        lig, tail = _symmetry_stages(system, perms, orbs, edges)
        raw = tail + lig
    elif variant == "one_step":
        raw = [{"kind": "all_at_once", "set": None, "edges": sorted(edges), "narration": "Switch on every coupling at once."}]
    elif variant == "atom_by_atom":
        raw = []
        for a in range(1, len(system.atoms)):
            es = sorted(e for e in edges if max(e) == a)
            raw.append({"kind": "add_atom", "set": [a], "edges": es,
                        "narration": f"Add atom {a + 1} ({system.atoms[a].element}) and its couplings to atoms already present."})
    else:
        raise ValueError(f"schedule: unknown variant {variant!r}")
    stages, skipped = [], []
    for st in raw:
        mask = np.zeros_like(H, dtype=bool)
        for a, b in st["edges"]:
            ia = [i for i, x in enumerate(atom_of) if x == a]
            ib = [i for i, x in enumerate(atom_of) if x == b]
            mask[np.ix_(ia, ib)] = mask[np.ix_(ib, ia)] = True
        Hk = np.where(mask, H, 0.0)
        if not np.any(Hk):
            note = (SKIP_NOTE if system.central_atom is not None else SKIP_NOTE_NO_CENTRE) if st["kind"] in ("within_set", "combine") \
                else "No couplings switch on at this step in this model."
            skipped.append({**st, "note": note})
            continue
        stages.append({**st, "matrix": Hk})
    covered = sorted(e for st in stages + skipped for e in st["edges"])
    if covered != sorted(edges):
        raise AssertionError(f"{system.name}: stages do not partition the edges ({covered} vs {sorted(edges)})")
    H0 = np.diag(np.diag(H))
    if np.max(abs(H0 + sum((st["matrix"] for st in stages), np.zeros_like(H)) - H)) > 1e-12:
        raise AssertionError(f"{system.name}: stage matrices do not sum to H")
    return {"H0": H0, "stages": stages, "skipped": skipped, "variant": variant}


def _signature(sched):
    return tuple(tuple(map(tuple, st["edges"])) for st in sched["stages"])


def schedule_options(system, group, H):
    """P5: the schedules on offer for a case, with the ones that give the same stages as an earlier option removed
    (with ligand-ligand couplings off, CH4's other orders collapse onto the lecture one)."""
    # The options depend only on geometry, symmetry, the edges and which couplings are non-zero: cache on those
    # (the view is recomputed on every slider move, D10).
    key = (system.name, system.group, system.positions().round(12).tobytes(), tuple(sorted(map(tuple, system.edges))),
           (abs(H) > 0).tobytes())
    if key in _OPTIONS_CACHE:
        return [dict(o) for o in _OPTIONS_CACHE[key]]
    out = _schedule_options(system, group, H)
    _OPTIONS_CACHE[key] = out
    return [dict(o) for o in out]


def _schedule_options(system, group, H):
    perms = atom_permutations(system, group) if group is not None else [(np.eye(3), list(range(len(system.atoms))))]
    orbs = orbits(system, perms)
    cand = ["lecture"]
    lect = _signature(build_schedule(system, group, H, "lecture"))
    for k, (_, sig) in enumerate(_pairing_ops(system, perms, orbs)):
        cand.append(f"pair:{k}")
    cand += ["central_first", "one_step", "atom_by_atom"]
    out, seen = [], set()
    for v in cand:
        sc = build_schedule(system, group, H, v)
        sig = _signature(sc)
        if sig in seen:
            continue
        seen.add(sig)
        if v.startswith("pair:"):
            pairs = [st["narration"].split(": pair atoms ")[1].rstrip(".") for st in sc["stages"] + sc["skipped"] if st["kind"] == "within_set"]
            title = "Pair differently: " + "; ".join(pairs)
        else:
            title = VARIANT_TITLES[v]
        out.append({"id": v, "title": title, "stages": len(sc["stages"])})
    return out


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
