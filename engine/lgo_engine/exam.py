"""Practice and Exam modes (Release 2, Milestone B): question generation and answer checking.

Checking runs in the engine, like everything else (spec section 10: the main thread computes no
chemistry). Answers are checked against this app's own model with the Miessler on-sites (DECISIONS.md M8).
Every mistake gets a reason from symmetry, the on-site table or the contact ledger, not just a mark.

Question kinds, modelled on the owner's past paper:
- "build": draw the MO diagram. Steps: order the atomic orbitals; say which orbitals combine; order the
  MOs by symmetry label; mark each MO bonding / nonbonding / antibonding and fill electrons; give the
  spin multiplicity and the bond order.
- "bond_order": count bonding and antibonding electrons and the bonds, give the bond order, and see why.
"""
import numpy as np

QUIZ_EXCLUDE = {"six H Oh", "four H Td"}   # diagnostic ligand-only sets, not molecules
NEAR_NONBONDING_EV = 0.5                    # |contact sum| below this: "nonbonding" also accepted
CHARACTERS = ("bonding", "nonbonding", "antibonding")


def _key(lv):
    return lv.get("linear_label") or lv["irrep"]


GREEK = {"sigma": "σ", "pi": "π", "delta": "δ"}


def _olabel(lv, number=True):
    """Orbital-style label for feedback text: 1a1g, 2σg, 1e' (lower-case Mulliken; Greek for linear)."""
    n = str(lv.get("number", "")) if number else ""
    if lv.get("linear_label"):
        base, _, par = lv["linear_label"].partition("_")
        return f"{n}{GREEK.get(base, base)}{par}"
    ir = lv["irrep"]
    return f"{n}{ir[:1].lower()}{ir[1:]}"


def _contacts(H, V):
    """Contact (off-diagonal) energy of an orbital, averaged over a degenerate set V (columns)."""
    V = np.atleast_2d(np.asarray(V).T).T
    off = H - np.diag(np.diag(H))
    return float(np.trace(V.T @ off @ V) / V.shape[1])


def _multiplicity(levels):
    unpaired = sum(min(lv["electrons"], 2 * lv["multiplicity"] - lv["electrons"]) for lv in levels)
    return unpaired + 1, unpaired


MULT_NAME = {1: "singlet", 2: "doublet", 3: "triplet", 4: "quartet", 5: "quintet"}


def _ao_items(system, params):
    """Distinct atomic-orbital types in the basis (element + orbital), with their on-site energies."""
    seen = {}
    for f in system.basis:
        el = system.atoms[f.atom].element
        name = params.orbital_name.get((el, f.kind), f.kind)
        seen.setdefault(f"{el} {name}", params.onsite_energy(el, f.kind))
    return [{"key": k, "label": k, "energy_eV": e} for k, e in seen.items()]


def _fragment_items(system, params, H, diagram, frag_vectors):
    """Items for the 'which orbitals combine' step, each with a vector in the full basis.

    Central-atom layouts use the diagram's symmetry-adapted columns (central-atom orbitals and LGOs).
    Diatomics use each atom's orbitals split by symmetry about the bond: s (sigma), p along the bond
    (sigma) and the two perpendicular p (pi). Two-column layouts (no central atom) have no such step."""
    if "left" not in diagram["columns"]:
        return None
    lin = system.linear
    if not lin:
        side = lambda col: [{"id": x["id"], "label": x["label"], "sym": x["irrep"], "multiplicity": x["multiplicity"],  # noqa: E731
                             "V": np.array(frag_vectors[x["id"]]), "kind": None} for x in diagram[col]]
        return side("left"), side("right")
    out = ([], [])
    N = len(system.basis)
    for a in (0, 1):
        el = system.atoms[a].element
        groups = {"s": [], "pz": [], "pxy": []}
        for i, f in enumerate(system.basis):
            if f.atom != a:
                continue
            if f.kind == "s":
                groups["s"].append(i)
            elif abs(f.direction[2]) > 0.5:
                groups["pz"].append(i)
            else:
                groups["pxy"].append(i)
        for g, ix in groups.items():
            if not ix:
                continue
            V = np.zeros((N, len(ix)))
            for j, i in enumerate(ix):
                V[i, j] = 1.0
            on = lambda k: params.orbital_name.get((el, k), k)  # noqa: E731
            nm = {"s": f"{el} {on('s')}", "pz": f"{el} {on('p')}z", "pxy": f"{el} {on('p')}x, {on('p')}y"}[g]
            out[a].append({"id": f"{'LR'[a]}{g}", "label": nm, "sym": "σ" if g != "pxy" else "π",
                           "multiplicity": len(ix), "V": V, "kind": g})
    return out


def _partner_truth(system, H, left, right):
    """For each left item: required partners (same symmetry and coupled in this preset), optional ones
    (diatomic s with the other atom's p along the bond: symmetry-allowed s-p mixing), and the rest."""
    truth = {}
    for L in left:
        req, opt, why = [], [], {}
        for R in right:
            same = L["sym"] == R["sym"]
            coup = float(np.linalg.norm(L["V"].T @ H @ R["V"]))
            if same and coup > 1e-9:
                if system.linear and {L["kind"], R["kind"]} == {"s", "pz"}:
                    opt.append(R["id"])
                else:
                    req.append(R["id"])
            why[R["id"]] = (same, coup)
        truth[L["id"]] = {"required": req, "optional": opt, "why": why}
    return truth


def task(system, params, view, H, frag_vectors, kind="build"):
    if system.name in QUIZ_EXCLUDE:
        raise ValueError(f"{system.name!r} is a diagnostic set, not offered as a question")
    lv = view["result"]["levels"]
    centre = view["diagram"]["centre"]
    t = {"kind": kind, "case": system.name, "electrons": system.electrons,
         "question": None, "mo_multiset": [], "ao_items": [], "partners": None}
    if kind == "bond_order":
        t["question"] = (f"Using the MO diagram of {system.name}, count the bonding and antibonding electrons and the "
                         f"equivalent bonds, and give the bond order per bond. If it is not a whole number, say why.")
        return t
    t["question"] = f"Draw the MO diagram of {system.name} ({system.electrons} valence electrons) and give its spin multiplicity and bond order."
    ao = _ao_items(system, params)
    t["ao_items"] = [{"key": a["key"], "label": a["label"]} for a in ao] if len(ao) > 1 else []
    fr = _fragment_items(system, params, H, view["diagram"], frag_vectors)
    if fr:
        t["partners"] = {"left": [{k: x[k] for k in ("id", "label", "sym", "multiplicity")} for x in fr[0]],
                         "right": [{k: x[k] for k in ("id", "label", "sym", "multiplicity")} for x in fr[1]],
                         "left_title": view["diagram"]["column_titles"]["left"], "right_title": view["diagram"]["column_titles"]["right"]}
    counts = {}
    for x in centre:
        k = _key(x)
        if k not in counts:
            counts[k] = {"key": k, "irrep": x["irrep"], "linear_label": x.get("linear_label"), "multiplicity": x["multiplicity"], "count": 0}
        counts[k]["count"] += 1
    t["mo_multiset"] = list(counts.values())
    return t


def check(system, params, view, H, frag_vectors, answers, kind="build"):
    lv = view["result"]["levels"]
    C = np.array(view["_eigenvectors"])
    steps = []
    bo = view["result"]["bond_orders"]["counting"]
    fmt = lambda e: f"{e:+.2f} eV".replace("-", "−")  # noqa: E731

    def step(name, ok_items, total, feedback):
        steps.append({"step": name, "score": ok_items, "max": total, "correct": ok_items == total, "feedback": feedback})

    if kind == "bond_order":
        fb, ok = [], 0
        want = {"bonding_electrons": bo["bonding_electrons"], "antibonding_electrons": bo["antibonding_electrons"], "bonds": bo["bonds"]}
        for k, w in want.items():
            got = answers.get(k)
            if got is not None and int(got) == w:
                ok += 1
            else:
                fb.append(f"{k.replace('_', ' ').capitalize()}: {w}, not {got}.")
        bo_ok = answers.get("bond_order") is not None and abs(float(answers["bond_order"]) - bo["per_bond"]) < 0.01
        ok += bo_ok
        if not bo_ok:
            fb.append(f"Bond order: ½ × ({bo['bonding_electrons']} − {bo['antibonding_electrons']}) / {bo['bonds']} = {bo['per_bond']:.3f}.")
        nb = [x for x in view["diagram"]["centre"] if x["character"] == "nonbonding" and x["electrons"] > 0]
        why = []
        for x in view["diagram"]["centre"]:
            if x["electrons"]:
                why.append(f"{_olabel(x)} ({x['electrons']} e⁻) {x['character']}")
        expl = "Occupied MOs: " + "; ".join(why) + "."
        if nb:
            expl += (" " + ", ".join(_olabel(x) for x in nb) +
                     f" hold{'s' if len(nb) == 1 else ''} {sum(x['electrons'] for x in nb)} electrons that count neither way: "
                     "they have no partner of their symmetry on the other side, so they stay nonbonding. That is why the bond order "
                     + ("is not a whole number." if abs(bo['per_bond'] - round(bo['per_bond'])) > 1e-9 else "is lower than the number of electron pairs suggests."))
        step("Bond order", ok, 4, fb + [expl])
        return _summary(steps)

    # 1. Atomic-orbital order (on-site energies, Miessler table).
    ao = {a["key"]: a for a in _ao_items(system, params)}
    if len(ao) > 1:
        order = answers.get("ao_order") or []
        fb = []
        good = sorted(order) == sorted(ao) and all(ao[a]["energy_eV"] <= ao[b]["energy_eV"] + 1e-9 for a, b in zip(order, order[1:]))
        if not good:
            for a, b in zip(order, order[1:]):
                if a in ao and b in ao and ao[a]["energy_eV"] > ao[b]["energy_eV"] + 1e-9:
                    fb.append(f"{b} ({fmt(ao[b]['energy_eV'])}) lies below {a} ({fmt(ao[a]['energy_eV'])}) in the Miessler table.")
                    break
            fb.append("Lowest first: " + ", ".join(f"{k} {fmt(v['energy_eV'])}" for k, v in sorted(ao.items(), key=lambda kv: kv[1]["energy_eV"])) + ".")
        step("Order the atomic orbitals", int(good), 1, fb)

    # 2. Which orbitals combine.
    fr = _fragment_items(system, params, H, view["diagram"], frag_vectors)
    if fr:
        left, right = fr
        names = {x["id"]: x["label"] for x in left + right}
        syms = {x["id"]: x["sym"] for x in left + right}
        truth = _partner_truth(system, H, left, right)
        got_all = answers.get("partners") or {}
        ok, fb = 0, []
        for L in left:
            T = truth[L["id"]]
            got = set(got_all.get(L["id"], []))
            req, opt = set(T["required"]), set(T["optional"])
            if req <= got <= req | opt:
                ok += 1
                if opt and not (opt & got):
                    fb.append(f"{L['label']}: right. It may also mix with {', '.join(names[r] for r in opt)} (same σ symmetry: s–p mixing).")
                continue
            for r in sorted(req - got):
                fb.append(f"{L['label']} ({L['sym']}) combines with {names[r]} ({syms[r]}): same symmetry, and they are coupled "
                          f"(coupling {T['why'][r][1]:.2f} eV).")
            for r in sorted(got - req - opt):
                same, coup = T["why"][r]
                if not same:
                    fb.append(f"{L['label']} is {L['sym']} and {names[r]} is {syms[r]}: different symmetry, so they cannot combine.")
                else:
                    fb.append(f"{L['label']} and {names[r]} share {L['sym']} symmetry, but their coupling is switched off in this preset "
                              "(the lecture's energy-gap simplification), so they do not combine.")
        step("Which orbitals combine", ok, len(left), fb)

    # 3. MO order, by symmetry label, lowest first.
    seq = answers.get("mo_order") or []
    truth_seq = [_key(x) for x in lv]
    ok_pos = []
    for k, key in enumerate(seq):
        if k < len(lv):
            acceptable = {_key(x) for x in lv if abs(x["energy_eV"] - lv[k]["energy_eV"]) < 1e-6}
            ok_pos.append(key in acceptable)
    good = len(seq) == len(truth_seq) and all(ok_pos) and sorted(seq) == sorted(truth_seq)
    fb = []
    if not good:
        bad = next((k for k, o in enumerate(ok_pos) if not o), None)
        if bad is not None:
            want = view["diagram"]["centre"][bad]
            got_lv = next((x for x in view["diagram"]["centre"] if _key(x) == seq[bad]), {"irrep": seq[bad]})
            fb.append(f"Position {bad + 1} from the bottom is {_olabel(want, False)} ({fmt(lv[bad]['energy_eV'])}), not {_olabel(got_lv, False)}.")
        fb.append("Lowest first: " + ", ".join(_olabel(x) for x in view["diagram"]["centre"]) + ".")
    step("Order the MOs", int(good), 1, fb)

    # 4. Character and electrons, per MO in the student's own order (positions that match the true order).
    rows = answers.get("rows") or []
    ok_c = ok_e = n = 0
    fb = []
    for k, key in enumerate(seq):
        if k >= len(lv) or not ok_pos[k] or k >= len(rows):
            continue
        n += 1
        x = lv[k]
        cx = view["diagram"]["centre"][k]
        Vk = C[:, x["indices"]]
        cont = _contacts(H, Vk)
        allowed = {cx["character"]} | ({"nonbonding"} if abs(cont) < NEAR_NONBONDING_EV else set())
        got = rows[k] or {}
        if got.get("character") in allowed:
            ok_c += 1
        else:
            reason = "stabilising" if cont < 0 else "destabilising"
            fb.append(f"{_olabel(cx)}: {cx['character']}. Its contacts add up to {fmt(cont)} ({reason}"
                      + (", small enough that 'nonbonding' is also accepted" if abs(cont) < NEAR_NONBONDING_EV else "") + ").")
        if got.get("electrons") is not None and int(got["electrons"]) == x["electrons"]:
            ok_e += 1
        else:
            fb.append(f"{_olabel(cx)}: {x['electrons']} electrons (filling from the bottom, two per orbital; this level has "
                      f"{x['multiplicity']} orbital{'s' if x['multiplicity'] > 1 else ''}, room for {2 * x['multiplicity']}).")
    skipped = len(lv) - n
    if skipped:
        fb.append(f"{skipped} MO{'s' if skipped > 1 else ''} not marked here because {'they are' if skipped > 1 else 'it is'} in the wrong position in the previous step.")
    step("Bonding character", ok_c, len(lv), [f for f in fb if "contacts" in f] + [f for f in fb if "not marked" in f])
    step("Fill the electrons", ok_e, len(lv), [f for f in fb if "electrons (" in f])

    # 5. Spin multiplicity (Hund's rule) and bond order (counting).
    mult, unpaired = _multiplicity(lv)
    got = answers.get("multiplicity")
    ok = got is not None and int(got) == mult
    step("Spin multiplicity", int(ok), 1, [] if ok else [
        f"{unpaired} unpaired electron{'s' if unpaired != 1 else ''} (Hund's rule: a partly filled degenerate level keeps its electrons "
        f"parallel), so 2S + 1 = {mult}, a {MULT_NAME.get(mult, str(mult))}."])
    got = answers.get("bond_order")
    ok = got is not None and abs(float(got) - bo["per_bond"]) < 0.01
    step("Bond order", int(ok), 1, [] if ok else [
        f"½ × ({bo['bonding_electrons']} bonding − {bo['antibonding_electrons']} antibonding electrons) over {bo['bonds']} "
        f"equivalent bond{'s' if bo['bonds'] > 1 else ''} = {bo['per_bond']:.3f}."])
    return _summary(steps)


def _summary(steps):
    score = sum(s["score"] / s["max"] for s in steps if s["max"])
    return {"steps": steps, "score": score, "max": len(steps), "percent": round(100 * score / max(len(steps), 1))}
