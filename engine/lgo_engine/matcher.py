"""Case matcher M1-M5 (spec section 3; Release 3 Milestone B). Rule-based and deterministic: no language model.

M1 reads a name, a formula with charge ("C3H5+", "SiF6 2-", "NH4^+") or a SMILES string. M2 finds a structure:
the curated library (the CONVENTIONS fixtures in matcher/library.json and the molecule packs) first; otherwise
VSEPR from the formula, for central-atom molecules, or connectivity from SMILES. VSEPR gives the shape, point
group and template as evidence only: owner decision R3-3 draws nothing without a sourced geometry. M3/M4 give
one of three verdicts with the evidence: drawn, not built yet (with the reason code and what is missing),
outside the method (with why). M5 names the nearest drawn case, what differs, and a morph path if one exists.

Reason codes, checked in this order (the first that applies is the verdict; later ones are listed as also
applying): outside the method (single atom, metal, salt lattice, network solid); not built yet (element without
parameters, d-block element, pi from terminal atoms, multicentre bonding, point group not in library, no
template for this shape, too large for step-by-step, no sourced geometry, structure not determined).
"""
import re

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
           "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
           "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg "
           "Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()
Z = {s: i + 1 for i, s in enumerate(SYMBOLS)}
# Main-group valence electrons (group number for s and p blocks); the d and f blocks are marked separately.
_GROUPS = {1: "H Li Na K Rb Cs Fr", 2: "Be Mg Ca Sr Ba Ra", 3: "B Al Ga In Tl", 4: "C Si Ge Sn Pb", 5: "N P As Sb Bi",
           6: "O S Se Te Po", 7: "F Cl Br I At", 8: "He Ne Ar Kr Xe Rn"}
VALENCE = {e: g for g, s in _GROUPS.items() for e in s.split()}
VALENCE["He"] = 2
D_BLOCK = {s for s in SYMBOLS if s not in VALENCE and Z[s] <= 112}
METALS_12 = set(_GROUPS[1].split()[1:]) | set(_GROUPS[2].split())
HALOGENS = set(_GROUPS[7].split())
NOBLE = set(_GROUPS[8].split())
CHALCOGEN = set(_GROUPS[6].split())
NETWORK = {"O2Si": "silica (SiO2) is a covalent network solid", "BN": "boron nitride is a covalent network solid",
           "CSi": "silicon carbide is a covalent network solid", "C": "carbon as an element is diamond or graphite"}
S_ONLY = {"H", "He"}   # valence shell has no p
STEP_LIMIT = 20   # basis functions; spec: "too large for step-by-step (more than about 20 basis functions)"

# VSEPR shapes AX_nE_m (all X equal): shape name, point group, template.
VSEPR = {(2, 0): ("linear", "Dinfh", "linear AX2"), (2, 1): ("bent", "C2v", "T1"), (2, 2): ("bent", "C2v", "T1"),
         (2, 3): ("linear", "Dinfh", "linear AX2"), (3, 0): ("trigonal planar", "D3h", "T2"),
         (3, 1): ("trigonal pyramidal", "C3v", "AX3E"), (3, 2): ("T-shaped", "C2v", None),
         (4, 0): ("tetrahedral", "Td", "T3"), (4, 1): ("seesaw", "C2v", None), (4, 2): ("square planar", "D4h", None),
         (5, 0): ("trigonal bipyramidal", "D3h", None), (5, 1): ("square pyramidal", "C4v", None),
         (6, 0): ("octahedral", "Oh", "T4")}
# Mixed terminals: the point group the VSEPR shape keeps (None: depends on the isomer).
MIXED = {("AX2E2", (1, 1)): "Cs", ("AX2E1", (1, 1)): "Cs", ("AX4E0", (3, 1)): "C3v", ("AX4E0", (2, 2)): "C2v",
         ("AX3E0", (2, 1)): "C2v", ("AX3E1", (2, 1)): "Cs", ("AX6E0", (5, 1)): "C4v", ("AX6E0", (4, 2)): None}
# Morph paths as case-map edges (D1): path id -> the other case it leads to (None: a frame, not a case).
MORPH_PARTNER = {"BH3: pyramidalise": "NH3", "NH3: flatten": "BH3", "H3+: assemble": "H2"}


# ---------------------------------------------------------------------------------------------- M1: parsing
def hill(counts):
    """Hill-order formula string: C first, then H, then alphabetical; without C, all alphabetical."""
    keys = sorted(counts)
    if "C" in counts:
        keys = ["C"] + (["H"] if "H" in counts else []) + [k for k in keys if k not in ("C", "H")]
    return "".join(f"{k}{counts[k] if counts[k] > 1 else ''}" for k in keys)


def _split_charge(t):
    """Trailing charge. The size of the charge needs an explicit form ("SiF6 2-", "SO4^2-", "SO4(2-)", "SO4-2",
    "Fe+++"); digits run straight into a sign are a count ("NH4+" is NH4 with charge +1, not NH with +4)."""
    for pat in (r"^(.*?)\s*[\^\s]\s*(\d+)\s*([+-])$", r"^(.*?)\(\s*(\d*)\s*([+-])\s*\)$"):
        m = re.match(pat, t)
        if m and m.group(1):
            return m.group(1), (int(m.group(2)) if m.group(2) else 1) * (1 if m.group(3) == "+" else -1)
    m = re.match(r"^(.*?[A-Za-z)\d])([+-])(\d+)$", t)
    if m:
        return m.group(1), int(m.group(3)) * (1 if m.group(2) == "+" else -1)
    m = re.match(r"^(.*?)\s*\^?([+-]+)$", t)
    if m and m.group(1) and len(set(m.group(2))) == 1:
        return m.group(1), len(m.group(2)) * (1 if m.group(2)[0] == "+" else -1)
    return t, 0


def parse_formula(text):
    """'C3H5+' -> ({'C': 3, 'H': 5}, 1). Accepts 'SiF6 2-', 'SO4^2-', 'NH4(+)', parentheses groups ('Al(CH3)3')."""
    t, charge = _split_charge(text.strip())
    t = t.strip()
    if not t or not re.fullmatch(r"[A-Za-z0-9()]+", t):
        return None

    def group(s, i):
        counts = {}
        while i < len(s):
            if s[i] == "(":
                inner, i = group(s, i + 1)
                if i >= len(s) or s[i] != ")":
                    raise ValueError("unbalanced parenthesis")
                i += 1
                n = re.match(r"\d*", s[i:]).group(0)
                i += len(n)
                for k, v in inner.items():
                    counts[k] = counts.get(k, 0) + v * (int(n) if n else 1)
            elif s[i] == ")":
                return counts, i
            else:
                m2 = re.match(r"([A-Z][a-z]?)(\d*)", s[i:])
                if not m2 or m2.group(1) not in Z:
                    raise ValueError(f"unknown element at {s[i:]!r}")
                counts[m2.group(1)] = counts.get(m2.group(1), 0) + (int(m2.group(2)) if m2.group(2) else 1)
                i += len(m2.group(0))
        return counts, i
    try:
        counts, i = group(t, 0)
    except ValueError:
        return None
    if i != len(t) or not counts:
        return None
    return counts, charge


ORGANIC = {"B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5), "S": (2, 4, 6), "F": (1,), "Cl": (1,), "Br": (1,), "I": (1,)}
AROMATIC = {"b": "B", "c": "C", "n": "N", "o": "O", "p": "P", "s": "S"}


def parse_smiles(text):
    """Small SMILES reader (owner decision R3-2): atoms (organic subset, aromatic lowercase, bracket atoms with
    H count and charge), bonds - = # :, branches, ring closures (digits and %nn). Returns dict(atoms, bonds) with
    implicit hydrogens filled in, or None if the string is not SMILES in this subset."""
    s = text.strip()
    if not s or " " in s:
        return None
    atoms, bonds = [], []
    stack, ring, prev, order = [], {}, None, None
    i = 0
    while i < len(s):
        c = s[i]
        if c == "(":
            stack.append(prev); i += 1; continue
        if c == ")":
            if not stack:
                return None
            prev = stack.pop(); i += 1; continue
        if c in "-=#:":
            order = {"-": 1, "=": 2, "#": 3, ":": 1.5}[c]; i += 1; continue
        if c.isdigit() or c == "%":
            if c == "%":
                num, i = s[i + 1:i + 3], i + 3
            else:
                num, i = c, i + 1
            if prev is None:
                return None
            if num in ring:
                a, o2 = ring.pop(num)
                bo = order or o2 or (1.5 if atoms[a]["aromatic"] and atoms[prev]["aromatic"] else 1)
                bonds.append((a, prev, bo))
            else:
                ring[num] = (prev, order)
            order = None
            continue
        if c == "[":
            j = s.find("]", i)
            if j < 0:
                return None
            m = re.fullmatch(r"(\d*)([A-Z][a-z]?|[bcnops])(@*)(H\d*)?([+-]+\d*|[+-]\d+)?", s[i + 1:j])
            if not m:
                return None
            el = m.group(2)
            arom = el in AROMATIC
            el = AROMATIC.get(el, el)
            if el not in Z:
                return None
            h = 0 if not m.group(4) else (int(m.group(4)[1:]) if len(m.group(4)) > 1 else 1)
            ch = 0
            if m.group(5):
                g = m.group(5)
                ch = (int(g[1:]) if g[1:].isdigit() else len(g)) * (1 if g[0] == "+" else -1)
            atoms.append({"element": el, "h": h, "charge": ch, "bracket": True, "aromatic": arom})
            i = j + 1
        else:
            m = re.match(r"Cl|Br|[BCNOPSFI]|[bcnops]", s[i:])
            if not m:
                return None
            el = m.group(0)
            arom = el in AROMATIC
            atoms.append({"element": AROMATIC.get(el, el), "h": None, "charge": 0, "bracket": False, "aromatic": arom})
            i += len(el)
        k = len(atoms) - 1
        if prev is not None:
            bo = order or (1.5 if atoms[prev]["aromatic"] and atoms[k]["aromatic"] else 1)
            bonds.append((prev, k, bo))
        prev, order = k, None
    if stack or ring or not atoms:
        return None
    for k, a in enumerate(atoms):
        if a["h"] is None:
            used = sum(bo for x, y, bo in bonds if k in (x, y))
            if a["aromatic"]:
                used = sum(1 for x, y, bo in bonds if k in (x, y)) + 1
            used = int(round(used))
            ok = [v for v in ORGANIC[a["element"]] if v >= used]
            a["h"] = (ok[0] - used) if ok else 0
    return {"atoms": atoms, "bonds": bonds}


def smiles_formula(g):
    counts = {}
    for a in g["atoms"]:
        counts[a["element"]] = counts.get(a["element"], 0) + 1
        if a["h"]:
            counts["H"] = counts.get("H", 0) + a["h"]
    return counts, sum(a["charge"] for a in g["atoms"])


def charge_text(q):
    return "" if q == 0 else (("+" if q > 0 else "−") if abs(q) == 1 else f"{abs(q)}{'+' if q > 0 else '−'}")


# ---------------------------------------------------------------------------------------------- library
def library(lib, packs):
    """Entries: dict(case, electrons, key, names, smiles, template, shape, source)."""
    out = []
    for e in lib["entries"]:
        key = None
        if e["formula"]:
            f = parse_formula(e["formula"])
            key = f"{hill(f[0])}{e['charge']:+d}"
        out.append({**e, "key": key, "source": "fixture"})
    for p in packs:
        f = parse_formula(p["identity"]["formula"])
        out.append({"case": p["identity"]["name"], "key": f"{hill(f[0])}{p['identity']['charge']:+d}",
                    "names": [p["identity"]["name"]] + p["identity"].get("aliases", []), "smiles": p["identity"].get("smiles"),
                    "template": p["template"], "shape": p["structure"]["shape"], "source": "pack", "formula": p["identity"]["formula"],
                    "charge": p["identity"]["charge"]})
    return out


def _norm(t):
    return re.sub(r"[\s_\-]+", " ", t.strip().lower())


# ---------------------------------------------------------------------------------------------- M2/M3: structure
def central_and_terminals(counts):
    """Central atom for a central-atom formula: the element appearing once that is not H, preferring
    non-halogens, then non-oxygen (H2O: O; HOCl: O; ClF3: Cl; SO3: S). None if no unique choice."""
    once = [e for e, n in counts.items() if n == 1 and e != "H"]
    if not once:
        return None
    pref = ([e for e in once if e not in HALOGENS and e != "O"] or [e for e in once if e not in HALOGENS] or once)
    if len(pref) > 1:
        # Several candidates of the same kind (e.g. HOCl gives O over Cl; ClF3 gives Cl): pick the less
        # electronegative by period position, the heavier halogen first among halogens.
        pref = sorted(pref, key=lambda e: (e in HALOGENS, -Z[e]))
    c = pref[0]
    terms = {e: n for e, n in counts.items() if e != c}
    terms = {e: n for e, n in terms.items() if n}
    return c, terms


def vsepr(c, terms, charge):
    """AXnEm from the formula, single bonds to H and halogens; None with a reason when VSEPR cannot decide."""
    n = sum(terms.values())
    multiple = [e for e in terms if e in CHALCOGEN or e in ("N", "C")]
    if multiple:
        return None, "pi_terminal", f"terminal {', '.join(sorted(multiple))} would form multiple (σ + π) bonds"
    if c not in VALENCE:
        return None, "d_block", f"{c} is a d- or f-block element"
    lp2 = VALENCE[c] - charge - n
    if lp2 < 0 or lp2 % 2:
        return None, "electron_count", f"{VALENCE[c]} valence electrons less charge {charge:+d} cannot make {n} single bonds with paired lone pairs"
    m = lp2 // 2
    if (n, m) not in VSEPR:
        return None, "no_shape", f"AX{n}E{m} has no shape in this matcher"
    return (n, m), None, None


def basis_size(c, terms):
    return 4 + sum(terms.values())


# ---------------------------------------------------------------------------------------------- M4/M5: verdict
def _evidence(**kw):
    return {k: v for k, v in kw.items() if v is not None}


def match(params, groups, text, mode="auto", lib=None, morphs=None):
    """Main entry. Returns dict(input, read_as, formula, verdict, code, reason, also, evidence, case, electrons,
    nearest, choices, alternatives). Verdicts: drawn, not_built, outside, choose (isomers), unread, empty."""
    lib = lib or {"entries": [], "isomers": {}}
    entries = library(lib, getattr(params, "packs", []))
    onsite = {e for e, _ in params.onsite}
    has_p = {e for e, k in params.onsite if k == "p"}
    t = text.strip()
    out = {"input": t, "mode": mode, "choices": [], "also": [], "nearest": None, "alternatives": [], "evidence": {}}
    if not t:
        return {**out, "verdict": "empty", "reason": "Type a name, a formula (with charge, e.g. NH4+) or a SMILES string."}

    # --- M1: name, formula, SMILES (auto tries them in that order; other readings are offered).
    readings = []
    if mode in ("auto", "name"):
        nn = _norm(t)
        hits = [e for e in entries if any(_norm(x) == nn for x in e["names"] + [e["case"]])]
        if hits:
            readings.append(("name", hits[0]))
    f = parse_formula(t) if mode in ("auto", "formula") else None
    g = parse_smiles(t) if mode in ("auto", "smiles") else None
    # A formula that repeats an element symbol ("OO", "CCO") is almost certainly SMILES: read it as SMILES first.
    toks = re.findall(r"[A-Z][a-z]?", re.sub(r"[\^\s(]*\d*[+-]+\)?$", "", t))
    repeats = len(toks) != len(set(toks))
    order = [("smiles", g), ("formula", f)] if (repeats and g) else [("formula", f), ("smiles", g)]
    readings += [(k, d) for k, d in order if d]
    if not readings:
        return {**out, "verdict": "unread", "reason": "Not recognised as a name, a formula or SMILES (subset: organic atoms, bracket atoms, bonds, branches, rings)."}

    def formula_of(kind, d):
        if kind == "name":
            return (d["key"] and (hill(parse_formula(re.match(r"([A-Za-z0-9]+)", d["key"]).group(1))[0])
                                  + charge_text(int(re.search(r"[+-]\d+$", d["key"]).group(0))))) or d["case"]
        if kind == "formula":
            return hill(d[0]) + charge_text(d[1])
        c, q = smiles_formula(d)
        return hill(c) + charge_text(q)
    kind, data = readings[0]
    first = formula_of(kind, data)
    seen = {first}
    for k, d in readings[1:]:
        fx = formula_of(k, d)
        if fx not in seen:
            seen.add(fx)
            out["alternatives"].append({"read_as": k, "formula": fx, "mode": k})
    out["read_as"] = kind
    if kind == "name":
        out["formula"] = first
        return _drawn(out, data, morphs)
    graph = data if kind == "smiles" else None
    counts, charge = (data if kind == "formula" else smiles_formula(graph))
    key = f"{hill(counts)}{charge:+d}"
    ftxt = hill(counts) + charge_text(charge)
    out["formula"] = ftxt
    # --- M2: curated library first; a formula with several structures asks for a choice.
    iso = lib.get("isomers", {}).get(key)
    if iso and graph is None:
        return {**out, "verdict": "choose", "choices": iso, "reason": f"{ftxt} has more than one structure: pick one, or give SMILES.",
                "evidence": {"formula": ftxt}}
    hits = [e for e in entries if e["key"] == key]
    if graph is not None and iso:
        ring = _ring(graph)
        hits = [h for h in hits if (h.get("template") == "T7") == ring and h.get("template") in ("T5", "T6", "T7")] or \
               ([] if ring else hits)
    if hits:
        return _drawn(out, hits[0], morphs)
    ev = {"formula": ftxt}
    elements = set(counts)
    natoms = sum(counts.values())
    # --- outside the method
    if natoms == 1:
        el = next(iter(elements))
        if el in D_BLOCK or el in METALS_12 or el in ("Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi"):
            return {**out, "verdict": "outside", "code": "metal", "reason": f"{el} as an element is a metal: it needs band theory, not a finite molecule.", "evidence": ev}
        if hill(counts) in NETWORK:
            return {**out, "verdict": "outside", "code": "network_solid", "reason": f"{NETWORK[hill(counts)]}: it needs band theory, not a finite molecule.", "evidence": ev}
        return {**out, "verdict": "outside", "code": "single_atom", "reason": "A single atom or ion has atomic orbitals, not molecular orbitals.", "evidence": ev}
    if hill(counts) in NETWORK and charge == 0:
        return {**out, "verdict": "outside", "code": "network_solid", "reason": f"{NETWORK[hill(counts)]}: it needs band theory, not a finite molecule.", "evidence": ev}
    alkali = set(_GROUPS[1].split()[1:])
    earth = {"Mg", "Ca", "Sr", "Ba"}
    if charge == 0 and len(elements) == 2 and elements & (alkali | earth) and elements & (HALOGENS | CHALCOGEN | {"N"}):
        if not (elements & earth and elements & HALOGENS and natoms == 3):
            return {**out, "verdict": "outside", "code": "salt_lattice",
                    "reason": "An ionic salt: as a solid it is a lattice, which needs band theory. (Gas-phase alkaline-earth dihalides such as BaF2 are molecules and are matched as such.)",
                    "evidence": ev}
    # --- M3: structure evidence (template, shape, group) before deciding why it is not built.
    reasons = []
    diatomic = natoms == 2
    c = terms = None
    if diatomic:
        a_, b_ = (list(counts) * 2)[:2] if len(counts) == 1 else sorted(counts)
        ev.update(shape="diatomic", template="T0", group="D∞h" if a_ == b_ else "C∞v", basis_functions=sum(1 if e in S_ONLY else 4 for e in (a_, b_)))
    elif graph is not None and _pi_system(graph):
        n = _heavy(graph)
        tmpl = "T7" if _ring(graph) else "T6"
        ev.update(shape=f"{'ring' if tmpl == 'T7' else 'open chain'} of {n} pπ orbitals", template=tmpl, basis_functions=n)
    else:
        ct = central_and_terminals(counts)
        multi = graph is not None and _heavy(graph) - sum(1 for a in graph["atoms"] if a["element"] in HALOGENS) > 1
        if ct and not multi:
            c, terms = ct
            ev.update(central=c, terminals=", ".join(f"{e}×{k}" for e, k in sorted(terms.items())), equivalent_terminals=len(terms) == 1,
                      basis_functions=basis_size(c, terms))
            shape, code, why = vsepr(c, terms, charge)
            if shape:
                n, m = shape
                name, group, tmpl = VSEPR[shape]
                label = f"AX{n}E{m}"
                if len(terms) > 1:
                    group = MIXED.get((label, tuple(sorted(terms.values(), reverse=True))), "—")
                ev.update(vsepr=f"{label} {name}", template=tmpl or "none", group=group if group else "depends on the isomer")
                if group is None:
                    reasons.append(("isomer", "The point group depends on the isomer (e.g. cis or trans): give SMILES or open a related case."))
                elif {"Dinfh": "D2h"}.get(group, group) not in groups and group != "—":
                    reasons.append(("point_group", f"Point group not in library: {group} (a group data file would add it, spec section 14)."))
                if tmpl is None:
                    reasons.append(("no_template", f"No template for the {name} shape: drawing it would need method steps beyond the taught ones."))
                if group == "—":
                    reasons.append(("point_group", "Point group not determined for this mix of terminal atoms."))
            elif code == "pi_terminal":
                reasons.append(("pi_terminal", f"π from terminal atoms: {why}. A σ-only drawing is not offered automatically yet."))
            elif code:
                reasons.append((code, why))
        else:
            heavy = [e for e, k in counts.items() if e != "H"]
            if len(heavy) == 1 and counts[heavy[0]] == 2 and heavy[0] in ("B", "Al") and natoms - 2 == 6:
                reasons.append(("multicentre", "Multicentre bonding: bridging atoms shared by two centres, which no template covers (e.g. B2H6)."))
            elif graph is not None:
                reasons.append(("no_template", "No template covers this connectivity: it has more than one central atom."))
            else:
                reasons.append(("no_structure", "No central atom can be picked from the formula: give SMILES for the connectivity, or open a related case."))
    # --- parameters (owner decision R3-1) come first among the reasons.
    params_r = []
    dblock = sorted(e for e in elements if e in D_BLOCK)
    missing = sorted(e for e in elements if e not in onsite and e not in D_BLOCK)
    if dblock:
        params_r.append(("d_block", f"{', '.join(dblock)}: d-block element{'s' if len(dblock) > 1 else ''} need d orbitals, which arrive with the advanced complex mode (Release 4)."))
    if missing:
        params_r.append(("element_parameters", f"Element without parameters: {', '.join(missing)} {'is' if len(missing) == 1 else 'are'} beyond the Miessler table this app uses (H to Ca, owner decision R3-1)."))
    nop = sorted({e for e in ([c] if c else [e for e in elements if e not in S_ONLY] if diatomic else []) if e in onsite and e not in has_p and e not in S_ONLY})
    if nop:
        params_r.append(("element_parameters", f"Element without parameters: the table gives {', '.join(nop)} no p energy (owner decision R3-1)."))
    if ev.get("basis_functions", 0) > STEP_LIMIT:
        reasons.append(("too_large", f"Too large for step-by-step ({ev['basis_functions']} basis functions)."))
    reasons = params_r + reasons
    reasons.append(("no_geometry", "Not in the library, so there is no sourced geometry: a molecule pack with a CCCBDB geometry would add it (owner decision R3-3)."))
    code, why = reasons[0]
    near = nearest(counts, charge, ev, entries, morphs)
    return {**out, "verdict": "not_built", "code": code, "reason": why, "also": [{"code": k, "reason": r} for k, r in reasons[1:]],
            "evidence": ev, "nearest": near}


def _heavy(g):
    return sum(1 for a in g["atoms"] if a["element"] != "H")


def _ring(g):
    """True if the heavy-atom graph has a cycle (bonds >= atoms among heavy atoms)."""
    heavy = {k for k, a in enumerate(g["atoms"]) if a["element"] != "H"}
    nb = sum(1 for x, y, _ in g["bonds"] if x in heavy and y in heavy)
    return len(heavy) > 2 and nb >= len(heavy)


def _pi_system(g):
    heavy = [a for a in g["atoms"] if a["element"] != "H"]
    if len(heavy) < 3 or any(a["element"] != "C" for a in heavy):
        return False
    # every carbon has three sigma neighbours (counting H) or carries a charge/radical: sp2
    deg = {k: 0 for k in range(len(g["atoms"]))}
    for x, y, _ in g["bonds"]:
        deg[x] += 1; deg[y] += 1
    return all(deg[k] + (a["h"] or 0) <= 3 for k, a in enumerate(g["atoms"]))


def _drawn(out, hit, morphs):
    # Related cases through morph paths: paths from this case, and paths from another case that lead here.
    rel = [{"morph": m["id"], "title": m["title"], "from": m["case"], "to": MORPH_PARTNER.get(m["id"])} for m in (morphs or [])
           if m["case"] == hit["case"] or MORPH_PARTNER.get(m["id"]) == hit["case"]]
    ev = _evidence(template=hit.get("template"), shape=hit.get("shape"), source="CONVENTIONS fixture" if hit["source"] == "fixture" else "molecule pack",
                   electrons=hit.get("electrons"))
    return {**out, "verdict": "drawn", "case": hit["case"], "electrons": hit.get("electrons"),
            "reason": f"In the library as {hit['case']}" + (f", with {hit['electrons']} electrons" if hit.get("electrons") else "") + ".",
            "evidence": ev, "related": rel}


def nearest(counts, charge, ev, entries, morphs):
    """M5: the closest drawn case. Score: same template (8; 5 for another pi template), same terminal count (3),
    central atom in the same group (4) or the same element (2 more), same terminal elements (2), each shared
    element (2), s-only atoms matching (3 for diatomics), then the closest valence-electron count. Names what
    differs and any morph path from that case."""
    tmpl = ev.get("template")
    c = ev.get("central")
    terms = {}
    if ev.get("terminals"):
        for part in ev["terminals"].split(", "):
            e, k = part.split("×")
            terms[e] = int(k)
    nval = sum(VALENCE.get(x, 0) * n for x, n in counts.items()) - charge
    PI = ("T5", "T6", "T7")
    best, best_score = None, None
    for e in entries:
        if not e.get("key") or e.get("electrons"):
            continue
        m = re.match(r"([A-Za-z0-9]+)([+-]\d+)$", e["key"])
        ec, eq = parse_formula(m.group(1))[0], int(m.group(2))
        score = 0
        if tmpl and e.get("template") == tmpl:
            score += 8
        elif tmpl in PI and e.get("template") in PI:
            score += 5
        ect = central_and_terminals(ec) if sum(ec.values()) > 2 else None
        if c and ect:
            if VALENCE.get(ect[0]) == VALENCE.get(c):
                score += 4
                if ect[0] == c:
                    score += 2
            if sum(ect[1].values()) == sum(terms.values()):
                score += 3
            if set(ect[1]) == set(terms):
                score += 2
        score += 2 * len(set(ec) & set(counts))
        if ev.get("shape") == "diatomic" and sum(ec.values()) == 2:
            if all(x in S_ONLY for x in ec) == all(x in S_ONLY for x in counts):
                score += 3
        ev_e = sum(VALENCE.get(x, 0) * n for x, n in ec.items()) - eq
        key = (score, -abs(ev_e - nval))
        if best_score is None or key > best_score:
            best, best_score = e, key
    if best is None or best_score[0] <= 0:
        return None
    diffs = []
    bc = parse_formula(re.match(r"([A-Za-z0-9]+)", best["key"]).group(1))[0]
    bct = central_and_terminals(bc) if sum(bc.values()) > 2 else None
    if c and bct:
        if bct[0] != c:
            diffs.append(f"{bct[0]} → {c}")
        for x in sorted(set(terms) | set(bct[1])):
            if terms.get(x, 0) != bct[1].get(x, 0):
                diffs.append(f"{x}: {bct[1].get(x, 0)} → {terms.get(x, 0)}")
    else:
        for x in sorted(set(counts) | set(bc)):
            if counts.get(x, 0) != bc.get(x, 0):
                diffs.append(f"{x}: {bc.get(x, 0)} → {counts.get(x, 0)}")
    bq = int(re.search(r"[+-]\d+$", best["key"]).group(0))
    if bq != charge:
        diffs.append(f"charge {bq:+d} → {charge:+d}")
    # A morph is offered only where a path joins the nearest case to this input's own structure (case-map edge).
    path = next((m for m in (morphs or []) if m["case"] == best["case"] and MORPH_PARTNER.get(m["id"]) and
                 any(e["case"] == MORPH_PARTNER[m["id"]] and e.get("key") == f"{hill(counts)}{charge:+d}" for e in entries)), None)
    return {"case": best["case"], "differs": diffs, "morph": path and {"id": path["id"], "title": path["title"]}}
