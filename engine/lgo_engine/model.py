"""System definition and validation (F1.3, F1.8)."""
from dataclasses import dataclass, field
import numpy as np

ORBITAL_TYPES = ("s", "p")  # d and f are later releases


@dataclass(frozen=True)
class Atom:
    element: str
    position: tuple  # Angstrom


@dataclass(frozen=True)
class BasisFunction:
    """One real basis function. For p functions, `direction` is the unit vector d such that
    the function is d . (px, py, pz); inward sigma donors use d = -u (CONVENTIONS section 3)."""
    atom: int
    kind: str
    direction: tuple = None
    label: str = ""


@dataclass
class System:
    name: str
    atoms: list
    basis: list
    edges: list                      # list of (a, b) atom-index pairs with couplings on
    electrons: int
    group: str = "C1"
    linear: bool = False
    central_atom: int = None         # index of the central atom, or None
    preset: str = None
    counting_convention: str = "neutral_atom"
    overrides: dict = field(default_factory=dict)  # {(el_a, orb_a, el_b, orb_b): value_eV}

    def validate(self):
        errs = []
        if not self.atoms:
            errs.append("no atoms")
        if not self.basis:
            errs.append("empty basis")
        for i, f in enumerate(self.basis):
            if f.kind not in ORBITAL_TYPES:
                errs.append(f"basis[{i}]: unknown orbital kind {f.kind!r}")
            if not 0 <= f.atom < len(self.atoms):
                errs.append(f"basis[{i}]: atom index {f.atom} out of range")
            if f.kind == "p":
                d = np.asarray(f.direction, float)
                if d.shape != (3,) or abs(np.linalg.norm(d) - 1) > 1e-12:
                    errs.append(f"basis[{i}]: p direction must be a unit 3-vector")
        for a, b in self.edges:
            if a == b or not (0 <= a < len(self.atoms) and 0 <= b < len(self.atoms)):
                errs.append(f"edge ({a},{b}) invalid")
        if self.electrons < 0 or self.electrons > 2 * len(self.basis):
            errs.append(f"electron count {self.electrons} outside 0..{2 * len(self.basis)}")
        if errs:
            raise ValueError(f"{self.name}: invalid configuration: " + "; ".join(errs))
        return self

    def positions(self):
        return np.array([a.position for a in self.atoms], float)
