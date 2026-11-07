"""LGO Method Explorer engine: pure NumPy, no I/O (spec sections 8-9).

The engine takes parameters and point-group tables as plain data and returns a
JSON-serialisable Result. File loading lives outside the package (engine/load.py),
so the same code runs under CPython and Pyodide.

Status: Release 1 (version 1.0.0), built against CONVENTIONS.md v1.4, frozen on the owner's
sign-off (B6). The header of CONVENTIONS.md gives the version the running build uses.
"""
from .model import Atom, BasisFunction, System
from .params import ParameterSet
from .hamiltonian import build_hamiltonian
from .solve import solve, cluster_levels
from .symmetry import Group, label_levels
from .analysis import analyse
from .ledger import contact_ledger
from .result import run
from .schedule import build_schedule, hamiltonian_at, checkpoints, lambdas

__version__ = "1.0.0"
__all__ = ["Atom", "BasisFunction", "System", "ParameterSet", "build_hamiltonian", "solve",
           "cluster_levels", "Group", "label_levels", "analyse", "contact_ledger", "run",
           "build_schedule", "hamiltonian_at", "checkpoints", "lambdas", "__version__"]
