"""LGO Method Explorer engine: pure NumPy, no I/O (spec sections 8-9).

The engine takes parameters and point-group tables as plain data and returns a
JSON-serialisable Result. File loading lives outside the package (engine/load.py),
so the same code runs under CPython and Pyodide.

Status: provisional, built against CONVENTIONS v1.2 review candidate, which the owner
has not yet signed off.
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

__version__ = "0.1.0-provisional"
__all__ = ["Atom", "BasisFunction", "System", "ParameterSet", "build_hamiltonian", "solve",
           "cluster_levels", "Group", "label_levels", "analyse", "contact_ledger", "run",
           "build_schedule", "hamiltonian_at", "checkpoints", "lambdas", "__version__"]
