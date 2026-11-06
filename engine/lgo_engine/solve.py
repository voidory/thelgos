"""Diagonalisation and degeneracy grouping (F2.3; tolerances from CONVENTIONS section 10)."""
import numpy as np


def solve(H):
    e, c = np.linalg.eigh(H)
    resid = float(np.linalg.norm(H @ c - c * e, axis=0).max()) if len(e) else 0.0
    ortho = float(np.max(abs(c.T @ c - np.eye(len(e))))) if len(e) else 0.0
    if resid > 1e-10 or ortho > 1e-10:
        raise AssertionError(f"eigensolver check failed: residual {resid}, orthogonality {ortho}")
    # Fix the global sign of each eigenvector: first coefficient above threshold positive (section 4).
    for k in range(c.shape[1]):
        idx = np.flatnonzero(abs(c[:, k]) > 1e-8)
        if len(idx) and c[idx[0], k] < 0:
            c[:, k] *= -1
    return e, c, {"residual_max_eV": resid, "orthogonality_max_error": ortho}


def cluster_levels(e, H):
    """Group eigenvalues into degenerate sets by total span, never by chained near-equalities.
    Tolerance 1e-8 * max(1 eV, ||H||_2), as CONVENTIONS section 10 states."""
    tol = 1e-8 * max(1.0, float(np.linalg.norm(H, 2))) if H.size else 1e-8
    groups, start = [], 0
    for k in range(1, len(e) + 1):
        if k == len(e) or e[k] - e[start] > tol:
            groups.append(list(range(start, k)))
            start = k
    return groups
