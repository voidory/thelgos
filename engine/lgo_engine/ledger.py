"""Contact ledger (L1, CONVENTIONS section 8 and A.7).

For a normalised real vector c the energy is the diagonal sum c_mu^2 H_mu,mu plus one row per
unordered off-diagonal pair with contribution 2 c_mu c_nu H_mu,nu. The raw contact weight
c_mu c_nu H_mu,nu is half the pair's energy contribution; both are reported and labelled.
"""
import numpy as np


def contact_ledger(H, c_vec, basis=None, tol=0.0):
    c_vec = np.asarray(c_vec, float)
    rows = []
    total = 0.0
    for m in range(len(c_vec)):
        x = c_vec[m] ** 2 * H[m, m]
        total += x
        rows.append({"kind": "onsite", "mu": m, "nu": m, "energy_eV": float(x)})
    for m in range(len(c_vec)):
        for v in range(m + 1, len(c_vec)):
            if H[m, v] == 0 or abs(c_vec[m] * c_vec[v]) <= tol:
                continue
            w = c_vec[m] * c_vec[v] * H[m, v]
            x = 2 * w
            total += x
            rows.append({"kind": "contact", "mu": m, "nu": v, "phase": "in" if c_vec[m] * c_vec[v] > 0 else "out",
                         "contact_weight_eV": float(w), "energy_eV": float(x),
                         "effect": "stabilising" if x < 0 else "destabilising"})
    return {"rows": rows, "total_eV": float(total), "rayleigh_eV": float(c_vec @ H @ c_vec)}
