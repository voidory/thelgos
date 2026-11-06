"""Comparison of a result with curated reference data (groundwork for the spec section 7 fidelity bar).

Currently the diatomics: orbital ordering (CONVENTIONS A.6) and ground-state spin, the same rules
as scripts/check_release1.py, here computed by the engine so the app can show them.
"""

TEXT = {"pi_u_below_sigma_g_2p": "πu below σg(2p)", "sigma_g_2p_below_pi_u": "σg(2p) below πu"}


def compare(system, result, curated):
    ref = (curated or {}).get("diatomics", {}).get(system.name)
    if not ref:
        return []
    out = []
    levels = result["levels"]
    if ref.get("ordering"):
        sg = [lv["energy_eV"] for lv in levels if lv.get("linear_label") == "sigma_g"]
        pu = [lv["energy_eV"] for lv in levels if lv.get("linear_label") == "pi_u"]
        model = "pi_u_below_sigma_g_2p" if pu[0] < max(sg) else "sigma_g_2p_below_pi_u"
        out.append({"category": "Orbital ordering", "agrees": model == ref["ordering"],
                    "reference": f"Experiment: {TEXT[ref['ordering']]}", "model": f"Model: {TEXT[model]}",
                    "source": ref["source"]})
    open_shell = bool(result["warnings"])
    out.append({"category": "Ground-state spin", "agrees": open_shell == (ref["multiplicity"] > 1),
                "reference": f"Experiment: {ref['ground_state']} ({'open shell' if ref['multiplicity'] > 1 else 'closed shell'})",
                "model": f"Model: {'partly filled degenerate level (open shell)' if open_shell else 'closed shell'}",
                "source": ref["source"]})
    return out
