"""Parameter scheme (F1.5, F2.12), from data rows shaped like params/*.csv and presets.json.

Lookup precedence (CONVENTIONS section 9.2): preset matrix-element override -> element-pair
row for the channel -> explicit channel default -> error. Absent is never silently zero.
"""


class ParameterSet:
    def __init__(self, onsite_rows, pair_rows, presets, channel_defaults=()):
        self.onsite = {(r["element"], r["orbital"][-1]): float(r["energy_eV"]) for r in onsite_rows}
        self.orbital_name = {(r["element"], r["orbital"][-1]): r["orbital"] for r in onsite_rows}
        self.rows = list(pair_rows)
        self.defaults = list(channel_defaults)
        self.presets = {p["id"]: p for p in presets}
        seen = set()
        for r in self.rows:
            key = tuple(sorted((r["element_a"], r["element_b"]))) + (r["channel"],)
            if key in seen:
                raise ValueError(f"duplicate pair row {key}")
            seen.add(key)

    def onsite_energy(self, element, kind):
        try:
            return self.onsite[(element, kind)]
        except KeyError:
            raise ValueError(f"no on-site energy for {element} {kind}") from None

    def coupling(self, a, b, channel, d, preset):
        """Slater-Koster parameter V for an element pair and channel at distance d (Angstrom)."""
        p = self.presets[preset]
        hits = [r for r in self.rows if sorted((r["element_a"], r["element_b"])) == sorted((a, b))
                and r["channel"] == channel] or [r for r in self.defaults if r["channel"] == channel]
        if len(hits) != 1:
            raise ValueError(f"no unique pair/default row for {a}-{b} {channel}; absent is not zero")
        r = hits[0]
        scale = float(r["scale"])
        if r["s_carrier"] == "H" and abs(scale - p["hydrogen_scale"]) > 1e-15:
            raise ValueError("preset/row hydrogen scale conflict")
        law = r["distance_law"]
        if law == "off":
            return 0.0
        if law == "fixed":
            return float(r["value_eV"]) * scale
        if law == "harrison_inverse_square":
            if d <= 0:
                raise ValueError("nonpositive distance")
            return float(r["eta"]) * float(r["K_eV_A2"]) * scale / d ** 2
        raise ValueError(f"unknown distance law {law!r}")

    def override(self, preset, ea, oa, eb, ob):
        for ov in self.presets[preset]["coupling_overrides"]:
            key = (ov["element_a"], ov["orbital_a"], ov["element_b"], ov["orbital_b"])
            if key == (ea, oa, eb, ob) or key == (eb, ob, ea, oa):
                return float(ov["value_eV"])
        return None
