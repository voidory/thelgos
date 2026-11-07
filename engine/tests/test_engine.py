"""Engine tests: reference spectra, analytic oracles, bond orders, ledger, spec section 12 property
tests and a mutation check. Run with `make test-engine`."""
import json
import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1])]
from load import ROOT, load_params, load_groups, conventions_sha256  # noqa: E402
import lgo_engine as E  # noqa: E402
from lgo_engine import cases, hamiltonian, app_api  # noqa: E402

P = load_params()
G = load_groups()
REF = json.loads((ROOT / "evidence/reference_spectra.json").read_text())


def run(system):
    return E.run(system, P, G, conventions_sha256())


def E2(a, b, g):
    return list(np.linalg.eigvalsh([[a, g], [g, b]]))


class ReferenceSpectra(unittest.TestCase):
    def test_all_release1_fixtures_match_committed_reference(self):
        systems = cases.release1(P)
        self.assertEqual({s.name for s in systems}, set(REF))
        for s in systems:
            r = run(s)
            got = [(lv["irrep"], lv["multiplicity"], lv["electrons"]) for lv in r["levels"]]
            want = [(x["irrep"], x["multiplicity"], x["electrons"]) for x in REF[s.name]]
            self.assertEqual(got, want, s.name)
            for lv, x in zip(r["levels"], REF[s.name]):
                self.assertAlmostEqual(lv["energy_eV"], x["energy_eV"], delta=1e-9, msg=s.name)


class AnalyticOracles(unittest.TestCase):
    """Closed-form blocks from CONVENTIONS Appendix A, written independently of hamiltonian.py."""

    def spectrum(self, s):
        return sorted(l for lv in run(s)["levels"] for l in [lv["energy_eV"]] * lv["multiplicity"])

    def test_lecture_water_a1_b2_b1(self):
        g = P.presets["H2O_lecture"]["geometry_overrides"]
        R, t = g["bond_length_angstrom"], np.deg2rad(g["angle_degrees"]) / 2
        es, ep, eh = P.onsite_energy("O", "s"), P.onsite_energy("O", "p"), P.onsite_energy("H", "s")
        vsp = P.coupling("O", "H", "sp_sigma", R, "H2O_lecture")
        a1 = np.linalg.eigvalsh([[es, 0, 0], [0, ep, -np.sqrt(2) * vsp * np.cos(t)], [0, -np.sqrt(2) * vsp * np.cos(t), eh]])
        want = sorted(list(a1) + E2(ep, eh, -np.sqrt(2) * vsp * np.sin(t)) + [ep])
        np.testing.assert_allclose(self.spectrum(cases.water(P, "H2O_lecture")), want, atol=1e-10)

    def test_sf6_lambda_formulas(self):
        R = P.presets["SF6_model"]["geometry_overrides"]["bond_length_angstrom"]
        c = lambda ch, d: P.coupling("F", "F", ch, d, "SF6_model")
        tc = (c("pp_pi", np.sqrt(2) * R) - c("pp_sigma", np.sqrt(2) * R)) / 2
        tt = -c("pp_sigma", 2 * R)
        ef = P.onsite_energy("F", "p")
        la, lt, le = ef + 4 * tc + tt, ef - tt, ef - 2 * tc + tt
        vsp, vps = P.coupling("S", "F", "sp_sigma", R, "SF6_model"), P.coupling("S", "F", "pp_sigma", R, "SF6_model")
        want = sorted(E2(P.onsite_energy("S", "s"), la, -np.sqrt(6) * vsp) + E2(P.onsite_energy("S", "p"), lt, -np.sqrt(2) * vps) * 3 + [le] * 2)
        np.testing.assert_allclose(self.spectrum(cases.sf6(P)), want, atol=1e-10)

    def test_ligand_only_and_pi(self):
        eh, ec = P.onsite_energy("H", "s"), P.onsite_energy("C", "p")
        R = P.presets["sixH_lecture"]["geometry_overrides"]["radius_angstrom"]
        tc = P.coupling("H", "H", "ss_sigma", np.sqrt(2) * R, "sixH_lecture")
        tt = P.coupling("H", "H", "ss_sigma", 2 * R, "sixH_lecture")
        np.testing.assert_allclose(self.spectrum(cases.six_h(P)), sorted([eh + 4 * tc + tt] + [eh - tt] * 3 + [eh - 2 * tc + tt] * 2), atol=1e-10)
        t = P.coupling("H", "H", "ss_sigma", 1.0, "H3_fixture")
        np.testing.assert_allclose(self.spectrum(cases.ring_h(P, 2, "H3+")), sorted([eh + 2 * t, eh - t, eh - t]), atol=1e-10)
        t = P.coupling("C", "C", "pp_pi", 1.4, "allyl_pi")
        np.testing.assert_allclose(self.spectrum(cases.allyl(P)), sorted([ec + np.sqrt(2) * t, ec, ec - np.sqrt(2) * t]), atol=1e-10)

    def test_homonuclear_diatomic_blocks(self):
        for name in ["B2", "C2", "N2", "O2", "F2", "Al2", "Si2", "P2", "S2", "Cl2"]:
            el, R = name[:-1], P.presets[f"{name}_model"]["geometry_overrides"]["bond_length_angstrom"]
            V = {ch: P.coupling(el, el, ch, R, f"{name}_model") for ch in ["ss_sigma", "sp_sigma", "pp_sigma", "pp_pi"]}
            es, ep = P.onsite_energy(el, "s"), P.onsite_energy(el, "p")
            want = sorted(E2(es + V["ss_sigma"], ep - V["pp_sigma"], -V["sp_sigma"]) + E2(es - V["ss_sigma"], ep + V["pp_sigma"], V["sp_sigma"])
                          + [ep + V["pp_pi"]] * 2 + [ep - V["pp_pi"]] * 2)
            np.testing.assert_allclose(self.spectrum(cases.diatomic(P, name)), want, atol=1e-10, err_msg=name)


class BondOrdersAndFrontier(unittest.TestCase):
    def bo(self, s):
        return run(s)["bond_orders"]

    def test_counting_bond_orders(self):
        for s, want in [(cases.allyl(P), 0.5), (cases.ring_h(P, 2, "H3+"), 1 / 3), (cases.borane(P), 1.0),
                        (cases.methane(P), 1.0), (cases.water(P, "H2O_lecture"), 1.0), (cases.sf6(P), 2 / 3),
                        (cases.diatomic(P, "CO"), 3.0), (cases.diatomic(P, "O2"), 2.0), (cases.diatomic(P, "F2"), 1.0),
                        (cases.diatomic(P, "Si2"), 2.0), (cases.diatomic(P, "P2"), 3.0), (cases.diatomic(P, "S2"), 2.0),
                        (cases.diatomic(P, "Al2"), 1.0), (cases.diatomic(P, "Cl2"), 1.0)]:
            self.assertAlmostEqual(self.bo(s)["counting"]["per_bond"], want, delta=1e-12, msg=s.name)

    def test_coulson_bond_orders(self):
        for v in self.bo(cases.allyl(P))["density_matrix"]["per_bond"].values():
            self.assertAlmostEqual(v, 1 / np.sqrt(2), delta=1e-12)
        for v in self.bo(cases.ring_h(P, 2, "H3+"))["density_matrix"]["per_bond"].values():
            self.assertAlmostEqual(v, 2 / 3, delta=1e-12)

    def test_frontier(self):
        r = run(cases.allyl(P))
        self.assertEqual((r["levels"][r["frontier"]["HOMO"]]["irrep"], r["levels"][r["frontier"]["LUMO"]]["irrep"]), ("B1", "A2"))
        r = run(cases.borane(P))
        self.assertEqual(r["levels"][r["frontier"]["LUMO"]]["irrep"], "A2''")
        self.assertIn("nonbonding", r["frontier"]["LUMO_flags"])
        r = run(cases.diatomic(P, "CO"))
        h = r["levels"][r["frontier"]["HOMO"]]
        self.assertEqual(h["linear_label"], "sigma")
        self.assertGreater(h["composition_atom"]["0"], 0.5)

    def test_partial_shell_warnings(self):
        for s in [cases.ring_h(P, 3, "neutral H3"), cases.six_h(P), cases.four_h(P), cases.diatomic(P, "O2"), cases.diatomic(P, "C2")]:
            self.assertIn("Jahn_Teller_warning", run(s)["warnings"], s.name)
        for s in [cases.ring_h(P, 2, "H3+"), cases.methane(P), cases.diatomic(P, "F2")]:
            self.assertEqual(run(s)["warnings"], [], s.name)


class LedgerAndProperties(unittest.TestCase):
    def test_ledger_totals_equal_eigenvalues(self):
        for s in cases.release1(P):
            H = hamiltonian.build_hamiltonian(s, P)
            e, c, _ = E.solve(H)
            for k in range(len(e)):
                led = E.contact_ledger(H, c[:, k])
                self.assertAlmostEqual(led["total_eV"], e[k], delta=1e-10, msg=s.name)

    def test_zero_coupling_returns_onsites(self):
        for s in cases.release1(P):
            s.edges = []
            H = hamiltonian.build_hamiltonian(s, P)
            np.testing.assert_allclose(np.sort(np.linalg.eigvalsh(H)), np.sort(np.diag(H)))

    def test_mo_count_equals_basis_size(self):
        for s in cases.release1(P):
            r = run(s)
            self.assertEqual(sum(lv["multiplicity"] for lv in r["levels"]), len(s.basis), s.name)

    def test_degenerate_rotation_and_sign_flip_invariance(self):
        rng = np.random.default_rng(0)
        for s in [cases.methane(P), cases.sf6(P), cases.six_h(P), cases.diatomic(P, "O2")]:
            H = hamiltonian.build_hamiltonian(s, P)
            e, c, _ = E.solve(H)
            lv = E.cluster_levels(e, H)
            base, _, _ = E.analyse(s, H, e, c, lv)
            c2 = c.copy()
            for ix in lv:
                Q, _ = np.linalg.qr(rng.normal(size=(len(ix), len(ix))))
                c2[:, ix] = c[:, ix] @ Q
            c2[:, 0] *= -1
            rot, _, _ = E.analyse(s, H, e, c2, lv)
            for a, b in zip(base, rot):
                self.assertEqual(a["character"], b["character"])
                np.testing.assert_allclose(a["composition_ao"], b["composition_ao"], atol=1e-12)
                self.assertAlmostEqual(a["bond_energy_eV"], b["bond_energy_eV"], delta=1e-10)

    def test_result_json_round_trip(self):
        r = run(cases.sf6(P))
        self.assertEqual(json.loads(json.dumps(r, allow_nan=False)), r)

    def test_invalid_config_rejected(self):
        s = cases.methane(P)
        s.basis = []
        with self.assertRaises(ValueError):
            s.validate()
        s = cases.methane(P)
        s.electrons = 99
        with self.assertRaises(ValueError):
            s.validate()


class Schedule(unittest.TestCase):
    """Construction slider (P1, P2, P5; F2.13) and the B5 stage-skipping rule."""

    def sched(self, s):
        H = hamiltonian.build_hamiltonian(s, P)
        return H, E.build_schedule(s, E.Group(G[s.group]), H)

    def test_every_fixture_partitions_edges_and_reaches_full_H(self):
        for s in cases.release1(P):
            H, sc = self.sched(s)
            np.testing.assert_allclose(E.hamiltonian_at(sc, 1.0), H, atol=1e-12, err_msg=s.name)
            np.testing.assert_allclose(E.hamiltonian_at(sc, 0.0), np.diag(np.diag(H)), atol=1e-12, err_msg=s.name)

    def test_stage_counts_and_b5_skips(self):
        want = {"CH4": (1, 2), "H2O model": (1, 1), "BH3": (1, 2), "SF6 model": (4, 0), "H3+": (2, 0),
                "six H Oh": (3, 0), "four H Td": (2, 0), "allyl cation pi": (1, 1), "CO": (1, 0), "O2": (1, 0)}
        for s in cases.release1(P):
            if s.name in want:
                _, sc = self.sched(s)
                self.assertEqual((len(sc["stages"]), len(sc["skipped"])), want[s.name], s.name)
        _, sc = self.sched(cases.methane(P))
        self.assertTrue(all(x["note"] == "No ligand–ligand coupling in this model." for x in sc["skipped"]))

    def test_ch4_ligand_plateau_pairs_levels(self):
        """P1: after fixing and flipping, LGO1' is level with LGO3' and LGO2' with LGO4'."""
        s = cases.four_h(P)
        H, sc = self.sched(s)
        e = np.linalg.eigvalsh(E.hamiltonian_at(sc, 1 / len(sc["stages"])))
        t = P.coupling("H", "H", "ss_sigma", np.sqrt(8 / 3) * 1.087, "fourH_diagnostic")
        eh = P.onsite_energy("H", "s")
        np.testing.assert_allclose(e, sorted([eh + t, eh + t, eh - t, eh - t]), atol=1e-10)

    def test_sf6_first_stage_is_trans_axis_pairs(self):
        s = cases.sf6(P)
        H, sc = self.sched(s)
        self.assertEqual(sc["stages"][0]["edges"], [(1, 2), (3, 4), (5, 6)])
        self.assertEqual(sc["stages"][-1]["kind"], "between_sets")

    def test_checkpoints_land_on_plateaus(self):
        for s in [cases.sf6(P), cases.ring_h(P, 2, "H3+"), cases.methane(P)]:
            _, sc = self.sched(s)
            K = len(sc["stages"])
            for tick in E.checkpoints(sc):
                self.assertTrue(all(l in (0.0, 1.0) for l in E.lambdas(tick["s"], K)), s.name)

    def test_alternative_order_same_endpoint(self):
        """P5: switching the stages on in another order changes the path, not the endpoint."""
        s = cases.sf6(P)
        H, sc = self.sched(s)
        rev = {"H0": sc["H0"], "stages": list(reversed(sc["stages"])), "skipped": []}
        np.testing.assert_allclose(np.linalg.eigvalsh(E.hamiltonian_at(rev, 1.0)), np.linalg.eigvalsh(H), atol=1e-10)


class DiagramSpec(unittest.TestCase):
    """F3.1-F3.4: every level once at its stored energy; the two line layers never mixed."""

    def test_views_for_all_cases(self):
        from lgo_engine import app_api
        for c in app_api.case_list(P):
            v = app_api.view(P, G, c["name"], s=1.0, level=0)
            json.dumps(v, allow_nan=False)
            d = v["diagram"]
            self.assertEqual([x["energy_eV"] for x in d["centre"]], [x["energy_eV"] for x in v["result"]["levels"]])
            self.assertEqual(sum(x["electrons"] for x in d["centre"]), d["electrons"], c["name"])
            for ln in d["lines"]:
                self.assertIn((ln["layer"], ln["style"]), {("mixing", "solid"), ("allowed", "dashed")})
                if ln["layer"] == "mixing":
                    self.assertGreaterEqual(ln["weight"], 0.05)
            frag = d["left"] + d["right"]
            if d["columns"][0] == "left" and d["column_titles"]["right"] == "Ligand group orbitals":
                self.assertTrue(all(x["irrep"] != "—" for x in frag), c["name"])
            self.assertAlmostEqual(v["ledger"]["total_eV"], v["ledger"]["eigenvalue_eV"], delta=1e-10)

    def test_diatomic_textbook_layout(self):
        """Atom 1 | MOs | atom 2; each MO's atomic weights over both sides sum to one."""
        from lgo_engine import app_api
        for name in ["H2", "B2", "C2", "N2", "O2", "F2", "CO", "Al2", "Si2", "P2", "S2", "Cl2"]:
            d = app_api.view(P, G, name)["diagram"]
            self.assertEqual(d["columns"], ["left", "centre", "right"], name)
            self.assertEqual(sum(x["multiplicity"] for x in d["left"] + d["right"]), sum(x["multiplicity"] for x in d["centre"]), name)
        d = app_api.view(P, G, "CO")["diagram"]
        self.assertEqual([x["label"] for x in d["left"]], ["C 2s", "C 2p"])
        self.assertEqual([x["label"] for x in d["right"]], ["O 2s", "O 2p"])
        self.assertEqual((d["column_titles"]["left"], d["column_titles"]["right"]), ("C (atom 1)", "O (atom 2)"))
        d = app_api.view(P, G, "Si2")["diagram"]
        self.assertEqual([x["label"] for x in d["left"]], ["Si 3s", "Si 3p"])

    def test_ch4_lgos_are_symmetry_adapted(self):
        from lgo_engine import app_api
        d = app_api.view(P, G, "CH4")["diagram"]
        self.assertEqual([(x["irrep"], x["multiplicity"]) for x in d["right"]], [("A1", 1), ("T2", 3)])


class ReadingPanel(unittest.TestCase):
    """The app's comparison panel uses the same rules as the checker (scripts/check_release1.py)."""

    def test_diatomic_comparisons_match_checker(self):
        from load import load_extras
        from lgo_engine import app_api
        X = load_extras()
        val = json.loads((ROOT / "checks/release1_validation.json").read_text())
        for r in val["diatomics"]:
            got = {c["category"]: c["agrees"] for c in app_api.view(P, G, r["case"], extras=X)["comparisons"]}
            self.assertEqual(got["Ground-state spin"], r["spin_state"]["agrees"], r["case"])
            if "orbital_ordering" in r:
                self.assertEqual(got["Orbital ordering"], r["orbital_ordering"]["agrees"], r["case"])

    def test_fragment_orbital_selection(self):
        """F3.7: every fragment orbital (LGO or central/atom orbital) is selectable; its ledger sums to its energy."""
        from lgo_engine import app_api
        n = 0
        for s in cases.release1(P):
            d = app_api.view(P, G, s.name)["diagram"]
            for f in d.get("left", []) + d["right"]:
                v = app_api.view(P, G, s.name, fragment=f["id"])
                self.assertAlmostEqual(v["ledger"]["total_eV"], f["energy_eV"], delta=1e-10, msg=(s.name, f["id"]))
                self.assertEqual(len(v["structure"]["contacts"]), sum(1 for r in v["ledger"]["rows"] if r["kind"] == "contact" and abs(r["energy_eV"]) >= 1e-9))
                n += 1
        self.assertGreater(n, 50)

    def test_photoelectron_spectrum(self):
        """Spec Release 1, CH4: two bands with degeneracies 3:1, badged approximate; measured values cited."""
        from load import load_extras
        from lgo_engine import app_api
        X = load_extras()
        p = app_api.view(P, G, "CH4", extras=X)["pes"]
        self.assertEqual(sorted(b["electrons"] for b in p["bands"]), [2, 6])
        self.assertTrue(any("cross-sections" in c for c in p["caveats"]) and any("He I" in c for c in p["caveats"]))
        self.assertAlmostEqual(p["reference"]["first_band_vertical_eV"], 13.60)
        for s in cases.release1(P):
            r = app_api.view(P, G, s.name, extras=X)
            occ = [lv for lv in r["result"]["levels"] if lv["electrons"] > 0]
            self.assertEqual(len(r["pes"]["bands"]), len(occ), s.name)
            self.assertEqual(sum(b["electrons"] for b in r["pes"]["bands"]), s.electrons, s.name)

    def test_every_case_name_in_notes_exists(self):
        from load import load_extras
        names = {s.name for s in cases.release1(P)}
        self.assertTrue(set(load_extras()["notes"]["cases"]) <= names)


class ReleaseTwoDials(unittest.TestCase):
    """Release 2 gate: invariants hold along every dial (D2, D4, D5, D6)."""

    def setUp(self):
        from lgo_engine import app_api
        self.api = app_api

    def test_electron_stepper(self):
        """D6: allyl cation -> radical -> anion, H3+ -> H3 -> H3-; Jahn-Teller warning at neutral H3."""
        api = self.api
        for name, seq in (("allyl cation pi", (2, 3, 4)), ("H3+", (2, 3, 4))):
            charges = []
            for n in seq:
                v = api.view(P, G, name, electrons=n)
                self.assertEqual(sum(x["electrons"] for x in v["result"]["levels"]), n)
                charges.append(v["electrons"]["charge"])
            self.assertEqual(charges, [1, 0, -1], name)
        w = api.view(P, G, "H3+", electrons=3)["result"]["warnings"]
        self.assertIn("Jahn_Teller_warning", w)
        self.assertEqual(api.view(P, G, "H3+", electrons=3)["result"]["levels"][1]["electrons"],
                         api.view(P, G, "neutral H3")["result"]["levels"][1]["electrons"])
        # allyl pi bond order per C-C: cation 0.5, radical 0.5 (the extra electron is nonbonding), anion 0.5
        for n in (2, 3, 4):
            self.assertAlmostEqual(api.view(P, G, "allyl cation pi", electrons=n)["result"]["bond_orders"]["counting"]["per_bond"], 0.5)
        with self.assertRaises(ValueError):
            api.view(P, G, "H2", electrons=5)

    def test_approximation_slider_endpoints(self):
        """D2: at 0 the lecture preset's spectrum, at 1 the model's; SF6 (identical presets) has no slider."""
        api = self.api
        e = lambda v: [x["energy_eV"] for x in v["result"]["levels"]]  # noqa: E731
        lec, mod = api.view(P, G, "H2O lecture"), api.view(P, G, "H2O model")
        np.testing.assert_allclose(e(api.view(P, G, "H2O model", approx=0.0)), e(lec), atol=1e-12)
        np.testing.assert_allclose(e(api.view(P, G, "H2O model", approx=1.0)), e(mod), atol=1e-12)
        mid = e(api.view(P, G, "H2O model", approx=0.5))
        self.assertTrue(e(mod)[0] < mid[0] < e(lec)[0])  # the O 2s level drops steadily as it is allowed to mix
        self.assertFalse(api.view(P, G, "SF6 model")["approx"]["available"])
        with self.assertRaises(ValueError):
            api.view(P, G, "SF6 model", approx=0.5)

    def test_degenerate_dial_invariants(self):
        """D4: rotating within any degenerate set keeps the energy and the set's summed atom weights."""
        api = self.api
        n = 0
        for s in cases.release1(P):
            v0 = api.view(P, G, s.name)
            items = [("level", k) for k, lv in enumerate(v0["result"]["levels"]) if lv["multiplicity"] > 1]
            items += [("fragment", f["id"]) for f in v0["diagram"].get("left", []) + v0["diagram"]["right"] if f["multiplicity"] > 1]
            for kind, key in items:
                ref = None
                for ang in (0, 30, 60, 90):
                    v = api.view(P, G, s.name, dial=ang, **{kind: key})
                    got = (v["ledger"]["total_eV"], v["dial"]["summed_weights"])
                    if ref is None:
                        ref = got
                    self.assertAlmostEqual(got[0], ref[0], delta=1e-10, msg=(s.name, key, ang))
                    np.testing.assert_allclose(got[1], ref[1], atol=1e-12, err_msg=f"{s.name} {key} {ang}")
                n += 1
        self.assertGreater(n, 20)

    def test_lone_pair_dial(self):
        """D5: water's nonbonding pair rotates into two equivalent lone pairs; sum and density unchanged."""
        api = self.api
        for name in ("H2O model", "H2O lecture"):
            for ang in (0, 15, 30, 45):
                lp = api.view(P, G, name, lone_pair_angle=ang)["lone_pairs"]
                self.assertAlmostEqual(lp["sum_eV"], lp["canonical_sum_eV"], delta=1e-10)
                self.assertLess(lp["density_change_max"], 1e-12)
            lp = api.view(P, G, name, lone_pair_angle=45)["lone_pairs"]
            self.assertAlmostEqual(lp["orbitals"][0]["energy_eV"], lp["orbitals"][1]["energy_eV"], delta=1e-10)
        self.assertIsNone(api.view(P, G, "CH4", lone_pair_angle=45)["lone_pairs"])


class NodesRulesPainting(unittest.TestCase):
    """Release 2, Milestone C: L4 rule check, L5 node markers, D3 phase painting (spec 'done when')."""

    def setUp(self):
        from lgo_engine import app_api
        self.api = app_api

    def test_l4_rule_check(self):
        api = self.api
        r = api.view(P, G, "allyl cation pi", level=0)["rule_check"]
        self.assertEqual(r["R6"]["badge"], "approximate")
        self.assertIn("unequal coefficients (1, 1.41, 1)", r["R6"]["reason"])
        self.assertEqual(r["R9"]["badge"], "exact"); self.assertTrue(r["R9"]["agrees"])
        for name, kw in (("CH4", {"fragment": "R0"}), ("CH4", {"fragment": "R1"}), ("H3+", {"level": 0})):
            r6 = api.view(P, G, name, **kw)["rule_check"]["R6"]
            self.assertEqual(r6["badge"], "exact", (name, kw))
            self.assertAlmostEqual(r6["prediction_eV"], r6["exact_eV"], delta=1e-9)
        for name, kw in (("SF6 model", {"fragment": "R2"}), ("six H Oh", {"level": 2})):
            self.assertIn("(2, 2, −1, −1, −1, −1)", api.view(P, G, name, **kw)["rule_check"]["R6"]["reason"], name)

    def test_l5_node_markers(self):
        api = self.api
        v = api.view(P, G, "allyl cation pi", level=1)
        self.assertEqual(v["nodes"]["through_atoms"], [1]); self.assertEqual(v["nodes"]["count"], 1)
        rows = v["structure"]["nodes"]["through_rows"]
        self.assertTrue(rows and all(abs(v["ledger"]["rows"][i]["energy_eV"]) < 1e-9 and v["ledger"]["rows"][i]["node_atom"] == "C2" for i in rows))
        v = api.view(P, G, "allyl cation pi", level=2)
        self.assertEqual(v["nodes"]["count"], 2); self.assertEqual(len(v["structure"]["nodes"]["between"]), 2)
        for b in v["structure"]["nodes"]["between"]:
            self.assertEqual(v["ledger"]["rows"][b["ledger_row"]]["phase"], "out")
        # Every node count equals the sign changes of the engine's eigenvector (chains: along the chain).
        for name in ("allyl cation pi", "H2", "H3+", "neutral H3"):
            v = api.view(P, G, name)
            C = np.array(api.run(api._system(P, name), P, G)["eigenvectors"])
            for k, (x, lv) in enumerate(zip(v["diagram"]["centre"], v["result"]["levels"])):
                self.assertIn("nodes", x, name)
                c = api.canonical_components(C[:, lv["indices"]])[:, 0]
                z = np.abs(c) < 1e-7 * np.max(np.abs(c))
                if name in ("allyl cation pi", "H2"):
                    want = int(z.sum()) + sum(1 for i in range(len(c) - 1) if not z[i] and not z[i + 1] and c[i] * c[i + 1] < 0)
                else:
                    want = (2 * int(z.sum()) + sum(1 for i in range(3) if not z[i] and not z[(i + 1) % 3] and c[i] * c[(i + 1) % 3] < 0)) // 2
                self.assertEqual(x["nodes"], want, (name, k))

    def test_d3_phase_painting(self):
        api = self.api
        for name, signs, label in (("H3+", [1, 1, 1], "A1'"), ("H3+", [1, -1, 0], "E'"), ("allyl cation pi", [1, 0, -1], "A2")):
            p = api.view(P, G, name, paint=signs)["paint"]
            self.assertTrue(p["pure"], (name, signs))
            self.assertEqual([d["irrep"] for d in p["decomposition"] if d["weight"] > 1e-9], [label])
        p = api.view(P, G, "H3+", paint=[1, 0, 0])["paint"]
        self.assertFalse(p["pure"])
        self.assertAlmostEqual(sum(d["weight"] for d in p["decomposition"]), 1.0, delta=1e-12)
        p = api.view(P, G, "SF6 model", paint=[1, 1, 1, 1, 1, 1])["paint"]
        self.assertTrue(p["pure"]); self.assertEqual(p["decomposition"][0]["label"], "LGO1")
        self.assertIsNone(api.view(P, G, "N2")["paint_set"])
        with self.assertRaises(ValueError):
            api.view(P, G, "H3+", paint=[0, 0, 0])


class MorphsAndWalsh(unittest.TestCase):
    """Release 2, Milestone D: morph paths (D1) and correlation plots (F3.8)."""

    def setUp(self):
        from lgo_engine import app_api
        self.api = app_api

    def test_every_path_reproduces_its_case(self):
        api = self.api
        paths = api.morph_paths(P)
        self.assertGreaterEqual(len(paths), 10)
        for p in paths:
            ref = [x["energy_eV"] for x in api.view(P, G, p["case"])["result"]["levels"]]
            v = api.view(P, G, p["case"], morph=p["id"], morph_u=p["case_u"])
            got = [x["energy_eV"] for x in v["result"]["levels"]]
            np.testing.assert_allclose(got, ref, atol=1e-10, err_msg=p["id"])
            self.assertTrue(v["morph"]["at_case"])
            for u in np.linspace(0, 1, 6):
                api.view(P, G, p["case"], morph=p["id"], morph_u=float(u))   # every frame runs

    def test_walsh_tracking(self):
        api = self.api
        w = api.walsh(P, G, "H2O model: bend")
        by = {}
        for L in w["lines"]:
            by.setdefault(L["irrep"], []).append(L)
        for irrep, ls in by.items():   # the same irrep never crosses: rank order holds at every sample
            ls.sort(key=lambda L: L["rank"])
            for a, b in zip(ls, ls[1:]):
                self.assertTrue(all(x <= y + 1e-12 for x, y in zip(a["energies"], b["energies"])), irrep)
        ends = sorted((L["irrep"], L["energies"][-1]) for L in w["lines"])
        b1 = [e for i, e in ends if i == "B1"][0]
        self.assertTrue(any(i == "A1" and abs(e - b1) < 1e-9 for i, e in ends))   # 3a1 and 1b1 meet as the πu pair at 180°
        e = np.array([L["energies"] for L in api.walsh(P, G, "H3+: assemble")["lines"]])
        gaps = np.sort(e, axis=0)[2] - np.sort(e, axis=0)[1]
        self.assertLess(gaps[-1], 1e-9); self.assertTrue(np.all(gaps[:-1] > 1e-6))   # e' degenerate only as a triangle

    def test_path_endpoints(self):
        api = self.api
        t = api.run(api._system(P, "allyl cation pi"), P, G)["hamiltonian_eV"][0][1]
        ec = P.onsite_energy("C", "p")
        got = [x["energy_eV"] for x in api.view(P, G, "allyl cation pi", morph="allyl cation pi: close the ring", morph_u=1)["result"]["levels"]]
        np.testing.assert_allclose(got, [ec + 2 * t, ec - t], atol=1e-10)
        v = api.view(P, G, "CH4", morph="CH4: flatten", morph_u=1)
        lv = v["result"]["levels"]
        pz = [x for x in lv if abs(x["energy_eV"] - P.onsite_energy("C", "p")) < 1e-9]
        hh = [x for x in lv if abs(x["energy_eV"] - P.onsite_energy("H", "s")) < 1e-9]
        self.assertEqual(len(pz), 1); self.assertEqual(len(hh), 1)   # planar: nonbonding C 2pz and a nonbonding H combination
        self.assertEqual(lv[v["result"]["frontier"]["HOMO"]]["energy_eV"], hh[0]["energy_eV"])   # model: H one occupied (no H–H coupling)
        w = api.walsh(P, G, "SF6 model: insert sulfur")
        eg = [L for L in w["lines"] if L["irrep"] == "Eg"][0]["energies"]
        self.assertLess(max(eg) - min(eg), 1e-9)   # the eg LGOs have no sulfur partner and stay put


class PracticeAndExam(unittest.TestCase):
    """Release 2, Milestone B: questions are generated and checked in the engine, against the Miessler levels (M8)."""

    def perfect(self, name):
        from lgo_engine import app_api, exam
        _, system, ctx, H, frags = app_api._quiz_context(P, G, name)
        lv = ctx["result"]["levels"]
        ans = {"ao_order": [a["key"] for a in sorted(exam._ao_items(system, P), key=lambda a: a["energy_eV"])],
               "mo_order": [exam._key(x) for x in lv],
               "rows": [{"character": c["character"], "electrons": x["electrons"]} for x, c in zip(lv, ctx["diagram"]["centre"])],
               "multiplicity": exam._multiplicity(lv)[0],
               "bond_order": ctx["result"]["bond_orders"]["counting"]["per_bond"]}
        fr = exam._fragment_items(system, P, H, ctx["diagram"], frags)
        if fr:
            truth = exam._partner_truth(system, H, *fr)
            ans["partners"] = {k: v["required"] for k, v in truth.items()}
        return ans

    def test_perfect_answers_score_full_marks_on_every_case(self):
        from lgo_engine import app_api
        names = app_api.quiz_cases(P)
        self.assertEqual(len(names), 30)   # 32 fixtures (CONVENTIONS v1.5) less the two ligand-only sets
        for n in names:
            r = app_api.exam_check(P, G, n, self.perfect(n))
            self.assertEqual(r["percent"], 100, (n, [s for s in r["steps"] if not s["correct"]]))

    def test_mistakes_get_reasons(self):
        from lgo_engine import app_api
        a = self.perfect("H2O lecture")
        a["partners"]["L0"] = ["R0"]           # lecture: O 2s decoupled
        a["partners"]["L2"] = ["R0"]           # B1 with A1
        r = app_api.exam_check(P, G, "H2O lecture", a)
        fb = " ".join(f for s in r["steps"] for f in s["feedback"])
        self.assertIn("switched off in this preset", fb)
        self.assertIn("different symmetry", fb)
        a = self.perfect("O2"); a["multiplicity"] = 1
        fb = " ".join(f for s in app_api.exam_check(P, G, "O2", a)["steps"] for f in s["feedback"])
        self.assertIn("2 unpaired electrons", fb)
        self.assertIn("triplet", fb)
        a = self.perfect("CO")
        homo = app_api.view(P, G, "CO")["result"]["frontier"]["HOMO"]
        a["rows"][homo]["character"] = "nonbonding"   # textbook label for the carbon lone pair, accepted
        self.assertEqual(app_api.exam_check(P, G, "CO", a)["percent"], 100)
        a = self.perfect("N2"); a["ao_order"] = a["ao_order"][::-1]
        r = app_api.exam_check(P, G, "N2", a)
        self.assertFalse(r["steps"][0]["correct"]); self.assertIn("Miessler", r["steps"][0]["feedback"][0])
        a = self.perfect("N2"); a["partners"]["Ls"] = ["Rs", "Rpz"]   # s-p mixing accepted
        self.assertEqual(app_api.exam_check(P, G, "N2", a)["percent"], 100)

    def test_bond_order_question_sf6(self):
        """The owner's paper, Q2(d): why the S-F bond order is not one."""
        from lgo_engine import app_api
        r = app_api.exam_check(P, G, "SF6 lecture", {"bonding_electrons": 8, "antibonding_electrons": 0, "bonds": 6, "bond_order": 2 / 3}, kind="bond_order")
        self.assertEqual(r["percent"], 100)
        r = app_api.exam_check(P, G, "SF6 lecture", {"bonding_electrons": 12, "antibonding_electrons": 0, "bonds": 6, "bond_order": 1}, kind="bond_order")
        fb = " ".join(r["steps"][0]["feedback"])
        self.assertIn("nonbonding", fb); self.assertIn("not a whole number", fb); self.assertIn("0.667", fb)
        with self.assertRaises(ValueError):
            app_api.exam_task(P, G, "six H Oh")


class DoubtsAndSchedules(unittest.TestCase):
    """Release 2, Milestone E: alternative schedules (P5) and doubt toggles (P4), each spec done-when."""

    def test_every_schedule_reaches_the_same_endpoint(self):
        """P5: endpoint energies identical to 1e-10 for every schedule; the paths themselves differ."""
        from lgo_engine.schedule import schedule_options, build_schedule, hamiltonian_at
        from lgo_engine.symmetry import Group
        for s in cases.release1(P):
            H = hamiltonian.build_hamiltonian(s, P)
            g = Group(G[s.group]) if s.group != "C1" else None
            e = np.linalg.eigvalsh(H)
            opts = schedule_options(s, g, H)
            self.assertEqual(opts[0]["id"], "lecture")
            sigs = set()
            for o in opts:
                sc = build_schedule(s, g, H, o["id"])
                self.assertLess(np.max(abs(np.linalg.eigvalsh(hamiltonian_at(sc, 1.0)) - e)), 1e-10, (s.name, o["id"]))
                sig = tuple(tuple(map(tuple, st["edges"])) for st in sc["stages"])
                self.assertNotIn(sig, sigs, f"{s.name}: duplicate schedule {o['id']}")
                sigs.add(sig)

    def test_sf6_cis_pairing_takes_another_path(self):
        from lgo_engine.schedule import build_schedule, hamiltonian_at
        from lgo_engine.symmetry import Group
        s = cases.sf6(P)
        H = hamiltonian.build_hamiltonian(s, P)
        g = Group(G["Oh"])
        v = E.app_api.view(P, G, "SF6 model", schedule="pair:1", s=0.25)
        self.assertIn("2-6, 3-7, 4-5", v["schedule"]["stages"][0]["narration"])
        a = np.linalg.eigvalsh(hamiltonian_at(build_schedule(s, g, H), 0.25))
        b = np.linalg.eigvalsh(hamiltonian_at(build_schedule(s, g, H, "pair:1"), 0.25))
        self.assertGreater(np.max(abs(a - b)), 0.1)
        self.assertLess(v["schedule"]["endpoint_max_deviation_eV"], 1e-10)
        with self.assertRaises(ValueError):
            E.app_api.view(P, G, "H2", schedule="central_first")

    def test_expanded_checkpoints_equal_doubted_rules_used(self):
        """P4: expanded checkpoints equal the distrusted rules the case uses, each on a stage boundary."""
        from lgo_engine.views import RULES
        allr = list(RULES)
        for c in E.app_api.case_list(P):
            used = E.app_api.view(P, G, c["name"])["doubt"]["used"]
            self.assertNotIn("R10", used)
            for doubt in ([], ["R2"], ["R3", "R6"], allr):
                v = E.app_api.view(P, G, c["name"], doubt=doubt)
                cps = v["doubt"]["checkpoints"]
                self.assertEqual(sorted(x["rule"] for x in cps), sorted(set(doubt) & set(used)), (c["name"], doubt))
                K = len(v["schedule"]["stages"])
                for x in cps:
                    self.assertTrue(x["lines"] and x["verdict"], (c["name"], x["rule"]))
                    self.assertAlmostEqual(x["s"] * K, round(x["s"] * K), 9)

    def test_doubt_checks_carry_the_engine_numbers(self):
        def cp(name, rule, **kw):
            return next(x for x in E.app_api.view(P, G, name, doubt=[rule], **kw)["doubt"]["checkpoints"] if x["rule"] == rule)
        # R2: water's O 2s-LGO1 pair, |t|/gap > 0.5, so the perturbative estimate is reported as failing.
        r2 = cp("H2O model", "R2")
        row = next(r for r in r2["table"]["rows"] if r[0] == "O 2s–LGO1")
        self.assertGreater(float(row[3]), 0.5)
        self.assertTrue(row[4].startswith("fails"))
        self.assertEqual(r2["alternative"]["kind"], "approx")
        # R3: allyl's missing C1-C3 contact, with the ring-closing morph as the thing to try.
        r3 = cp("allyl cation pi", "R3")
        self.assertIn("C1–C3", r3["lines"][0])
        self.assertEqual(r3["alternative"]["target"], "allyl cation pi: close the ring")
        # R6: exact for H3+ (a ±1 member of e'), approximate for allyl (unequal coefficients); L4's eg reason kept.
        self.assertEqual(cp("H3+", "R6")["verdict"], "exact here")
        self.assertEqual(cp("allyl cation pi", "R6")["verdict"], "approximate here")
        eg = next(r for r in cp("SF6 model", "R6")["table"]["rows"] if r[0] == "LGO3")
        self.assertIn("(2, 2, −1, −1, −1, −1)", eg[3])
        # R8: counting vs density matrix side by side; R4 exact count.
        self.assertIn("0.667", " ".join(cp("SF6 model", "R8")["lines"]))
        self.assertIn("10 valence basis functions in, 10 molecular orbitals out", cp("SF6 model", "R4")["lines"][0])
        # On another schedule R5 and R7 are not part of the route.
        used = E.app_api.view(P, G, "SF6 model", schedule="one_step")["doubt"]["used"]
        self.assertNotIn("R5", used)
        self.assertNotIn("R7", used)


class ConventionsV15(unittest.TestCase):
    """CONVENTIONS v1.5: C3v, NH3, H2S, SiH4, the hydrides NH, OH, HF, and the AH3 morphs."""

    def test_nh3_block_and_labels(self):
        s = cases.ammonia(P)
        H = hamiltonian.build_hamiltonian(s, P)
        g = P.presets["NH3_model"]["geometry_overrides"]
        R, c = g["bond_length_angstrom"], np.cos(np.deg2rad(g["angle_degrees"]))
        b = np.arccos(np.sqrt((2 * c + 1) / 3))
        vss, vsp = 0.6 * -1.40 * 7.62 / R ** 2, 0.6 * 1.84 * 7.62 / R ** 2
        es, ep, eh = P.onsite_energy("N", "s"), P.onsite_energy("N", "p"), P.onsite_energy("H", "s")
        A = [[es, 0, np.sqrt(3) * vss], [0, ep, np.sqrt(3) * vsp * np.cos(b)], [np.sqrt(3) * vss, np.sqrt(3) * vsp * np.cos(b), eh]]
        Ee = E2(ep, eh, -np.sqrt(1.5) * vsp * np.sin(b))
        want = np.sort(np.r_[np.linalg.eigvalsh(A), Ee, Ee])
        self.assertLess(np.max(abs(np.linalg.eigvalsh(H) - want)), 1e-10)
        r = run(s)
        self.assertEqual([lv["irrep"] for lv in r["levels"]], ["A1", "E", "A1", "A1", "E"])
        self.assertEqual(r["levels"][r["frontier"]["HOMO"]]["irrep"], "A1")

    def test_hydride_bond_orders_and_spin(self):
        for name, open_shell in (("NH", True), ("OH", True), ("HF", False)):
            for var, bo in (("model", 2.0), ("lecture", 1.0)):
                v = app_api.view(P, G, f"{name} {var}", extras=__import__("load").load_extras())
                self.assertEqual(v["result"]["bond_orders"]["counting"]["per_bond"], bo, (name, var))
                spin = next(c for c in v["comparisons"] if c["category"] == "Ground-state spin")
                self.assertTrue(spin["agrees"], (name, var))
                self.assertEqual(bool(v["result"]["warnings"]), open_shell)
                if var == "lecture":
                    self.assertEqual(v["result"]["levels"][0]["character"], "nonbonding")
        self.assertEqual(app_api.view(P, G, "NH3")["result"]["bond_orders"]["counting"]["per_bond"], 4 / 3)
        self.assertEqual(app_api.view(P, G, "H2S")["result"]["bond_orders"]["counting"]["per_bond"], 1.5)

    def test_ah3_walsh_bh3_planar_nh3_pyramidal(self):
        from lgo_engine import morph
        def total(pid, k):
            W = morph.walsh(P, G, pid)
            n = W["electrons"]
            es = sorted(((l["energies"][k], l["multiplicity"]) for l in W["lines"]))
            t = 0.0
            for e, m in es:
                put = min(2 * m, n); t += put * e; n -= put
            return t
        self.assertLess(total("BH3: pyramidalise", 0), total("BH3: pyramidalise", -1))   # planar BH3 lower
        self.assertLess(total("NH3: flatten", 0), total("NH3: flatten", -1))             # pyramidal NH3 lower
        W = morph.walsh(P, G, "NH3: flatten")
        lp = next(l for l in W["lines"] if l["irrep"] == "A1" and l["rank"] == 1)
        self.assertAlmostEqual(lp["energies"][-1], P.onsite_energy("N", "p"), places=9)   # pure N 2pz when planar


class InspectorViews(unittest.TestCase):
    """L2, L3, P3, P6, F5.1, F5.2, F5.4-F5.5 on the engine side."""

    def setUp(self):
        from load import load_extras
        from lgo_engine import app_api
        self.api, self.X = app_api, load_extras()

    def test_structure_lines_match_ledger_rows(self):
        for c in self.api.case_list(P):
            v = self.api.view(P, G, c["name"], level=0)
            rows = {k for k, r in enumerate(v["ledger"]["rows"]) if r["kind"] == "contact" and abs(r["energy_eV"]) >= 1e-9}
            self.assertEqual({x["ledger_row"] for x in v["structure"]["contacts"]}, rows, c["name"])

    def test_matching_zero_pattern_follows_irreps(self):
        for c in self.api.case_list(P):
            m = self.api.view(P, G, c["name"])["diagram"]["matching"]
            if not m or m["rows"][0]["irrep"] == "—":
                continue
            for i, r in enumerate(m["rows"]):
                for j, col in enumerate(m["cols"]):
                    if r["irrep"] != col["irrep"]:
                        self.assertEqual(m["tally_eV"][i][j], 0.0, (c["name"], r, col))
        m = self.api.view(P, G, "H2O lecture")["diagram"]["matching"]
        self.assertEqual(m["tally_eV"][0], [0.0, 0.0])  # O 2s decoupled in the lecture preset

    def test_origins(self):
        v = self.api.view(P, G, "SF6 model")
        eg = [k for k, lv in enumerate(v["result"]["levels"]) if lv["irrep"] == "Eg"][0]
        self.assertEqual(v["origins"][eg]["checkpoint"], 3)  # eg is complete once all F-F couplings are on
        self.assertTrue(all(o["checkpoint"] == 4 for k, o in enumerate(v["origins"]) if k != eg))
        v = self.api.view(P, G, "H2O lecture")
        self.assertEqual([o["checkpoint"] for o in v["origins"]], [0, 1, 1, 0, 1, 1])

    def test_explanation_is_deterministic(self):
        a = self.api.view(P, G, "SF6 model")["explanation"]
        self.assertEqual(a, self.api.view(P, G, "SF6 model")["explanation"])
        self.assertIn("Eg have no partner", a)

    def test_stage_badges(self):
        v = self.api.view(P, G, "SF6 model")
        self.assertEqual([st["rules"] for st in v["schedule"]["stages"]], [["R5"], ["R5"], ["R5"], ["R7", "R2"]])
        self.assertEqual(set(v["rules"]), {f"R{i}" for i in range(1, 11)})

    def test_export_and_import(self):
        d = self.api.export(P, G, "CO", s=1.0, level=3, conventions_sha256=conventions_sha256())
        json.loads(json.dumps(d, allow_nan=False))
        ok = self.api.validate_import(P, d, conventions_sha256())
        self.assertEqual(ok["config"], {"case": "CO", "s": 1.0, "level": 3})
        bad = json.loads(json.dumps(d)); bad["config"]["case"] = "XeF2"
        with self.assertRaisesRegex(ValueError, "config.case"):
            self.api.validate_import(P, bad)


class Mutation(unittest.TestCase):
    """Spec section 12: flipping one angular-factor sign must make at least one check fail."""

    def test_sp_sign_flip_is_detected(self):
        original = hamiltonian.element

        def mutated(fa, fb, n, V):
            x = original(fa, fb, n, V)
            return -x if (fa.kind, fb.kind) == ("s", "p") else x
        hamiltonian.element = mutated
        try:
            with self.assertRaises(Exception):
                ReferenceSpectra("test_all_release1_fixtures_match_committed_reference").test_all_release1_fixtures_match_committed_reference()
        finally:
            hamiltonian.element = original


if __name__ == "__main__":
    unittest.main()
