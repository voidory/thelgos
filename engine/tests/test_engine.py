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
from lgo_engine import cases, hamiltonian  # noqa: E402

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
        for name in ["B2", "C2", "N2", "O2", "F2"]:
            el, R = name[0], P.presets[f"{name}_model"]["geometry_overrides"]["bond_length_angstrom"]
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
                        (cases.diatomic(P, "CO"), 3.0), (cases.diatomic(P, "O2"), 2.0), (cases.diatomic(P, "F2"), 1.0)]:
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
