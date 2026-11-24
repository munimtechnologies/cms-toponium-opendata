# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Unit tests for cmsfit. Stdlib tests run anywhere; numeric tests need numpy and scipy on POSIX
(they run inside the compete job)."""
import copy
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmsfit as cf  # noqa: E402

try:
    import numpy as np
    import scipy  # noqa: F401
except ImportError:  # pragma: no cover
    np = None

HAVE_NP = np is not None and os.name == "posix"


def mini_skim(sample, kind, recs, part=0, nparts=1):
    return {"version": cf.SKIM_VERSION, "sample": sample, "kind": kind, "part": part, "nparts": nparts,
            "max_files": 0, "max_entries_per_file": 0, "constants": cf.skim_constants(), "records": recs}


def mini_rec(xsec, ncls, n_ent, sumw, variations=("nom",), fill=1.0, none=0, uri="f"):
    h = [[fill] * cf.NBIN for _ in range(ncls)]
    return {"title": "t", "xsec_pb": xsec, "files_assigned": [uri],
            "files_done": [{"uri": uri, "entries": n_ent, "sumw_events": sumw, "sumw_runs": sumw}],
            "files_failed": [], "files_skipped_time": [], "h": {v: copy.deepcopy(h) for v in variations},
            "h2_nom": copy.deepcopy(h), "gen_sumw": [fill] * ncls, "n_entries": n_ent, "sumw_events": sumw,
            "sumw_runs": sumw, "gen_source": {"hard": n_ent - none, "fallback": 0, "none": none}}


def good_mc(n_tt=9000000):
    tt = [mini_skim("tt", "tt", {cf.REC_TT: mini_rec(cf.XS_TT2L, cf.NCLS_TT, n_tt // 2, 1.0e6, cf.TT_VARS,
                                                     uri="tt%d" % p)}, part=p, nparts=2) for p in (0, 1)]
    tw = mini_skim("tw", "bkg", {cf.REC_TW: mini_rec(cf.XS_TW, 1, 10, 5.0, cf.BKG_VARS),
                                 cf.REC_TWB: mini_rec(cf.XS_TW, 1, 10, 5.0, cf.BKG_VARS)})
    dy = mini_skim("dy", "bkg", {cf.REC_DY: mini_rec(cf.XS_DY, 1, 10, 5.0, cf.BKG_VARS)})
    return tt + [tw, dy]


class TestStdlib(unittest.TestCase):
    def test_git_blob_sha1(self):
        self.assertEqual(cf.git_blob_sha1(b""), "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391")
        self.assertEqual(cf.git_blob_sha1(b"hello\n"), "ce013625030ba8dba906f756967f9e9ca394464a")

    def test_constants(self):
        self.assertAlmostEqual(cf.XS_TT2L, 88.51, places=2)
        self.assertAlmostEqual(cf.BR_2L2NU, 0.10614564, places=10)
        self.assertEqual(len(cf.scan_grid()), 401)
        self.assertEqual(cf.scan_grid()[88], 8.8)
        self.assertEqual(len(cf.BOUNDS), cf.NPAR)

    def test_class_centres(self):
        c = cf.class_centres()
        self.assertEqual(len(c), 40)
        self.assertAlmostEqual(c[0], -0.975)
        self.assertAlmostEqual(c[-1], 0.975)
        self.assertAlmostEqual(sum(c), 0.0, places=12)

    def test_d_win(self):
        c = cf.class_centres()
        for D in (0.4, -0.7, 0.0):
            G = [(1.0 - D * ck) / 2.0 for ck in c]
            self.assertAlmostEqual(cf.d_win(G), D * (1.0 - 0.000625), places=9)

    def test_spin_weights(self):
        c = cf.class_centres()
        G = [(1.0 - 0.35 * ck) / 2.0 for ck in c]
        self.assertEqual(cf.spin_weights(G, None), [1.0] * 40)
        for tgt in (-1.0, 1.0 / 3.0):
            w = cf.spin_weights(G, tgt)
            self.assertAlmostEqual(cf.d_win([a * b for a, b in zip(w, G)]), tgt, delta=2e-3)
        with self.assertRaises(ValueError):
            cf.spin_weights([0.0] * 39 + [1.0], -1.0)   # D_win = -2.925: reweighting undefined

    def test_gamma_bb(self):
        self.assertEqual(cf.gamma_bb(10.0, 7.0, 0.0), 1.0)
        for nu, n, d in ((10.0, 10.0, 0.1), (1000.0, 1000.0, 0.3), (5.0, 9.0, 0.2), (400.0, 2.0, 0.5),
                         (3.0, 0.0, 0.1), (300.0, 0.0, 0.2)):
            g = cf.gamma_bb(nu, n, d)
            if n == nu:
                self.assertAlmostEqual(g, 1.0, places=12)
            if n > 0:   # stationarity of nu g - n ln(nu g) + (g - 1)^2 / (2 d^2)
                self.assertAlmostEqual(nu - n / g + (g - 1.0) / (d * d), 0.0, delta=1e-6 * nu)
            else:
                self.assertAlmostEqual(g, max(0.0, 1.0 - nu * d * d), places=12)

    def test_grades(self):
        self.assertEqual(cf.grade_h4(3.0), "HIGH")
        self.assertEqual(cf.grade_h4(2.999), "MEDIUM")
        self.assertEqual(cf.grade_h4(1.5), "MEDIUM")
        self.assertEqual(cf.grade_h4(1.49), "LOW")
        self.assertEqual(cf.grade_h1(3.0), "EVIDENCE")
        self.assertEqual(cf.grade_h1(2.0), "HINT")
        self.assertEqual(cf.grade_h1(1.99), "NOT SEEN")
        lim = 2.0 * math.sqrt(1.0 + 1.69)
        self.assertEqual(cf.grade_h2(8.8, 7.8, 9.8), "CONSISTENT")
        self.assertEqual(cf.grade_h2(8.8 + lim - 1e-6, 7.8, 9.8), "CONSISTENT")
        self.assertEqual(cf.grade_h2(8.8 + lim + 1e-6, 7.8, 9.8), "HIGHER")
        self.assertEqual(cf.grade_h2(8.8 - lim - 1e-6, 7.8, 9.8), "LOWER")
        self.assertEqual(cf.grade_h2(8.8, None, 9.8), "NOT GRADABLE")
        self.assertEqual(cf.grade_h3(10.0, 10.5, 10.6), "SINGLET PREFERRED")
        self.assertEqual(cf.grade_h3(10.0, 10.5, 10.4), "UNRESOLVED")
        self.assertEqual(cf.grade_h3(10.0, 9.5, 10.6), "NOT SINGLET")
        self.assertEqual(cf.grade_h3(10.0, 10.0, 10.0), "UNRESOLVED")
        self.assertTrue(cf.gate_g3(0.7) and cf.gate_g3(1.3))
        self.assertFalse(cf.gate_g3(1.31) or cf.gate_g3(0.69))
        self.assertTrue(cf.gate_g4([1.0, 1.005, 1.009], [2.0, 2.0]))
        self.assertFalse(cf.gate_g4([1.0, 1.02], [2.0, 2.0]))

    def test_interval(self):
        mus = [0.0, 0.1, 0.2, 0.3, 0.4]
        self.assertEqual(cf.interval(mus, [3.0, 0.9, 0.0, 1.0, 1.1]), (0.1, 0.3))
        self.assertEqual(cf.interval(mus, [2.0] * 5), (None, None))

    def test_merge_and_g1(self):
        m = cf.merge_skims(good_mc())
        R = m["tt"]["records"][cf.REC_TT]
        self.assertEqual(R["n_entries"], 9000000)
        self.assertEqual(R["sumw_runs"], 2.0e6)
        self.assertEqual(R["h"]["nom"][5][7], 2.0)
        self.assertEqual(R["gen"][3], 2.0)
        self.assertEqual(len(R["files_done"]), 2)
        ok, why = cf.gate_g1(m)
        self.assertTrue(ok, why)

    def test_g1_failures(self):
        ok, why = cf.gate_g1(cf.merge_skims(good_mc(n_tt=7000000)))
        self.assertFalse(ok)
        sk = good_mc()
        sk[0]["records"][cf.REC_TT]["files_done"][0]["sumw_runs"] = 1.0e6 * (1 + 1e-5)
        self.assertFalse(cf.gate_g1(cf.merge_skims(sk))[0])
        sk = good_mc()
        del sk[2]["records"][cf.REC_TWB]
        self.assertFalse(cf.gate_g1(cf.merge_skims(sk))[0])
        sk = good_mc()
        sk[1]["records"][cf.REC_TT]["gen_source"]["none"] = 400000   # has_gen 0.911 in part 1
        ok, why = cf.gate_g1(cf.merge_skims(sk))
        self.assertFalse(ok)
        self.assertTrue(any("has_gen" in w for w in why))

    def test_merge_errors(self):
        sk = good_mc()
        with self.assertRaises(ValueError):
            cf.merge_skims(sk + [copy.deepcopy(sk[0])])       # part given twice
        sk = good_mc()
        sk[0]["version"] = "0.0.9"
        with self.assertRaises(ValueError):
            cf.merge_skims(sk)
        sk = good_mc()
        sk[0]["constants"]["BTAG_WP"] = 0.3
        with self.assertRaises(ValueError):
            cf.merge_skims(sk)
        sk = good_mc()
        sk[0]["constants"]["M_EDGES"] = sk[0]["constants"]["M_EDGES"][:-1] + [1000.0]
        with self.assertRaises(ValueError):
            cf.merge_skims(sk)

    def test_g2(self):
        dG = mini_skim("data_G", "data", {cf.REC_G: mini_rec(None, 1, 5, 0.0)})
        dH = mini_skim("data_H", "data", {cf.REC_H: mini_rec(None, 1, 5, 0.0)})
        self.assertTrue(cf.gate_g2(cf.merge_skims([dG, dH]))[0])
        self.assertFalse(cf.gate_g2(cf.merge_skims([dG]))[0])
        dH2 = copy.deepcopy(dH)
        dH2["records"][cf.REC_H]["files_assigned"] = ["a", "b", "c"]
        self.assertFalse(cf.gate_g2(cf.merge_skims([dG, dH2]))[0])
        dH3 = copy.deepcopy(dH)
        dH3["nparts"] = 2
        self.assertFalse(cf.gate_g2(cf.merge_skims([dG, dH3]))[0])


@unittest.skipUnless(HAVE_NP, "needs numpy and scipy on POSIX")
class TestNumeric(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mc, cls.data = cf.synthetic_skims()
        cls.merged = cf.merge_skims(cls.mc + cls.data)
        cls.tpl = cf.build_templates(cls.merged)

    def test_templates(self):
        t = self.tpl
        tot0 = float((t["low"]["nom"] + t["rest"]["nom"]).sum())
        for name in cf.SHAPE_ONLY:
            for d in ("_up", "_dn"):
                v = name + d
                self.assertAlmostEqual(float((t["low"][v] + t["rest"][v]).sum()) / tot0, 1.0, places=12)
        jes = float((t["low"]["jes_up"] + t["rest"]["jes_up"]).sum()) / tot0
        self.assertGreater(jes, 1.005)                     # JES keeps its rate effect
        self.assertTrue(30000 < tot0 < 120000, tot0)
        for tgt in cf.TARGETS:
            s = t["S"][tgt]["nom"]
            self.assertTrue((s >= 0).all() and s.sum() > 0)
        self.assertTrue((t["delta"] > 0).all() and (t["delta"] < 0.2).all())
        self.assertEqual(t["data"].shape, (cf.NBIN,))
        self.assertAlmostEqual(t["D_win"], 0.3 * (1 - 0.000625), delta=0.02)

    def test_signal_norm(self):
        # S per pb (unchanged template) = L x BR x (selected window / generated window)
        R = self.merged["tt"]["records"][cf.REC_TT]
        H = np.asarray(R["h"]["nom"])[:cf.NWIN].sum()
        G = np.asarray(R["gen"])[:cf.NWIN].sum()
        self.assertAlmostEqual(float(self.tpl["S"]["unchanged"]["nom"].sum()),
                               cf.L_PB * cf.BR_2L2NU * H / G, places=6)

    def test_gradient(self):
        m = cf.Model(self.tpl, "singlet")
        n = m.asimov(8.8) * 1.03
        x = np.array([5.0, 1.04, 0.3, -0.4, 0.5, -0.2, 0.35, -0.6, 0.7, -0.3])
        v, g = m.nll(x, n)
        for i in range(cf.NPAR):
            e = 1e-6
            xp, xm = x.copy(), x.copy()
            xp[i] += e
            xm[i] -= e
            num = (m.nll(xp, n)[0] - m.nll(xm, n)[0]) / (2 * e)
            self.assertAlmostEqual(g[i], num, delta=1e-4 * max(1.0, abs(num)), msg=cf.PARS[i])

    def test_asimov_recovers(self):
        m = cf.Model(self.tpl, "singlet")
        a = cf.analyse(m, m.asimov(8.8), scan=False)
        self.assertAlmostEqual(a["mu_hat"], 8.8, delta=0.2)
        self.assertAlmostEqual(a["k_tt"], 1.0, delta=0.005)
        self.assertLess(a["nll_sb"], 1e-5)
        self.assertGreater(a["Z"], 0.5)
        self.assertTrue(a["G4"], (a["nll_sb_starts"], a["nll_b_starts"]))

    def test_no_signal(self):
        m = cf.Model(self.tpl, "singlet")
        a = cf.analyse(m, m.asimov(0.0), scan=False)
        self.assertLess(a["q0"], 1e-3)
        self.assertLess(a["mu_hat"], 0.5)

    def test_scan_interval(self):
        m = cf.Model(self.tpl, "singlet")
        a = cf.analyse(m, m.asimov(8.8), scan=True)
        self.assertEqual(len(a["scan_q"]), 401)
        self.assertTrue(a["lo68"] <= 8.8 <= a["hi68"], (a["lo68"], a["hi68"]))
        self.assertTrue(min(a["scan_q"]) >= -1e-9)

    def test_pseudo_data(self):
        a = cf.analyse(cf.Model(self.tpl, "singlet"), self.tpl["data"], scan=False)
        self.assertAlmostEqual(a["mu_hat"], 8.8, delta=2.0)     # rounded Asimov: close to injected
        self.assertTrue(cf.gate_g3(a["k_tt"]))

    def test_stat_only_model(self):
        m = cf.Model(self.tpl, "singlet", bb=False)
        fixed = {p: 0.0 for p in cf.PARS[2:]}
        a = cf.analyse(m, m.asimov(8.8), scan=False, fixed=fixed)
        full = cf.analyse(cf.Model(self.tpl, "singlet"), cf.Model(self.tpl, "singlet").asimov(8.8), scan=False)
        self.assertGreaterEqual(a["Z"] + 1e-6, full["Z"])      # nuisances can only reduce Z_A

    def test_main_selftest(self):
        import tempfile
        d = tempfile.mkdtemp()
        rc = cf.main(["--mode", "selftest", "--out-dir", d, "--no-scan"])
        self.assertEqual(rc, 0)
        self.assertTrue(os.path.exists(os.path.join(d, "cmsfit-selftest.json")))

    def test_asimov_mode_refuses_data(self):
        import json
        import tempfile
        d = tempfile.mkdtemp()
        for i, s in enumerate(self.mc + self.data):
            with open(os.path.join(d, "skim-%d.json" % i), "w") as f:
                json.dump(s, f)
        self.assertEqual(cf.main(["--mode", "asimov", "--skim-dir", d, "--out-dir", d, "--no-scan"]), 7)


if __name__ == "__main__":
    unittest.main()
