# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Unit tests for cmsskim. String/stdlib tests run anywhere; numeric tests need numpy on POSIX
(they run inside the compete job); selection tests also need awkward."""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmsskim as cs  # noqa: E402

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None
try:
    import awkward as ak
except ImportError:  # pragma: no cover
    ak = None

HAVE_NP = np is not None and os.name == "posix"
HAVE_AK = HAVE_NP and ak is not None

U = ("root://eospublic.cern.ch//eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/"
     "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/120000/A.root")

PS_V9 = ("PS weights (w_var / w_nominal);   [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2"
         "[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;")
PS_ALT = "PS weights; [0] is ISR=0.5 FSR=1; [1] is ISR=1 FSR=0.5; [2] is ISR=2 FSR=1; [3] is ISR=1 FSR=2"
SC_RENS = ("LHE scale variation weights (w_var / w_nominal); "
           "[0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; "
           "[2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; "
           "[4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; "
           "[6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; "
           "[8] is renscfact=2d0 facscfact=2d0 ")
SC_MUF = ('LHE scale variation weights (w_var / w_nominal); [0] is MUF="0.5" MUR="0.5"; '
          '[1] is MUF="1.0" MUR="0.5"; [2] is MUF="2.0" MUR="0.5"; [3] is MUF="0.5" MUR="1.0"; '
          '[4] is MUF="2.0" MUR="1.0"; [5] is MUF="0.5" MUR="2.0"; [6] is MUF="1.0" MUR="2.0"; '
          '[7] is MUF="2.0" MUR="2.0"')


class TestStdlib(unittest.TestCase):
    def test_extract(self):
        obj = {"m": {"_file_indices": [{"files": [{"uri": U}, {"uri": U.replace("A.root", "B.root")}]},
                                       {"files": [{"uri": U}, {"uri": "root://eospublic.cern.ch//eos/x.json"}]}]}}
        self.assertEqual(cs.extract_root_uris(obj), [U, U.replace("A.root", "B.root")])

    def test_url_forms(self):
        f = cs.url_forms(U)
        self.assertEqual(f[0], "https://opendata.cern.ch/eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/"
                               "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/"
                               "106X_mcRun2_asymptotic_v17-v1/120000/A.root")
        self.assertEqual(f[2], U)
        with self.assertRaises(ValueError):
            cs.url_forms("root://other.cern.ch//eos/x.root")

    def test_split_files(self):
        u = ["c", "a", "b", "d", "e"]
        self.assertEqual(cs.split_files(u, 0, 2), ["a", "c", "e"])
        self.assertEqual(cs.split_files(u, 1, 2), ["b", "d"])
        self.assertEqual(cs.split_files(u, 0, 2, max_files=2), ["a", "c"])
        self.assertEqual(cs.split_files(u, 0, 1), ["a", "b", "c", "d", "e"])
        with self.assertRaises(ValueError):
            cs.split_files(u, 2, 2)

    def test_parse_ps(self):
        self.assertEqual(cs.parse_ps_doc(PS_V9), {"isr_up": 0, "fsr_up": 1, "isr_dn": 2, "fsr_dn": 3})
        self.assertEqual(cs.parse_ps_doc(PS_ALT), {"isr_dn": 0, "fsr_dn": 1, "isr_up": 2, "fsr_up": 3})
        with self.assertRaises(ValueError):
            cs.parse_ps_doc("PS weights (w_var / w_nominal); [0] is ISR=2 FSR=1")
        with self.assertRaises(ValueError):
            cs.parse_ps_doc(None)

    def test_parse_scale(self):
        self.assertEqual(cs.parse_scale_doc(SC_RENS), {"mur_up": 7, "mur_dn": 1, "muf_up": 5, "muf_dn": 3})
        self.assertEqual(cs.parse_scale_doc(SC_MUF), {"mur_up": 6, "mur_dn": 1, "muf_up": 4, "muf_dn": 3})
        with self.assertRaises(ValueError):
            cs.parse_scale_doc("nothing here")

    def test_constants(self):
        self.assertEqual(cs.NRECO, 36)
        self.assertEqual(len(cs.WIN_EDGES), 41)
        self.assertAlmostEqual(cs.WIN_EDGES[-1], 1.0)
        self.assertTrue(all(a < b for a, b in zip(cs.M_EDGES, cs.M_EDGES[1:])))
        self.assertEqual(cs.NCLS_TT, 43)
        self.assertAlmostEqual(cs.BR_2L2NU, 0.10614564, places=7)
        self.assertEqual(len(set(cs.TT_VARS)), 11)


def _to_ptetaphim(p):
    px, py, pz = p[1], p[2], p[3]
    pt = math.hypot(px, py)
    return pt, math.asinh(pz / pt), math.atan2(py, px), math.sqrt(max(p[0] ** 2 - px * px - py * py - pz * pz, 0.0))


def _spin_system():
    """Top pair with known lepton directions in the parent rest frames, boosted to a lab frame."""
    mt = 172.5
    q = np.array([30.0, 40.0, 0.0])
    et = math.sqrt(q @ q + mt * mt)
    t = np.array([et, *q])
    tb = np.array([et, *(-q)])
    u1 = np.array([1.0, 2.0, 2.0]) / 3.0
    u2 = np.array([0.0, -0.6, 0.8])
    lp = cs.boost(np.array([30.0, *(30.0 * u1)]), q / et)
    lm = cs.boost(np.array([25.0, *(25.0 * u2)]), -q / et)
    B = np.array([0.1, -0.2, 0.6])
    return [cs.boost(x, B) for x in (t, tb, lp, lm)], float(u1 @ u2), 2.0 * et


@unittest.skipUnless(HAVE_NP, "numpy on POSIX only")
class TestNumeric(unittest.TestCase):
    def test_p4_mass(self):
        p = cs.p4(3.0, 0.0, 0.0, 4.0)
        self.assertAlmostEqual(p[0], 5.0)
        self.assertAlmostEqual(float(cs.mass(p)), 4.0)

    def test_rest_and_inverse(self):
        P = cs.p4(50.0, 1.2, -0.7, 300.0)
        r = cs.to_rest(P, P)
        self.assertAlmostEqual(r[0], 300.0, places=8)
        self.assertLess(float(np.abs(r[1:]).max()), 1e-8)
        p = np.array([10.0, 1.0, 2.0, 3.0])
        b = np.array([0.3, -0.1, 0.5])
        self.assertLess(float(np.abs(cs.boost(cs.boost(p, b), -b) - p).max()), 1e-10)

    def test_cos_phi_gen(self):
        (t, tb, lp, lm), expect, _ = _spin_system()
        self.assertAlmostEqual(expect, 2.0 / 15.0, places=12)
        self.assertAlmostEqual(float(cs.cos_phi_gen(t, tb, lp, lm)), expect, places=9)
        # vectorised over two events
        T = np.stack([t, t]); TB = np.stack([tb, tb]); LP = np.stack([lp, lp]); LM = np.stack([lm, lm])
        self.assertTrue(np.allclose(cs.cos_phi_gen(T, TB, LP, LM), [expect, expect], atol=1e-9))

    def test_cos_star(self):
        u1 = np.array([1.0, 0.0, 0.0]); u2 = np.array([0.6, 0.8, 0.0])
        X = np.array([300.0, 0.0, 0.0, 0.0])
        l1 = np.array([50.0, *(50.0 * u1)]); l2 = np.array([40.0, *(40.0 * u2)])
        B = np.array([0.1, -0.2, 0.6])
        got = cs.cos_star(cs.boost(l1, B), cs.boost(l2, B), cs.boost(X, B))
        self.assertAlmostEqual(float(got), 0.6, places=9)

    def test_dphi(self):
        self.assertAlmostEqual(float(cs.dphi(3.0, -3.0)), 2 * math.pi - 6.0, places=12)
        self.assertAlmostEqual(float(cs.dphi(0.0, 2.0)), 2.0, places=12)

    def test_bins(self):
        self.assertEqual(list(cs.bin_index([0.0, 159.9, 160.0, 1000.0, 650.0], cs.M_EDGES)), [0, 0, 1, 8, 8])
        self.assertEqual(list(cs.bin_index([-1.0, -0.5, 0.999, 1.0], cs.COS_EDGES)), [0, 1, 3, 3])
        self.assertEqual(int(cs.reco_index(np.array([250.0]), np.array([0.2]))[0]), 14)

    def test_gen_class(self):
        got = cs.gen_class(np.array([340.0, 349.9, 360.0, 500.0, np.nan]),
                           np.array([-1.0, 0.99, 0.0, 0.0, np.nan]),
                           np.array([True, True, True, True, False]))
        self.assertEqual(list(got), [0, 39, 40, 41, 42])

    def test_lumi_mask(self):
        golden = {"1": [[1, 5], [10, 12]], "2": [[3, 3]]}
        got = cs.lumi_mask([1, 1, 1, 2, 2, 3], [1, 6, 11, 3, 4, 1], golden)
        self.assertEqual(list(got), [True, False, True, True, False, False])


def _events():
    e = dict(pt=[40.0], eta=[0.5], phi=[0.0], q=[-1])
    base_j = dict(pt=[60.0, 45.0, 35.0], eta=[0.1, 1.0, 0.5], phi=[-2.0, 1.0, 0.1], m=[10.0, 8.0, 5.0],
                  b=[0.9, 0.1, 0.1])
    rows = []
    for k in range(6):
        j = {kk: list(v) for kk, v in base_j.items()}
        mq, trig = 1, True
        if k == 1:
            mq = -1                       # same sign
        if k == 2:
            j["b"] = [0.1, 0.1, 0.9]      # only the overlapping jet is b-tagged
        if k == 3:
            j["pt"] = [60.0, 30.5, 35.0]  # second jet fails after -3 % JES
        if k == 4:
            j = dict(pt=[80.0, 50.0, 40.0], eta=[-1.0, 1.5, 2.0], phi=[-1.0, -2.5, 0.8], m=[9.0, 7.0, 6.0],
                     b=[0.1, 0.9, 0.1])
        if k == 5:
            trig = False
        rows.append(dict(e=e, mq=mq, j=j, trig=trig))
    d = {
        "PV_npvsGood": [5] * 6,
        "Electron_pt": [r["e"]["pt"] for r in rows], "Electron_eta": [r["e"]["eta"] for r in rows],
        "Electron_phi": [r["e"]["phi"] for r in rows], "Electron_charge": [r["e"]["q"] for r in rows],
        "Electron_cutBased": [[4]] * 6, "Electron_deltaEtaSC": [[0.0]] * 6,
        "Muon_pt": [[30.0]] * 6, "Muon_eta": [[-0.3]] * 6, "Muon_phi": [[2.0]] * 6,
        "Muon_charge": [[r["mq"]] for r in rows], "Muon_tightId": [[True]] * 6,
        "Muon_pfRelIso04_all": [[0.05]] * 6,
        "Jet_pt": [r["j"]["pt"] for r in rows], "Jet_eta": [r["j"]["eta"] for r in rows],
        "Jet_phi": [r["j"]["phi"] for r in rows], "Jet_mass": [r["j"]["m"] for r in rows],
        "Jet_jetId": [[6, 6, 6]] * 6, "Jet_btagDeepFlavB": [r["j"]["b"] for r in rows],
        "HLT_A": [r["trig"] for r in rows], "Flag_B": [True] * 6,
    }
    return ak.Array(d)


@unittest.skipUnless(HAVE_AK, "awkward + numpy on POSIX only")
class TestSelection(unittest.TestCase):
    def test_masks(self):
        ev = _events()
        nom = cs.select_reco(ev, 1.0, ["HLT_A"], ["Flag_B"])
        self.assertEqual(list(nom["mask"]), [True, False, False, True, True, False])
        dn = cs.select_reco(ev, 0.97, ["HLT_A"], ["Flag_B"])
        self.assertEqual(list(dn["mask"]), [True, False, False, False, True, False])
        up = cs.select_reco(ev, 1.03, ["HLT_A"], ["Flag_B"])
        self.assertEqual(list(up["mask"]), [True, False, False, True, True, False])
        self.assertTrue(np.isnan(nom["m_llbb"][1]))

    def test_observables(self):
        ev = _events()
        r = cs.select_reco(ev, 1.0, ["HLT_A"], ["Flag_B"])
        E = cs.p4(40.0, 0.5, 0.0, cs.ME)
        M = cs.p4(30.0, -0.3, 2.0, cs.MMU)
        X0 = E + M + cs.p4(60.0, 0.1, -2.0, 10.0) + cs.p4(45.0, 1.0, 1.0, 8.0)
        self.assertAlmostEqual(r["m_llbb"][0], float(cs.mass(X0)), places=6)
        self.assertAlmostEqual(r["cos_star"][0], float(cs.cos_star(E, M, X0)), places=9)
        self.assertAlmostEqual(r["dphi_ll"][0], 2.0, places=6)
        # b-tagged jet first, then the hardest other jet (80, not 40)
        X4 = E + M + cs.p4(50.0, 1.5, -2.5, 7.0) + cs.p4(80.0, -1.0, -1.0, 9.0)
        self.assertAlmostEqual(r["m_llbb"][4], float(cs.mass(X4)), places=6)


@unittest.skipUnless(HAVE_AK, "awkward + numpy on POSIX only")
class TestGen(unittest.TestCase):
    def test_gen_tt(self):
        (t, tb, lp, lm), expect, mtt = _spin_system()
        T, TB, LP, LM = (_to_ptetaphim(x) for x in (t, tb, lp, lm))
        W = (10.0, 0.0, 0.0, 80.4)
        decoy = (99.0, 0.3, 0.3, 172.5)
        ev0 = [(6, 128, -1, T), (-6, 128, -1, TB), (-13, 128, -1, LP), (11, 128, -1, LM), (6, 8192, -1, decoy)]
        ev1 = [(24, 0, -1, W), (-24, 0, -1, W), (6, 8192, -1, T), (-6, 8192, -1, TB),
               (-11, 4096, 0, LP), (13, 4096, 1, LM)]
        ev2 = [(21, 0, -1, W)]
        evs = [ev0, ev1, ev2]
        d = {"GenPart_pdgId": [[p[0] for p in e] for e in evs],
             "GenPart_statusFlags": [[p[1] for p in e] for e in evs],
             "GenPart_genPartIdxMother": [[p[2] for p in e] for e in evs],
             "GenPart_pt": [[p[3][0] for p in e] for e in evs],
             "GenPart_eta": [[p[3][1] for p in e] for e in evs],
             "GenPart_phi": [[p[3][2] for p in e] for e in evs],
             "GenPart_mass": [[p[3][3] for p in e] for e in evs]}
        m, c, has, src = cs.gen_tt(ak.Array(d))
        self.assertEqual(list(has), [True, True, False])
        self.assertEqual(src, {"hard": 1, "fallback": 1, "none": 1})
        self.assertAlmostEqual(m[0], mtt, places=5)
        self.assertAlmostEqual(m[1], mtt, places=5)
        self.assertAlmostEqual(c[0], expect, places=7)
        self.assertAlmostEqual(c[1], expect, places=7)
        self.assertEqual(list(cs.gen_class(m, c, has)), [cs.CLS_THR, cs.CLS_THR, cs.CLS_NOGEN])


if __name__ == "__main__":
    unittest.main()
