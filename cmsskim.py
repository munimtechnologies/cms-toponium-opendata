# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""cmsskim 0.1.0 -- histogram skim for PREREG-CMSTOPO-v0.1 (CMS 2016 G-H open data, e-mu channel).

For one sample (and one part of its file list) it streams the public NanoAOD files over HTTPS,
applies the pre-registered e-mu selection and fills the pre-registered 2D histogram
m(llbb) x cos(phi*) (9 x 4 = 36 bins), plus, for simulation, the systematic variations, and,
for the top-pair simulation, the generator-level class of every event:
  classes 0..39  generator m(tt) < 350 GeV, split in 40 bins of the generator spin angle
                 cos(phi) (lepton directions in their parent top rest frames), so that the
                 fit can build spin-singlet / triplet / unchanged threshold templates;
  class 40       350 <= m(tt) < 400 GeV;  class 41  m(tt) >= 400 GeV;  class 42  no generator tops.
It also stores generator-level (before selection) sums of weights per class, for normalisation.

No fit is done here; cmsfit (written before any data histogram is read) does that.

Usage: python cmsskim.py --sample tt --part 0 --nparts 4 --out-dir /out
"""
import argparse
import json
import math
import os
import re
import sys
import time
import urllib.parse
import urllib.request

try:  # numpy is only needed inside the job; stdlib-only tests must still import this module
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

VERSION = "0.1.0"
API = "https://opendata.cern.ch/api/records/"
GOLDEN_URL = ("https://opendata.cern.ch/record/14220/files/"
              "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")

BR_W_LNU = 0.1086                     # per lepton flavour (PDG average)
BR_2L2NU = (3 * BR_W_LNU) ** 2        # 0.10614..., e/mu/tau
SIGMA_TT_PB = 833.9                   # total top-pair cross section used by CMS (mt = 172.5 GeV)

# sample -> records and cross sections (pb) per record (None for data)
SAMPLES = {
    "tt": {"kind": "tt", "records": {67801: SIGMA_TT_PB * BR_2L2NU}},
    "tw": {"kind": "bkg", "records": {64895: 19.47, 64839: 19.47}},
    "dy": {"kind": "bkg", "records": {35669: 6077.22}},
    "data_G": {"kind": "data", "records": {30528: None}},
    "data_H": {"kind": "data", "records": {30561: None}},
}

# ---- pre-registered binning (PREREG-CMSTOPO-v0.1, section 4)
M_EDGES = [0.0, 160.0, 200.0, 240.0, 280.0, 330.0, 400.0, 500.0, 650.0, float("inf")]
COS_EDGES = [-1.0, -0.5, 0.0, 0.5, 1.0]
NM = len(M_EDGES) - 1
NCOS = len(COS_EDGES) - 1
NRECO = NM * NCOS
DPHI_EDGES = [math.pi * i / 10.0 for i in range(11)]
NDPHI = 10
NWIN = 40
WIN_EDGES = [-1.0 + 2.0 * i / NWIN for i in range(NWIN + 1)]
M_WIN = 350.0
M_THR = 400.0
CLS_THR, CLS_REST, CLS_NOGEN = NWIN, NWIN + 1, NWIN + 2
NCLS_TT = NWIN + 3

# ---- pre-registered selection constants
BTAG_WP = 0.2489       # DeepJet medium, UL2016 post-VFP
JES = 0.03             # jet-energy-scale proxy variation
ME, MMU = 0.000511, 0.10566
TRIGGERS = ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"]
FLAGS = ["Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
         "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
         "Flag_BadPFMuonFilter", "Flag_BadPFMuonDzFilter", "Flag_eeBadScFilter"]
RECO_BRANCHES = ["PV_npvsGood",
                 "Electron_pt", "Electron_eta", "Electron_phi", "Electron_charge",
                 "Electron_cutBased", "Electron_deltaEtaSC",
                 "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_tightId", "Muon_pfRelIso04_all",
                 "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_jetId", "Jet_btagDeepFlavB"]
DATA_BRANCHES = ["run", "luminosityBlock"]
MC_BRANCHES = ["genWeight"]
TT_BRANCHES = ["LHEScaleWeight", "PSWeight", "GenPart_pdgId", "GenPart_statusFlags",
               "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"]
TT_VARS = ["nom", "mur_up", "mur_dn", "muf_up", "muf_dn", "isr_up", "isr_dn", "fsr_up", "fsr_dn",
           "jes_up", "jes_dn"]
BKG_VARS = ["nom", "jes_up", "jes_dn"]
DATA_VARS = ["nom"]


# ======================================================================= stdlib helpers (tested)
def extract_root_uris(obj, out=None):
    """All CMS open-data ROOT file URIs in a JSON object, in order, unique."""
    if out is None:
        out = []
    if isinstance(obj, dict):
        for v in obj.values():
            extract_root_uris(v, out)
    elif isinstance(obj, list):
        for v in obj:
            extract_root_uris(v, out)
    elif isinstance(obj, str):
        if obj.startswith("root://eospublic.cern.ch//eos/opendata/") and obj.endswith(".root"):
            if obj not in out:
                out.append(obj)
    return out


def url_forms(uri):
    """HTTPS via the portal, HTTPS straight to EOS, then the XRootD URI."""
    prefix = "root://eospublic.cern.ch/"
    if not uri.startswith(prefix):
        raise ValueError("not an eospublic URI: %r" % uri)
    path = "/" + uri[len(prefix):].lstrip("/")
    return ["https://opendata.cern.ch" + path, "https://eospublic.cern.ch:443/" + path, uri]


def split_files(uris, part, nparts, max_files=0):
    """Sorted URIs, every nparts-th starting at part; optionally only the first max_files of those."""
    if not 0 <= part < nparts:
        raise ValueError("part %d not in [0, %d)" % (part, nparts))
    sel = [u for i, u in enumerate(sorted(uris)) if i % nparts == part]
    return sel[:max_files] if max_files and max_files > 0 else sel


def parse_ps_doc(doc):
    """Map isr_up/isr_dn/fsr_up/fsr_dn to PSWeight indices from the branch docstring."""
    found = {}
    for m in re.finditer(r"\[(\d+)\]\s*is\s*ISR\s*=\s*([0-9.]+)\s*FSR\s*=\s*([0-9.]+)", doc or ""):
        found[(float(m.group(2)), float(m.group(3)))] = int(m.group(1))
    need = {"isr_up": (2.0, 1.0), "isr_dn": (0.5, 1.0), "fsr_up": (1.0, 2.0), "fsr_dn": (1.0, 0.5)}
    out = {}
    for k, key in need.items():
        if key not in found:
            raise ValueError("PSWeight docstring lacks ISR=%g FSR=%g: %r" % (key[0], key[1], doc))
        out[k] = found[key]
    return out


def parse_scale_doc(doc):
    """Map mur_up/mur_dn/muf_up/muf_dn to LHEScaleWeight indices from the branch docstring.
    Understands 'renscfact=2d0 facscfact=1d0' and 'MUR="2.0" MUF="1.0"' styles."""
    found = {}
    for seg in re.split(r"(?=\[\d+\]\s*is)", doc or ""):
        m = re.match(r"\[(\d+)\]\s*is(.*)", seg, re.S)
        if not m:
            continue
        body = m.group(2)
        r = re.search(r"(?:renscfact|MUR|muR)\s*=\s*\"?([0-9.]+)", body)
        f = re.search(r"(?:facscfact|MUF|muF)\s*=\s*\"?([0-9.]+)", body)
        if r and f:
            found[(float(r.group(1).rstrip(".")), float(f.group(1).rstrip(".")))] = int(m.group(1))
    need = {"mur_up": (2.0, 1.0), "mur_dn": (0.5, 1.0), "muf_up": (1.0, 2.0), "muf_dn": (1.0, 0.5)}
    out = {}
    for k, key in need.items():
        if key not in found:
            raise ValueError("LHEScaleWeight docstring lacks muR=%g muF=%g: %r" % (key[0], key[1], doc))
        out[k] = found[key]
    return out


# ======================================================================= numpy helpers (tested)
def p4(pt, eta, phi, m):
    """(E, px, py, pz) stacked on the last axis."""
    pt = np.asarray(pt, dtype=float)
    eta = np.asarray(eta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    e = np.sqrt(px * px + py * py + pz * pz + np.asarray(m, dtype=float) ** 2)
    return np.stack([e, px, py, pz], axis=-1)


def mass(p):
    m2 = p[..., 0] ** 2 - (p[..., 1:] ** 2).sum(axis=-1)
    return np.sqrt(np.maximum(m2, 0.0))


def boost(p, b):
    """Boost 4-vectors p by velocity b: a particle at rest ends up moving with velocity b."""
    b = np.asarray(b, dtype=float)
    b2 = (b * b).sum(axis=-1)
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = (b * p[..., 1:]).sum(axis=-1)
    safe = np.where(b2 > 0, b2, 1.0)
    coef = np.where(b2 > 0, (g - 1.0) * bp / safe, 0.0) + g * p[..., 0]
    e = g * (p[..., 0] + bp)
    v = p[..., 1:] + coef[..., None] * b
    return np.concatenate([e[..., None], v], axis=-1)


def to_rest(p, P):
    """p seen in the rest frame of P."""
    return boost(p, -P[..., 1:] / P[..., 0:1])


def unit3(p):
    v = p[..., 1:]
    n = np.sqrt((v * v).sum(axis=-1))
    return v / np.where(n > 0, n, 1.0)[..., None]


def cos_phi_gen(t, tb, lp, lm):
    """Spin-correlation angle: l+ in the top rest frame, l- in the antitop rest frame, both
    reached through the top-pair rest frame (helicity-basis convention of ATLAS and CMS)."""
    P = t + tb
    t1, tb1, lp1, lm1 = to_rest(t, P), to_rest(tb, P), to_rest(lp, P), to_rest(lm, P)
    a = unit3(to_rest(lp1, t1))
    b = unit3(to_rest(lm1, tb1))
    return np.clip((a * b).sum(axis=-1), -1.0, 1.0)


def cos_star(l1, l2, X):
    """Cosine of the opening angle of l1 and l2 in the rest frame of X."""
    a = unit3(to_rest(l1, X))
    b = unit3(to_rest(l2, X))
    return np.clip((a * b).sum(axis=-1), -1.0, 1.0)


def dphi(phi1, phi2):
    d = np.mod(np.asarray(phi1, dtype=float) - np.asarray(phi2, dtype=float) + math.pi, 2 * math.pi) - math.pi
    return np.abs(d)


def bin_index(x, edges):
    i = np.searchsorted(np.asarray(edges, dtype=float), np.asarray(x, dtype=float), side="right") - 1
    return np.clip(i, 0, len(edges) - 2)


def reco_index(m_llbb, cos_s):
    return bin_index(m_llbb, M_EDGES) * NCOS + bin_index(cos_s, COS_EDGES)


def gen_class(mtt, cphi, has_gen):
    mtt = np.where(has_gen, mtt, 0.0)
    cphi = np.where(has_gen, cphi, 0.0)
    cls = np.where(mtt < M_WIN, bin_index(cphi, WIN_EDGES), np.where(mtt < M_THR, CLS_THR, CLS_REST))
    return np.where(has_gen, cls, CLS_NOGEN).astype(int)


def lumi_mask(runs, lumis, golden):
    runs = np.asarray(runs)
    lumis = np.asarray(lumis)
    ok = np.zeros(len(runs), dtype=bool)
    for r in np.unique(runs):
        rng = golden.get(str(int(r)))
        if not rng:
            continue
        sel = runs == r
        lu = lumis[sel]
        m = np.zeros(len(lu), dtype=bool)
        for a, b in rng:
            m |= (lu >= a) & (lu <= b)
        ok[sel] = m
    return ok


# ======================================================================= awkward selection
def _first(x, ok, fill):
    import awkward as ak
    return np.asarray(ak.to_numpy(ak.fill_none(ak.firsts(x[ok]), fill)), dtype=float)


def _count(ok):
    import awkward as ak
    return np.asarray(ak.to_numpy(ak.sum(ok, axis=1)), dtype=int)


def select_reco(ev, jes=1.0, triggers=(), flags=()):
    """Pre-registered e-mu selection. Returns dict of numpy arrays over all events: mask and,
    where mask is True, m_llbb, cos_star, dphi_ll (NaN elsewhere)."""
    import awkward as ak
    n = len(ev)
    sc = abs(ev["Electron_eta"] + ev["Electron_deltaEtaSC"])
    e_ok = ((ev["Electron_cutBased"] >= 4) & (ev["Electron_pt"] > 20) & (abs(ev["Electron_eta"]) < 2.4)
            & ((sc <= 1.4442) | (sc >= 1.5660)))
    m_ok = ((ev["Muon_tightId"]) & (ev["Muon_pfRelIso04_all"] < 0.15) & (ev["Muon_pt"] > 20)
            & (abs(ev["Muon_eta"]) < 2.4))
    ne, nmu = _count(e_ok), _count(m_ok)
    e_pt, e_eta, e_phi = (_first(ev["Electron_pt"], e_ok, 0.0), _first(ev["Electron_eta"], e_ok, 99.0),
                          _first(ev["Electron_phi"], e_ok, 0.0))
    e_q = _first(ev["Electron_charge"], e_ok, 0.0)
    m_pt, m_eta, m_phi = (_first(ev["Muon_pt"], m_ok, 0.0), _first(ev["Muon_eta"], m_ok, 99.0),
                          _first(ev["Muon_phi"], m_ok, 0.0))
    m_q = _first(ev["Muon_charge"], m_ok, 0.0)

    jpt = ev["Jet_pt"] * jes
    jm = ev["Jet_mass"] * jes
    jeta, jphi = ev["Jet_eta"], ev["Jet_phi"]

    def dr(eta0, phi0):
        dp = (jphi - phi0 + math.pi) % (2 * math.pi) - math.pi
        return np.sqrt((jeta - eta0) ** 2 + dp ** 2)

    j_ok = (((ev["Jet_jetId"] & 2) > 0) & (jpt > 30) & (abs(jeta) < 2.4)
            & (dr(e_eta, e_phi) > 0.4) & (dr(m_eta, m_phi) > 0.4))
    isb = j_ok & (ev["Jet_btagDeepFlavB"] > BTAG_WP)
    nj, nb = _count(j_ok), _count(isb)
    key = (ak.where(isb, 1.0e5, 0.0) + jpt)[j_ok]
    order = ak.argsort(key, axis=1, ascending=False)

    def top2(x):
        y = ak.pad_none(x[j_ok][order], 2, axis=1, clip=True)
        return np.asarray(ak.to_numpy(ak.fill_none(y, 0.0)), dtype=float).reshape(n, 2)

    j_pt, j_eta, j_phi, j_m = top2(jpt), top2(jeta), top2(jphi), top2(jm)

    trig = np.zeros(n, dtype=bool)
    for name in triggers:
        trig |= np.asarray(ak.to_numpy(ev[name]), dtype=bool)
    flg = np.ones(n, dtype=bool)
    for name in flags:
        flg &= np.asarray(ak.to_numpy(ev[name]), dtype=bool)
    pv = np.asarray(ak.to_numpy(ev["PV_npvsGood"]), dtype=int) >= 1

    memu2 = 2.0 * e_pt * m_pt * (np.cosh(np.clip(e_eta - m_eta, -50, 50)) - np.cos(e_phi - m_phi))
    memu = np.sqrt(np.maximum(memu2, 0.0))
    mask = ((ne == 1) & (nmu == 1) & (e_q * m_q == -1) & (np.maximum(e_pt, m_pt) > 25) & (memu > 20)
            & (nj >= 2) & (nb >= 1) & trig & flg & pv)
    out = {"mask": mask, "m_llbb": np.full(n, np.nan), "cos_star": np.full(n, np.nan),
           "dphi_ll": np.full(n, np.nan)}
    if mask.any():
        E = p4(e_pt[mask], e_eta[mask], e_phi[mask], ME)
        M = p4(m_pt[mask], m_eta[mask], m_phi[mask], MMU)
        J1 = p4(j_pt[mask, 0], j_eta[mask, 0], j_phi[mask, 0], j_m[mask, 0])
        J2 = p4(j_pt[mask, 1], j_eta[mask, 1], j_phi[mask, 1], j_m[mask, 1])
        X = E + M + J1 + J2
        out["m_llbb"][mask] = mass(X)
        out["cos_star"][mask] = cos_star(E, M, X)
        out["dphi_ll"][mask] = dphi(e_phi[mask], m_phi[mask])
    return out


def gen_tt(ev):
    """Generator top, antitop, l+ and l- (from the W decays; e, mu or tau). Hard-process copies
    first; per-event fallback to last-copy tops and first-copy leptons whose mother is a W.
    Returns mtt, cos_phi, has_gen (numpy) and counts of the source used."""
    import awkward as ak
    n = len(ev)
    pid = ev["GenPart_pdgId"]
    fl = ev["GenPart_statusFlags"]
    mom = ev["GenPart_genPartIdxMother"]
    hard = (fl & 128) > 0
    first_copy = (fl & 4096) > 0
    last_copy = (fl & 8192) > 0
    apid = abs(pid)
    is_l = (apid == 11) | (apid == 13) | (apid == 15)
    safe = ak.where(mom >= 0, mom, 0)
    from_w = (mom >= 0) & (abs(pid[safe]) == 24)
    sets = {
        "hard": (hard & (pid == 6), hard & (pid == -6), hard & is_l & (pid < 0), hard & is_l & (pid > 0)),
        "fallback": (last_copy & (pid == 6), last_copy & (pid == -6),
                     first_copy & from_w & is_l & (pid < 0), first_copy & from_w & is_l & (pid > 0)),
    }
    res = {}
    for key, sels in sets.items():
        ok = np.ones(n, dtype=bool)
        vecs = []
        for s in sels:
            ok &= _count(s) >= 1
            vecs.append(p4(_first(ev["GenPart_pt"], s, 0.0), _first(ev["GenPart_eta"], s, 0.0),
                           _first(ev["GenPart_phi"], s, 0.0), _first(ev["GenPart_mass"], s, 0.0)))
        res[key] = (ok, vecs)
    ok_h, vh = res["hard"]
    ok_f, vf = res["fallback"]
    use_f = (~ok_h) & ok_f
    has = ok_h | use_f
    vec = [np.where(ok_h[:, None], a, b) for a, b in zip(vh, vf)]
    mtt = np.full(n, np.nan)
    cphi = np.full(n, np.nan)
    if has.any():
        t, tb, lp, lm = (v[has] for v in vec)
        mtt[has] = mass(t + tb)
        cphi[has] = cos_phi_gen(t, tb, lp, lm)
    src = {"hard": int(ok_h.sum()), "fallback": int(use_f.sum()), "none": int((~has).sum())}
    return mtt, cphi, has, src


def weight_column(x, i):
    import awkward as ak
    y = ak.pad_none(x, i + 1, axis=1)
    return np.asarray(ak.to_numpy(ak.fill_none(y[:, i], 1.0)), dtype=float)


# ======================================================================= accumulation
class Acc:
    def __init__(self, ncls, variations):
        self.ncls = ncls
        self.variations = list(variations)
        self.h = {v: np.zeros((ncls, NRECO)) for v in self.variations}
        self.h2 = np.zeros((ncls, NRECO))
        self.n = np.zeros((ncls, NRECO))
        self.dphi = np.zeros((ncls, NDPHI))
        self.gen = np.zeros(ncls)
        self.n_entries = 0
        self.sumw = 0.0
        self.sumw2 = 0.0
        self.sumw_runs = 0.0
        self.nsel = {v: 0 for v in self.variations}
        self.gen_source = {"hard": 0, "fallback": 0, "none": 0}

    def add(self, o):
        for v in self.variations:
            self.h[v] += o.h[v]
            self.nsel[v] += o.nsel[v]
        self.h2 += o.h2
        self.n += o.n
        self.dphi += o.dphi
        self.gen += o.gen
        self.n_entries += o.n_entries
        self.sumw += o.sumw
        self.sumw2 += o.sumw2
        self.sumw_runs += o.sumw_runs
        for k in self.gen_source:
            self.gen_source[k] += o.gen_source[k]

    def to_json(self):
        return {"ncls": self.ncls, "variations": self.variations,
                "h": {v: self.h[v].tolist() for v in self.variations}, "h2_nom": self.h2.tolist(),
                "n_nom": self.n.tolist(), "dphi_nom": self.dphi.tolist(), "gen_sumw": self.gen.tolist(),
                "n_entries": self.n_entries, "sumw_events": self.sumw, "sumw2_events": self.sumw2,
                "sumw_runs": self.sumw_runs, "nsel": self.nsel, "gen_source": self.gen_source}


def process_chunk(ev, kind, ctx, acc):
    import awkward as ak
    n = len(ev)
    acc.n_entries += n
    if kind == "data":
        w = np.ones(n)
        good = lumi_mask(ak.to_numpy(ev["run"]), ak.to_numpy(ev["luminosityBlock"]), ctx["golden"])
    else:
        w = np.asarray(ak.to_numpy(ev["genWeight"]), dtype=float)
        good = np.ones(n, dtype=bool)
        acc.sumw += float(w.sum())
        acc.sumw2 += float((w * w).sum())
    if kind == "tt":
        mtt, cphi, has, src = gen_tt(ev)
        cls = gen_class(mtt, cphi, has)
        np.add.at(acc.gen, cls, w)
        for k in src:
            acc.gen_source[k] += src[k]
    else:
        cls = np.zeros(n, dtype=int)
    wv = {"nom": w}
    if kind == "tt":
        for k, i in ctx["scale_idx"].items():
            wv[k] = w * weight_column(ev["LHEScaleWeight"], i)
        for k, i in ctx["ps_idx"].items():
            wv[k] = w * weight_column(ev["PSWeight"], i)
    jes_list = [("nom", 1.0)]
    if kind != "data":
        jes_list += [("jes_up", 1.0 + JES), ("jes_dn", 1.0 - JES)]
    for name, s in jes_list:
        r = select_reco(ev, s, ctx["triggers"], ctx["flags"])
        m = r["mask"] & good
        if not m.any():
            continue
        ridx = reco_index(r["m_llbb"][m], r["cos_star"][m])
        c = cls[m]
        if name == "nom":
            for v, wa in wv.items():
                np.add.at(acc.h[v], (c, ridx), wa[m])
                acc.nsel[v] += int(m.sum())
            np.add.at(acc.h2, (c, ridx), w[m] ** 2)
            np.add.at(acc.n, (c, ridx), 1.0)
            np.add.at(acc.dphi, (c, bin_index(r["dphi_ll"][m], DPHI_EDGES)), w[m])
        else:
            np.add.at(acc.h[name], (c, ridx), w[m])
            acc.nsel[name] += int(m.sum())


# ======================================================================= I/O
def get_json(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cmsskim/0.1 (open-data research)"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(10 * (i + 1))
    raise RuntimeError("GET failed %s: %s" % (url, last))


def file_uris(recid, rec):
    """ROOT URIs of a record; if the record JSON does not list them, read its file-index JSONs."""
    uris = extract_root_uris(rec)
    if uris:
        return uris
    meta = rec.get("metadata") or {}
    for idx in meta.get("_file_indices") or []:
        key = idx.get("key") or ""
        if key.endswith(".json"):
            try:
                for u in extract_root_uris(get_json("https://opendata.cern.ch/record/%d/files/%s" % (recid, key))):
                    if u not in uris:
                        uris.append(u)
            except Exception:  # noqa: BLE001
                pass
    return uris


def log(fh, line):
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def open_file(uri):
    import uproot
    errs = []
    for url in url_forms(uri):
        handlers = [None]
        if url.startswith("https://"):
            handlers = [None, uproot.source.http.HTTPSource]
        for h in handlers:
            try:
                kw = {"timeout": 300}
                if h is not None:
                    kw["handler"] = h
                f = uproot.open(url, **kw)
                f["Events"]
                return f, url, errs
            except Exception as exc:  # noqa: BLE001
                errs.append("%s [%s]: %s" % (url, getattr(h, "__name__", "default"), repr(exc)[:200]))
    return None, None, errs


def runs_sumw(f):
    try:
        runs = f["Runs"]
        for name in ("genEventSumw", "genEventSumw_"):
            if name in runs.keys():
                return float(np.sum(runs[name].array(library="np")))
    except Exception:  # noqa: BLE001
        pass
    return float("nan")


def process_file(uri, kind, ctx, ncls, variations, step, max_entries, fh):
    f, url, errs = open_file(uri)
    if f is None:
        raise RuntimeError("cannot open: " + " | ".join(errs)[:600])
    tree = f["Events"]
    names = set(tree.keys())
    ctx["triggers"] = [t for t in TRIGGERS if t in names]
    ctx["flags"] = [x for x in FLAGS if x in names]
    if not ctx["triggers"]:
        raise RuntimeError("no e-mu trigger branch present")
    branches = list(RECO_BRANCHES) + ctx["triggers"] + ctx["flags"]
    if kind == "data":
        branches += DATA_BRANCHES
    else:
        branches += MC_BRANCHES
    if kind == "tt":
        branches += TT_BRANCHES
        ps_doc, sc_doc = tree["PSWeight"].title, tree["LHEScaleWeight"].title
        ctx["ps_idx"] = parse_ps_doc(ps_doc)
        ctx["scale_idx"] = parse_scale_doc(sc_doc)
        ctx.setdefault("docs", {"PSWeight": ps_doc, "LHEScaleWeight": sc_doc})
    miss = [b for b in branches if b not in names]
    if miss:
        raise RuntimeError("missing branches: %s" % ",".join(miss))
    acc = Acc(ncls, variations)
    if kind != "data":
        acc.sumw_runs = runs_sumw(f)
    stop = min(max_entries, tree.num_entries) if max_entries else tree.num_entries
    first = True
    for ev in tree.iterate(branches, step_size=step, entry_stop=stop, library="ak"):
        process_chunk(ev, kind, ctx, acc)
        if first and kind == "tt":
            frac = 1.0 - acc.gen_source["none"] / max(1, acc.n_entries)
            log(fh, "CMSSKIM_GENCHECK has_gen=%.4f hard=%d fallback=%d none=%d" % (
                frac, acc.gen_source["hard"], acc.gen_source["fallback"], acc.gen_source["none"]))
            if frac < 0.95:
                raise SystemExit(5)
        first = False
    return acc, url


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, choices=sorted(SAMPLES))
    ap.add_argument("--part", type=int, default=0)
    ap.add_argument("--nparts", type=int, default=1)
    ap.add_argument("--max-files", type=int, default=0)
    ap.add_argument("--max-minutes", type=float, default=300.0)
    ap.add_argument("--step", type=int, default=25000)
    ap.add_argument("--max-entries-per-file", type=int, default=0)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    t_start = time.time()
    tag = "%s-p%d" % (args.sample, args.part)
    fh = open(os.path.join(args.out_dir, "skim-%s.txt" % tag), "w")
    out_path = os.path.join(args.out_dir, "skim-%s.json" % tag)
    import uproot
    import awkward
    log(fh, "CMSSKIM_START version=%s sample=%s part=%d/%d python=%s numpy=%s uproot=%s awkward=%s" % (
        VERSION, args.sample, args.part, args.nparts, sys.version.split()[0], np.__version__,
        uproot.__version__, awkward.__version__))
    spec = SAMPLES[args.sample]
    kind = spec["kind"]
    variations = TT_VARS if kind == "tt" else (BKG_VARS if kind == "bkg" else DATA_VARS)
    ncls = NCLS_TT if kind == "tt" else 1
    ctx = {}
    if kind == "data":
        ctx["golden"] = get_json(GOLDEN_URL)
        log(fh, "CMSSKIM_GOLDEN runs=%d" % len(ctx["golden"]))
    result = {"version": VERSION, "sample": args.sample, "kind": kind, "part": args.part,
              "nparts": args.nparts, "max_files": args.max_files,
              "max_entries_per_file": args.max_entries_per_file,
              "constants": {"M_EDGES": [x if math.isfinite(x) else "inf" for x in M_EDGES],
                            "COS_EDGES": COS_EDGES, "NWIN": NWIN, "M_WIN": M_WIN, "M_THR": M_THR,
                            "BTAG_WP": BTAG_WP, "JES": JES, "BR_2L2NU": BR_2L2NU,
                            "SIGMA_TT_PB": SIGMA_TT_PB},
              "records": {}}
    n_ok = 0
    for recid, xsec in spec["records"].items():
        rec = get_json(API + str(recid))
        uris = split_files(file_uris(recid, rec), args.part, args.nparts, args.max_files)
        title = (rec.get("metadata") or {}).get("title")
        log(fh, "CMSSKIM_RECORD %d %s files_assigned=%d" % (recid, title, len(uris)))
        racc = Acc(ncls, variations)
        rinfo = {"title": title, "xsec_pb": xsec, "files_assigned": uris, "files_done": [],
                 "files_failed": [], "files_skipped_time": []}
        for uri in uris:
            if (time.time() - t_start) / 60.0 > args.max_minutes:
                rinfo["files_skipped_time"].append(uri)
                log(fh, "CMSSKIM_FILE SKIP_TIME %s" % uri)
                continue
            done = False
            for attempt in range(3):
                t0 = time.time()
                try:
                    acc, url = process_file(uri, kind, ctx, ncls, variations, args.step,
                                            args.max_entries_per_file, fh)
                except SystemExit:
                    log(fh, "CMSSKIM_GENFAIL %s" % uri)
                    raise
                except Exception as exc:  # noqa: BLE001
                    msg = repr(exc)[:400]
                    log(fh, "CMSSKIM_FILE RETRY %d %s %s" % (attempt, uri, msg))
                    if "docstring lacks" in msg or "missing branches" in msg:
                        log(fh, "CMSSKIM_DOCFAIL %s" % msg)
                        sys.exit(6)
                    time.sleep(30 * (attempt + 1))
                    continue
                racc.add(acc)
                rinfo["files_done"].append({"uri": uri, "url": url, "entries": acc.n_entries,
                                            "sumw_events": acc.sumw, "sumw_runs": acc.sumw_runs,
                                            "seconds": round(time.time() - t0, 1)})
                log(fh, "CMSSKIM_FILE OK %s entries=%d sel=%d sumw=%.6g sumw_runs=%.6g %.0fs" % (
                    uri, acc.n_entries, acc.nsel["nom"], acc.sumw, acc.sumw_runs, time.time() - t0))
                done = True
                n_ok += 1
                break
            if not done:
                rinfo["files_failed"].append(uri)
                log(fh, "CMSSKIM_FILE FAIL %s" % uri)
            rinfo.update(racc.to_json())
            rinfo["docs"] = ctx.get("docs")
            rinfo["ps_idx"] = ctx.get("ps_idx")
            rinfo["scale_idx"] = ctx.get("scale_idx")
            rinfo["triggers_used"] = ctx.get("triggers")
            rinfo["flags_used"] = ctx.get("flags")
            result["records"][str(recid)] = rinfo
            with open(out_path + ".tmp", "w") as jf:
                json.dump(result, jf)
            os.replace(out_path + ".tmp", out_path)
        rinfo.update(racc.to_json())
        result["records"][str(recid)] = rinfo
        log(fh, "CMSSKIM_RECORD_DONE %d done=%d failed=%d skipped_time=%d entries=%d sel_nom=%d" % (
            recid, len(rinfo["files_done"]), len(rinfo["files_failed"]), len(rinfo["files_skipped_time"]),
            racc.n_entries, racc.nsel["nom"]))
    result["minutes"] = round((time.time() - t_start) / 60.0, 1)
    with open(out_path + ".tmp", "w") as jf:
        json.dump(result, jf)
    os.replace(out_path + ".tmp", out_path)
    log(fh, "CMSSKIM_DONE files_ok=%d minutes=%.1f" % (n_ok, result["minutes"]))
    fh.close()
    return 0 if n_ok > 0 else 3


if __name__ == "__main__":
    sys.exit(main())
