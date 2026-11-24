# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""cmsfit 0.1.0 -- the pre-registered likelihood fit of PREREG-CMSTOPO-v0.1, section 6
(with the implementation details fixed in AMENDMENT-CMSTOPO-0.1.1).

Input: skim JSON files written by cmsskim 0.1.0 (one per sample and part).
Modes:
  selftest  synthetic skims made here -> the whole chain (templates, Asimov, data-like fit);
            used to check the machinery, never graded.
  asimov    simulation skims only -> gate G1, expected sensitivity Z_A at 8.8 and 6.4 pb,
            gate G4 on the Asimov fits, grade H4. Refuses to read any data skim (blindness).
  data      simulation + data skims -> gates G1-G4, grades H1-H3 (+ H4 and information-only fits).

Expected count in fit bin i (36 bins = 9 m(llbb) x 4 cos(phi*)):
  nu_i = gamma_i * [ k_tt * ( T_i(theta) + mu * S_i(theta_JES) )
                     + 1.15^theta_tW * TW_i(theta_JES) + 1.5^theta_oth * O_i(theta_JES) ]
Usage:
  python cmsfit.py --mode selftest --out-dir /out
  python cmsfit.py --mode asimov --skim-dir skims --out-dir /out [--fetch manifest.json]
"""
import argparse
import glob
import hashlib
import json
import math
import os
import sys
import time
import urllib.request

try:  # numpy/scipy are only needed inside the job; stdlib-only tests must still import this module
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

VERSION = "0.1.0"
SKIM_VERSION = "0.1.0"
REPO_RAW = "https://raw.githubusercontent.com/munimtechnologies/cms-toponium-opendata/"

# ---- pre-registered constants (PREREG-CMSTOPO-v0.1 sections 2, 4, 5, 6, 8)
L_PB = 16393.4
BR_W_LNU = 0.1086
BR_2L2NU = (3 * BR_W_LNU) ** 2          # 0.10614564
SIGMA_TT_PB = 833.9
XS_TT2L = SIGMA_TT_PB * BR_2L2NU        # 88.51 pb
XS_TW = 19.47
XS_DY = 6077.22
REC_TT, REC_TW, REC_TWB, REC_DY, REC_G, REC_H = "67801", "64895", "64839", "35669", "30528", "30561"
M_EDGES = [0.0, 160.0, 200.0, 240.0, 280.0, 330.0, 400.0, 500.0, 650.0, "inf"]
COS_EDGES = [-1.0, -0.5, 0.0, 0.5, 1.0]
NM, NCOS = 9, 4
NBIN = NM * NCOS
NWIN = 40
M_WIN, M_THR = 350.0, 400.0
NCLS_TT = NWIN + 3                      # 0..39 window, 40 = 350-400, 41 = >=400, 42 = no generator tops
N_LOW = NWIN + 1                        # classes with generator m(tt) < 400 GeV (theta_thr acts on them)
BTAG_WP, JES = 0.2489, 0.03
TT_VARS = ["nom", "mur_up", "mur_dn", "muf_up", "muf_dn", "isr_up", "isr_dn", "fsr_up", "fsr_dn",
           "jes_up", "jes_dn"]
BKG_VARS = ["nom", "jes_up", "jes_dn"]
SHAPE_ONLY = ["mur", "muf", "isr", "fsr"]

MU_CMS, MU_NRQCD, CMS_SYST = 8.8, 6.4, 1.3
MU_MAX = 60.0
KTT_BOUNDS = (0.3, 3.0)
THETA_BOUND = 5.0                       # AMENDMENT-CMSTOPO-0.1.1 (a)
THR_SIZE, TW_SIZE, OTH_SIZE = 0.03, 1.15, 1.5
MU_STARTS = [0.0, 4.0, 8.8, 15.0, 25.0]
KTT_STARTS_B = [0.8, 0.9, 1.0, 1.1, 1.2]  # AMENDMENT-CMSTOPO-0.1.1 (b): background-only starts
SCAN_STEP, SCAN_MAX = 0.1, 40.0
YIELD_FLOOR = 1e-9
G4_TOL = 0.01
TARGETS = {"singlet": -1.0, "triplet": 1.0 / 3.0, "unchanged": None}

PARS = ["mu", "k_tt", "th_mur", "th_muf", "th_isr", "th_fsr", "th_jes", "th_thr", "th_tw", "th_oth"]
IDX = {p: i for i, p in enumerate(PARS)}
NPAR = len(PARS)


# ======================================================================= stdlib helpers (tested)
def git_blob_sha1(data):
    """Git's blob hash of a byte string (what GitHub trees list as 'sha')."""
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def class_centres(n=NWIN):
    return [-1.0 + 2.0 * (k + 0.5) / n for k in range(n)]


def d_win(G):
    """Spin-correlation coefficient D of a window: the distribution (1 - D c)/2 has mean -D/3."""
    c = class_centres(len(G))
    tot = float(sum(G))
    if tot <= 0:
        raise ValueError("empty window")
    return -3.0 * sum(g * ck for g, ck in zip(G, c)) / tot


def spin_weights(G, d_tgt):
    """w_k = (1 - D_tgt c_k) / (1 - D_win c_k); all ones for the unchanged template (d_tgt None)."""
    if d_tgt is None:
        return [1.0] * len(G)
    dw = d_win(G)
    out = []
    for ck in class_centres(len(G)):
        den = 1.0 - dw * ck
        num = 1.0 - d_tgt * ck
        if den <= 0 or num < 0:
            raise ValueError("spin reweighting undefined: D_win=%r D_tgt=%r" % (dw, d_tgt))
        out.append(num / den)
    return out


def gamma_bb(nu, n, delta):
    """Barlow-Beeston-lite gamma for one bin, profiled analytically (PREREG section 6)."""
    if delta <= 0:
        return 1.0
    d2 = delta * delta
    b = nu * d2 - 1.0
    sq = math.sqrt(b * b + 4.0 * n * d2)
    return 2.0 * n * d2 / (b + sq) if b > 0 else 0.5 * (-b + sq)


def grade_h4(z):
    return "HIGH" if z >= 3.0 else ("MEDIUM" if z >= 1.5 else "LOW")


def grade_h1(z):
    return "EVIDENCE" if z >= 3.0 else ("HINT" if z >= 2.0 else "NOT SEEN")


def grade_h2(mu_hat, lo, hi):
    if lo is None or hi is None:
        return "NOT GRADABLE"
    s = 0.5 * (hi - lo)
    if abs(mu_hat - MU_CMS) <= 2.0 * math.sqrt(s * s + CMS_SYST * CMS_SYST):
        return "CONSISTENT"
    return "HIGHER" if mu_hat > MU_CMS else "LOWER"


def grade_h3(nll_singlet, nll_unchanged, nll_triplet):
    du = 2.0 * (nll_unchanged - nll_singlet)
    dt = 2.0 * (nll_triplet - nll_singlet)
    if du >= 1.0 and dt >= 1.0:
        return "SINGLET PREFERRED"
    if du <= -1.0 or dt <= -1.0:
        return "NOT SINGLET"
    return "UNRESOLVED"


def gate_g3(k_tt):
    return 0.7 <= k_tt <= 1.3


def spread(values):
    return max(values) - min(values)


def gate_g4(nlls_sb, nlls_b):
    return spread(nlls_sb) <= G4_TOL and spread(nlls_b) <= G4_TOL


def scan_grid():
    n = int(round(SCAN_MAX / SCAN_STEP))
    return [round(SCAN_STEP * i, 10) for i in range(n + 1)]


def interval(mus, qs, level=1.0):
    """Lowest and highest scanned mu with q = 2 Delta NLL <= level (PREREG: 0.1 pb grid)."""
    ok = [m for m, q in zip(mus, qs) if q <= level]
    if not ok:
        return None, None
    return min(ok), max(ok)


def skim_constants():
    return {"M_EDGES": M_EDGES, "COS_EDGES": COS_EDGES, "NWIN": NWIN, "M_WIN": M_WIN, "M_THR": M_THR,
            "BTAG_WP": BTAG_WP, "JES": JES, "BR_2L2NU": BR_2L2NU, "SIGMA_TT_PB": SIGMA_TT_PB}


def check_constants(c):
    want = skim_constants()
    for k, v in want.items():
        got = c.get(k)
        if isinstance(v, list):
            if not isinstance(got, list) or len(got) != len(v):
                raise ValueError("skim constant %s differs: %r" % (k, got))
            for a, b in zip(got, v):
                if isinstance(b, str) or isinstance(a, str):
                    if str(a) != str(b):
                        raise ValueError("skim constant %s differs: %r" % (k, got))
                elif abs(float(a) - float(b)) > 1e-12:
                    raise ValueError("skim constant %s differs: %r" % (k, got))
        elif got is None or abs(float(got) - float(v)) > 1e-12:
            raise ValueError("skim constant %s differs: %r" % (k, got))


def _add(a, b):
    """Element-wise sum of equally shaped nested lists (None + b = copy of b)."""
    if a is None:
        return json.loads(json.dumps(b))
    if isinstance(b, list):
        if not isinstance(a, list) or len(a) != len(b):
            raise ValueError("shape mismatch while merging skims")
        return [_add(x, y) for x, y in zip(a, b)]
    return float(a) + float(b)


def merge_skims(skims):
    """Sum skim JSONs (all parts of each sample). Returns {sample: {...}}; pure Python."""
    out = {}
    for s in skims:
        if s.get("version") != SKIM_VERSION:
            raise ValueError("skim version %r, need %s" % (s.get("version"), SKIM_VERSION))
        check_constants(s.get("constants") or {})
        samp = s["sample"]
        m = out.setdefault(samp, {"kind": s["kind"], "nparts": s["nparts"], "parts": [],
                                  "part_has_gen": {}, "max_entries_per_file": [], "records": {}})
        if s["nparts"] != m["nparts"] or s["kind"] != m["kind"]:
            raise ValueError("inconsistent parts for sample %s" % samp)
        if s["part"] in m["parts"]:
            raise ValueError("part %d of %s given twice" % (s["part"], samp))
        m["parts"].append(s["part"])
        m["max_entries_per_file"].append(int(s.get("max_entries_per_file") or 0))
        n_ent, n_none = 0, 0
        for rid, r in s["records"].items():
            R = m["records"].setdefault(rid, {"title": r.get("title"), "xsec_pb": r.get("xsec_pb"), "h": {},
                                              "h2": None, "gen": None, "n_entries": 0, "sumw_events": 0.0,
                                              "sumw_runs": 0.0, "files_assigned": [], "files_done": [],
                                              "files_failed": [], "files_skipped_time": []})
            if R["xsec_pb"] != r.get("xsec_pb"):
                raise ValueError("cross section differs between parts of record %s" % rid)
            for v, h in (r.get("h") or {}).items():
                R["h"][v] = _add(R["h"].get(v), h)
            if r.get("h2_nom") is not None:
                R["h2"] = _add(R["h2"], r["h2_nom"])
            if r.get("gen_sumw") is not None:
                R["gen"] = _add(R["gen"], r["gen_sumw"])
            R["n_entries"] += int(r.get("n_entries") or 0)
            R["sumw_events"] += float(r.get("sumw_events") or 0.0)
            R["sumw_runs"] += float(r.get("sumw_runs") if r.get("sumw_runs") is not None else float("nan"))
            for key in ("files_assigned", "files_done", "files_failed", "files_skipped_time"):
                R[key] += list(r.get(key) or [])
            n_ent += int(r.get("n_entries") or 0)
            n_none += int((r.get("gen_source") or {}).get("none", 0))
        if s["kind"] == "tt":
            m["part_has_gen"][s["part"]] = 1.0 - n_none / max(1, n_ent)
    return out


def gate_g1(merged):
    """Simulation gate (PREREG section 7). Returns (pass, reasons)."""
    why = []
    tt = merged.get("tt")
    if not tt or REC_TT not in tt["records"]:
        return False, ["no top-pair simulation"]
    n_tt = tt["records"][REC_TT]["n_entries"]
    if n_tt < 8000000:
        why.append("top-pair events processed %d < 8,000,000" % n_tt)
    for part, f in sorted(tt["part_has_gen"].items()):
        if not f >= 0.95:
            why.append("tt part %s has_gen %.4f < 0.95" % (part, f))
    for samp in ("tt", "tw", "dy"):
        for rid, R in (merged.get(samp) or {"records": {}})["records"].items():
            for fd in R["files_done"]:
                se, sr = float(fd.get("sumw_events", float("nan"))), float(fd.get("sumw_runs", float("nan")))
                rel = abs(se - sr) / sr if sr and sr == sr else float("nan")
                if not rel < 1e-6:
                    why.append("sum genWeight vs genEventSumw off in %s (%r)" % (fd.get("uri"), rel))
        if samp in merged and any(merged[samp]["max_entries_per_file"]):
            why.append("%s skimmed with a per-file entry cap" % samp)
    for samp, rids in (("tw", (REC_TW, REC_TWB)), ("dy", (REC_DY,))):
        for rid in rids:
            R = (merged.get(samp) or {"records": {}})["records"].get(rid)
            if not R or not R["files_done"]:
                why.append("no file processed for record %s" % rid)
    return not why, why


def gate_g2(merged):
    """Data gate: every part present and >= 95 % of each MuonEG record's files processed."""
    why = []
    for samp, rid in (("data_G", REC_G), ("data_H", REC_H)):
        m = merged.get(samp)
        if not m:
            why.append("no %s skim" % samp)
            continue
        if sorted(m["parts"]) != list(range(m["nparts"])):
            why.append("%s parts %s of %d" % (samp, sorted(m["parts"]), m["nparts"]))
        R = m["records"].get(rid)
        if not R:
            why.append("record %s missing" % rid)
            continue
        n_ass, n_done = len(R["files_assigned"]), len(R["files_done"])
        if n_ass == 0 or n_done < 0.95 * n_ass:
            why.append("record %s: %d of %d files processed" % (rid, n_done, n_ass))
    return not why, why


# ======================================================================= templates (numpy)
def _arr(x):
    return np.asarray(x, dtype=float)


def build_templates(merged, lumi=L_PB):
    """Per-bin templates of PREREG section 5 (all in expected events for luminosity lumi)."""
    R = merged["tt"]["records"][REC_TT]
    if abs(float(R["xsec_pb"]) - XS_TT2L) > 1e-6 * XS_TT2L:
        raise ValueError("top-pair cross section %r != %r" % (R["xsec_pb"], XS_TT2L))
    if not R["sumw_runs"] > 0:
        raise ValueError("top-pair genEventSumw missing")
    norm = lumi * XS_TT2L / R["sumw_runs"]
    H = {v: _arr(R["h"][v]) for v in TT_VARS}
    for v in TT_VARS:
        if H[v].shape != (NCLS_TT, NBIN):
            raise ValueError("tt histogram %s has shape %s" % (v, H[v].shape))
    low = {v: norm * H[v][:N_LOW].sum(axis=0) for v in TT_VARS}
    rest = {v: norm * H[v][N_LOW:].sum(axis=0) for v in TT_VARS}
    tot0 = float((low["nom"] + rest["nom"]).sum())
    for name in SHAPE_ONLY:          # shape only: varied template rescaled to the nominal total
        for d in ("_up", "_dn"):
            v = name + d
            r = tot0 / float((low[v] + rest[v]).sum())
            low[v] = low[v] * r
            rest[v] = rest[v] * r
    G = _arr(R["gen"])[:NWIN]
    S, dwin = {}, d_win(G.tolist())
    for tgt, dt in TARGETS.items():
        w = _arr(spin_weights(G.tolist(), dt))
        den = float((w * G).sum())
        S[tgt] = {v: lumi * 1.0 * BR_2L2NU * (w @ H[v][:NWIN]) / den for v in BKG_VARS}
    sumw2 = norm * norm * _arr(R["h2"]).sum(axis=0)
    norms = {REC_TT: norm}

    def bkg(samp, recs, xs):
        nonlocal sumw2
        out = {v: np.zeros(NBIN) for v in BKG_VARS}
        for rid in recs:
            B = merged[samp]["records"][rid]
            if abs(float(B["xsec_pb"]) - xs) > 1e-9 * xs:
                raise ValueError("cross section of %s is %r, want %r" % (rid, B["xsec_pb"], xs))
            if not B["sumw_runs"] > 0:
                raise ValueError("genEventSumw missing for %s" % rid)
            nb = lumi * xs / B["sumw_runs"]
            norms[rid] = nb
            for v in BKG_VARS:
                out[v] = out[v] + nb * _arr(B["h"][v])[0]
            sumw2 = sumw2 + nb * nb * _arr(B["h2"])[0]
        return out

    TW = bkg("tw", (REC_TW, REC_TWB), XS_TW)
    O = bkg("dy", (REC_DY,), XS_DY)
    total = low["nom"] + rest["nom"] + TW["nom"] + O["nom"]
    delta = np.where(total > 0, np.sqrt(sumw2) / np.where(total > 0, total, 1.0), 0.0)
    tpl = {"low": low, "rest": rest, "S": S, "TW": TW, "O": O, "delta": delta, "D_win": dwin,
           "norms": norms, "data": None}
    data = [merged[s]["records"][r] for s, r in (("data_G", REC_G), ("data_H", REC_H)) if s in merged]
    if data:
        tpl["data"] = sum(_arr(d["h"]["nom"])[0] for d in data)
    return tpl


# ======================================================================= model and fit (numpy)
class Model:
    """nu(x) and its Jacobian for parameters PARS; per-bin piecewise-linear interpolation."""

    def __init__(self, tpl, target="singlet", bb=True):
        L, R = tpl["low"], tpl["rest"]
        self.target = target
        self.low0, self.rest0 = L["nom"], R["nom"]
        self.slopes = []
        for p, v in [("th_mur", "mur"), ("th_muf", "muf"), ("th_isr", "isr"), ("th_fsr", "fsr"),
                     ("th_jes", "jes")]:
            self.slopes.append((IDX[p], L[v + "_up"] - L["nom"], L["nom"] - L[v + "_dn"],
                                R[v + "_up"] - R["nom"], R["nom"] - R[v + "_dn"]))

        def sl(t):
            return t["nom"], t["jes_up"] - t["nom"], t["nom"] - t["jes_dn"]
        self.S0, self.S_up, self.S_dn = sl(tpl["S"][target])
        self.TW0, self.TW_up, self.TW_dn = sl(tpl["TW"])
        self.O0, self.O_up, self.O_dn = sl(tpl["O"])
        self.d2 = (tpl["delta"] ** 2) if bb else np.zeros(NBIN)
        self.bb = bb

    def nu_jac(self, x):
        mu, k = x[0], x[1]
        low, rest = self.low0.copy(), self.rest0.copy()
        J = np.zeros((NPAR, NBIN))
        dl, dr = {}, {}
        for i, ul, dn_l, ur, dn_r in self.slopes:
            th = x[i]
            a, b = (ul, ur) if th >= 0 else (dn_l, dn_r)
            low = low + th * a
            rest = rest + th * b
            dl[i], dr[i] = a, b
        fthr = 1.0 + THR_SIZE * x[IDX["th_thr"]]
        T = low * fthr + rest
        thj = x[IDX["th_jes"]]
        sS, sTW, sO = ((self.S_up, self.TW_up, self.O_up) if thj >= 0 else (self.S_dn, self.TW_dn, self.O_dn))
        S = self.S0 + thj * sS
        TW = self.TW0 + thj * sTW
        O = self.O0 + thj * sO
        ftw = TW_SIZE ** x[IDX["th_tw"]]
        fo = OTH_SIZE ** x[IDX["th_oth"]]
        nu = k * (T + mu * S) + ftw * TW + fo * O
        J[0] = k * S
        J[1] = T + mu * S
        for i in dl:
            J[i] = k * (dl[i] * fthr + dr[i])
        J[IDX["th_jes"]] += k * mu * sS + ftw * sTW + fo * sO
        J[IDX["th_thr"]] = k * THR_SIZE * low
        J[IDX["th_tw"]] = math.log(TW_SIZE) * ftw * TW
        J[IDX["th_oth"]] = math.log(OTH_SIZE) * fo * O
        floored = nu < YIELD_FLOOR
        if floored.any():
            nu = np.where(floored, YIELD_FLOOR, nu)
            J[:, floored] = 0.0
        return nu, J

    def gamma(self, nu, n):
        d2 = self.d2
        b = nu * d2 - 1.0
        sq = np.sqrt(b * b + 4.0 * n * d2)
        pos = b > 0
        g = np.where(pos, 2.0 * n * d2 / np.where(pos, b + sq, 1.0), 0.5 * (-b + sq))
        return np.where(d2 > 0, g, 1.0)

    def nll(self, x, n):
        """Profiled NLL (Poisson deviance form, i.e. up to a constant in n) and its gradient."""
        nu, J = self.nu_jac(x)
        g = self.gamma(nu, n)
        t = nu * g
        pos = n > 0
        safe_t = np.where(t > 0, t, YIELD_FLOOR)
        term = np.where(pos, t - n + n * np.log(np.where(pos, n, 1.0) / safe_t), t)
        d2 = self.d2
        pen = np.where(d2 > 0, (g - 1.0) ** 2 / (2.0 * np.where(d2 > 0, d2, 1.0)), 0.0)
        th = x[2:]
        val = float(term.sum() + pen.sum() + 0.5 * float(th @ th))
        dnu = g - n / nu
        grad = J @ dnu
        grad[2:] += th
        return val, grad

    def asimov(self, mu):
        x = np.zeros(NPAR)
        x[0], x[1] = mu, 1.0
        return self.nu_jac(x)[0]


BOUNDS = [(0.0, MU_MAX), KTT_BOUNDS] + [(-THETA_BOUND, THETA_BOUND)] * (NPAR - 2)


def fit(model, n, x0, fixed=None):
    from scipy.optimize import minimize
    fixed = fixed or {}
    base = np.array(x0, dtype=float)
    for p, v in fixed.items():
        base[IDX[p]] = v
    free = [i for i in range(NPAR) if PARS[i] not in fixed]

    def f(z):
        x = base.copy()
        x[free] = z
        v, g = model.nll(x, n)
        return v, g[free]

    res = minimize(f, base[free], jac=True, method="L-BFGS-B", bounds=[BOUNDS[i] for i in free],
                   options={"maxiter": 5000, "maxfun": 20000, "ftol": 1e-14, "gtol": 1e-9})
    x = base.copy()
    x[free] = res.x
    return {"nll": float(res.fun), "x": x, "ok": bool(res.success), "message": str(res.message)}


def start_point(mu=0.0, k=1.0):
    x = np.zeros(NPAR)
    x[0], x[1] = mu, k
    return x


def fit_multi(model, n, starts, fixed=None):
    fits = [fit(model, n, s, fixed) for s in starts]
    best = min(fits, key=lambda r: r["nll"])
    return best, [r["nll"] for r in fits]


def analyse(model, n, scan=True, fixed=None):
    """Signal+background and background-only fits, q0, Z, G4 numbers and (optionally) the profile scan."""
    fixed = dict(fixed or {})
    sb, sb_nlls = fit_multi(model, n, [start_point(m) for m in MU_STARTS], fixed)
    fb = dict(fixed)
    fb["mu"] = 0.0
    b, b_nlls = fit_multi(model, n, [start_point(0.0, k) for k in KTT_STARTS_B], fb)
    mu_hat = float(sb["x"][0])
    q0 = max(0.0, 2.0 * (b["nll"] - sb["nll"])) if mu_hat > 0 else 0.0
    out = {"target": model.target, "mu_hat": mu_hat, "k_tt": float(sb["x"][1]), "q0": q0, "Z": math.sqrt(q0),
           "nll_sb": sb["nll"], "nll_b": b["nll"], "nll_sb_starts": sb_nlls, "nll_b_starts": b_nlls,
           "G4": gate_g4(sb_nlls, b_nlls), "params_sb": dict(zip(PARS, map(float, sb["x"]))),
           "params_b": dict(zip(PARS, map(float, b["x"]))), "fixed": fixed}
    nu_hat = model.nu_jac(sb["x"])[0]
    out["pulls"] = ((n - nu_hat) / np.sqrt(nu_hat)).tolist()
    out["nu_hat"] = nu_hat.tolist()
    if scan:
        grid = scan_grid()
        prev = sb["x"].copy()
        vals = []
        for m in grid:
            fm = dict(fixed)
            fm["mu"] = m
            r1 = fit(model, n, prev, fm)
            r2 = fit(model, n, sb["x"], fm)
            r = r1 if r1["nll"] <= r2["nll"] else r2
            prev = r["x"]
            vals.append(r["nll"])
        nll_min = min(sb["nll"], min(vals))
        qs = [2.0 * (v - nll_min) for v in vals]
        lo, hi = interval(grid, qs)
        out.update({"scan_mu": grid, "scan_q": qs, "lo68": lo, "hi68": hi,
                    "interval_at_edge": hi is not None and hi >= SCAN_MAX - 1e-9})
    return out


# ======================================================================= synthetic skims (selftest)
def synthetic_skims(seed=7, mu_data=MU_CMS):
    """Skim-like JSON dicts with invented but plausible templates, for testing the machinery only."""
    rng = np.random.default_rng(seed)
    c = _arr(class_centres())
    xc = _arr([-0.75, -0.25, 0.25, 0.75])

    def prof(mp, slope):
        mp = _arr(mp) / float(sum(mp))
        cp = 1.0 + slope * xc
        return np.outer(mp, cp / cp.sum()).ravel()

    mb = np.repeat(np.arange(NM), NCOS)
    tilt = (mb - 4) / 4.0

    def variations(H, vars_):
        out = {"nom": H}
        shifts = {"mur": 0.03, "muf": 0.02, "isr": 0.04, "fsr": -0.03, "jes": 0.05}
        for name, a in shifts.items():
            if name + "_up" in vars_:
                rate = 0.02 if name == "jes" else 0.05
                out[name + "_up"] = H * (1.0 + a * tilt + rate)
                out[name + "_dn"] = H * (1.0 - a * tilt - rate * 0.8)
        return out

    def rec(xsec, hv, gen, sumw, n_ent):
        fd = [{"uri": "synthetic", "url": "synthetic", "entries": n_ent, "sumw_events": sumw,
               "sumw_runs": sumw, "seconds": 0.0}]
        return {"title": "synthetic", "xsec_pb": xsec, "files_assigned": ["synthetic"], "files_done": fd,
                "files_failed": [], "files_skipped_time": [], "h": {v: a.tolist() for v, a in hv.items()},
                "h2_nom": hv["nom"].tolist(), "gen_sumw": list(map(float, gen)), "n_entries": n_ent,
                "sumw_events": sumw, "sumw2_events": sumw, "sumw_runs": sumw,
                "gen_source": {"hard": n_ent, "fallback": 0, "none": 0}}

    def skim(sample, kind, recs):
        return {"version": SKIM_VERSION, "sample": sample, "kind": kind, "part": 0, "nparts": 1,
                "max_files": 0, "max_entries_per_file": 0, "constants": skim_constants(), "records": recs}

    N = 9.0e6
    G = np.zeros(NCLS_TT)
    G[:NWIN] = 0.06 * N * (1.0 - 0.3 * c) / 2.0 / 20.0
    G[40], G[41], G[42] = 0.12 * N, 0.80 * N, 0.02 * N
    eff = 0.05
    H = np.zeros((NCLS_TT, NBIN))
    win_m = [0.02, 0.25, 0.35, 0.22, 0.10, 0.04, 0.015, 0.004, 0.001]
    for k in range(NWIN):
        H[k] = G[k] * eff * prof(win_m, 0.6 * c[k])
    H[40] = G[40] * eff * prof([0.01, 0.12, 0.25, 0.25, 0.18, 0.12, 0.05, 0.015, 0.005], 0.1)
    H[41] = G[41] * eff * prof([0.005, 0.03, 0.08, 0.13, 0.18, 0.22, 0.18, 0.12, 0.055], -0.1)
    H[42] = G[42] * eff * prof([1.0] * NM, 0.0)
    H = H * (1.0 + 0.02 * rng.standard_normal(H.shape))
    tt = skim("tt", "tt", {REC_TT: rec(XS_TT2L, variations(H, TT_VARS), G, N, int(N))})
    tw_recs = {}
    for rid in (REC_TW, REC_TWB):
        Nw = 3.4e6
        Hw = (Nw * 0.012 * prof([0.01, 0.08, 0.15, 0.18, 0.18, 0.16, 0.12, 0.08, 0.04], 0.0))[None, :]
        tw_recs[rid] = rec(XS_TW, variations(Hw, BKG_VARS), [Nw], Nw, int(Nw))
    tw = skim("tw", "bkg", tw_recs)
    Nd = 1.0e7
    Hd = (Nd * 2e-5 * prof([0.3, 0.25, 0.15, 0.1, 0.08, 0.05, 0.04, 0.02, 0.01], 0.2))[None, :]
    dy = skim("dy", "bkg", {REC_DY: rec(XS_DY, variations(Hd, BKG_VARS), [Nd], Nd, int(Nd))})
    mc = [tt, tw, dy]
    tpl = build_templates(merge_skims(mc))
    n = np.round(Model(tpl, "singlet").asimov(mu_data))
    nG = np.floor(n / 2.0)
    data = []
    for samp, rid, arr in (("data_G", REC_G, nG), ("data_H", REC_H, n - nG)):
        d = rec(None, {"nom": arr[None, :]}, [0.0], 0.0, int(arr.sum()))
        d["files_done"][0].update({"sumw_events": None, "sumw_runs": None})
        data.append(skim(samp, "data", {rid: d}))
    return mc, data


# ======================================================================= modes
def log(fh, line):
    print(line, flush=True)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


def summarise(a):
    keys = ["target", "mu_hat", "k_tt", "q0", "Z", "nll_sb", "nll_b", "G4", "lo68", "hi68", "interval_at_edge"]
    return {k: a.get(k) for k in keys}


def run_asimov(tpl, fh, scan=True):
    res = {}
    for tgt in TARGETS:
        model = Model(tpl, tgt)
        res[tgt] = {}
        for mu_true in (MU_CMS, MU_NRQCD):
            nA = model.asimov(mu_true)
            a = analyse(model, nA, scan=scan and tgt == "singlet" and mu_true == MU_CMS)
            res[tgt]["%.1f" % mu_true] = a
            log(fh, "CMSFIT_ASIMOV target=%s mu_true=%.1f mu_hat=%.3f k_tt=%.4f Z_A=%.3f G4=%s lo68=%s hi68=%s" % (
                tgt, mu_true, a["mu_hat"], a["k_tt"], a["Z"], a["G4"], a.get("lo68"), a.get("hi68")))
    return res


def yields(tpl):
    return {"tt_low": float(tpl["low"]["nom"].sum()), "tt_rest": float(tpl["rest"]["nom"].sum()),
            "tW": float(tpl["TW"]["nom"].sum()), "other": float(tpl["O"]["nom"].sum()),
            "S_per_pb": {t: float(tpl["S"][t]["nom"].sum()) for t in TARGETS},
            "D_win": tpl["D_win"], "max_delta": float(tpl["delta"].max())}


def run_data(tpl, fh, scan=True):
    n = tpl["data"]
    out = {"n_total": float(n.sum()), "n": n.tolist(), "fits": {}}
    for tgt in TARGETS:
        a = analyse(Model(tpl, tgt), n, scan=scan and tgt == "singlet")
        out["fits"][tgt] = a
        log(fh, "CMSFIT_DATA target=%s mu_hat=%.3f k_tt=%.4f Z=%.3f nll_sb=%.4f nll_b=%.4f G4=%s lo68=%s hi68=%s" % (
            tgt, a["mu_hat"], a["k_tt"], a["Z"], a["nll_sb"], a["nll_b"], a["G4"], a.get("lo68"), a.get("hi68")))
    info = {}
    info["no_thr"] = summarise(analyse(Model(tpl, "singlet"), n, scan=False, fixed={"th_thr": 0.0}))
    stat_fixed = {p: 0.0 for p in PARS[2:]}
    info["stat_only"] = summarise(analyse(Model(tpl, "singlet", bb=False), n, scan=False, fixed=stat_fixed))
    out["info"] = info
    log(fh, "CMSFIT_INFO no_thr mu_hat=%.3f Z=%.3f | stat_only mu_hat=%.3f Z=%.3f" % (
        info["no_thr"]["mu_hat"], info["no_thr"]["Z"], info["stat_only"]["mu_hat"], info["stat_only"]["Z"]))
    return out


def load_dir(d):
    skims = []
    for p in sorted(glob.glob(os.path.join(d, "**", "skim-*.json"), recursive=True)):
        with open(p) as f:
            skims.append(json.load(f))
    return skims


def fetch(manifest, dest):
    """Download committed skim files, each pinned by its git blob sha1; abort on any mismatch."""
    with open(manifest) as f:
        items = json.load(f)
    os.makedirs(dest, exist_ok=True)
    for it in items:
        url = REPO_RAW + it["commit"] + "/" + it["path"]
        req = urllib.request.Request(url, headers={"User-Agent": "cmsfit/0.1 (open-data research)"})
        with urllib.request.urlopen(req, timeout=300) as r:
            data = r.read()
        sha = git_blob_sha1(data)
        if sha != it["git_sha1"]:
            raise SystemExit("PIN MISMATCH %s: %s != %s" % (it["path"], sha, it["git_sha1"]))
        with open(os.path.join(dest, os.path.basename(it["path"])), "wb") as f:
            f.write(data)
        print("CMSFIT_FETCH OK %s %s" % (it["path"], sha), flush=True)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if np is not None and isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if np is not None and isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return str(o)
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["selftest", "asimov", "data"])
    ap.add_argument("--skim-dir", default="skims")
    ap.add_argument("--fetch", default="")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--no-scan", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()
    fh = open(os.path.join(args.out_dir, "cmsfit-%s.txt" % args.mode), "w")
    import scipy
    log(fh, "CMSFIT_START version=%s mode=%s python=%s numpy=%s scipy=%s" % (
        VERSION, args.mode, sys.version.split()[0], np.__version__, scipy.__version__))
    scan = not args.no_scan
    if args.mode == "selftest":
        mc, data = synthetic_skims()
        skims = mc + data
    else:
        if args.fetch:
            fetch(args.fetch, args.skim_dir)
        skims = load_dir(args.skim_dir)
    kinds = sorted({s["kind"] for s in skims})
    log(fh, "CMSFIT_SKIMS n=%d kinds=%s samples=%s" % (len(skims), kinds, sorted({s["sample"] for s in skims})))
    if args.mode == "asimov" and "data" in kinds:
        log(fh, "CMSFIT_BLIND_ABORT data skims given in asimov mode")
        return 7
    merged = merge_skims(skims)
    g1, g1_why = gate_g1(merged)
    log(fh, "CMSFIT_GATE G1 %s %s" % ("PASS" if g1 else "FAIL", "; ".join(g1_why)))
    tpl = build_templates(merged)
    y = yields(tpl)
    log(fh, "CMSFIT_YIELDS " + json.dumps(y))
    result = {"version": VERSION, "mode": args.mode, "gates": {"G1": {"pass": g1, "why": g1_why}},
              "yields": y, "delta": tpl["delta"].tolist(),
              "templates_nom": {"tt_low": tpl["low"]["nom"].tolist(), "tt_rest": tpl["rest"]["nom"].tolist(),
                                "tW": tpl["TW"]["nom"].tolist(), "other": tpl["O"]["nom"].tolist(),
                                "S_singlet_per_pb": tpl["S"]["singlet"]["nom"].tolist()}}
    asi = run_asimov(tpl, fh, scan=scan)
    result["asimov"] = asi
    a88 = asi["singlet"]["%.1f" % MU_CMS]
    g4a = a88["G4"]
    h4 = grade_h4(a88["Z"]) if (g1 and g4a) else "NOT GRADABLE"
    result["H4"] = {"Z_A": a88["Z"], "Z_A_6.4": asi["singlet"]["%.1f" % MU_NRQCD]["Z"], "G4_asimov": g4a,
                    "grade": h4}
    log(fh, "CMSFIT_H4 Z_A=%.3f Z_A_6.4=%.3f G4=%s grade=%s" % (
        a88["Z"], asi["singlet"]["%.1f" % MU_NRQCD]["Z"], g4a, h4))
    if args.mode in ("data", "selftest") and tpl["data"] is not None:
        g2, g2_why = gate_g2(merged)
        log(fh, "CMSFIT_GATE G2 %s %s" % ("PASS" if g2 else "FAIL", "; ".join(g2_why)))
        d = run_data(tpl, fh, scan=scan)
        s = d["fits"]["singlet"]
        g3, g4 = gate_g3(s["k_tt"]), s["G4"]
        log(fh, "CMSFIT_GATE G3 %s k_tt=%.4f" % ("PASS" if g3 else "FAIL", s["k_tt"]))
        log(fh, "CMSFIT_GATE G4 %s sb_spread=%.5f b_spread=%.5f" % (
            "PASS" if g4 else "FAIL", spread(s["nll_sb_starts"]), spread(s["nll_b_starts"])))
        ok = g1 and g2 and g3 and g4
        h1 = grade_h1(s["Z"]) if ok else "NOT GRADABLE"
        h2 = grade_h2(s["mu_hat"], s.get("lo68"), s.get("hi68")) if ok else "NOT GRADABLE"
        h3 = (grade_h3(s["nll_sb"], d["fits"]["unchanged"]["nll_sb"], d["fits"]["triplet"]["nll_sb"])
              if ok else "NOT GRADABLE")
        d["gates"] = {"G2": {"pass": g2, "why": g2_why}, "G3": g3, "G4": g4}
        d["grades"] = {"H1": h1, "H2": h2, "H3": h3}
        result["data"] = d
        log(fh, "CMSFIT_GRADES H1=%s H2=%s H3=%s (Z=%.3f mu_hat=%.3f [%s, %s])" % (
            h1, h2, h3, s["Z"], s["mu_hat"], s.get("lo68"), s.get("hi68")))
        if args.mode == "selftest":
            log(fh, "CMSFIT_SELFTEST note=synthetic templates, never graded; injected mu=%.1f" % MU_CMS)
    result["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(args.out_dir, "cmsfit-%s.json" % args.mode), "w") as f:
        json.dump(jsonable(result), f)
    log(fh, "CMSFIT_DONE seconds=%.1f" % result["seconds"])
    fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
