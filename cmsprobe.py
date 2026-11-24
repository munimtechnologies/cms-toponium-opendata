# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""cmsprobe 0.1.0 -- feasibility probe for a CMS 2016 open-data toponium test.

It answers only plumbing questions, and looks at NO physics distribution:
  1. which CERN Open Data records hold the 2016 (UL, NanoAODv9) MuonEG data, the dilepton
     top-pair simulation and the main backgrounds, and how big they are;
  2. whether NanoAOD files can be read remotely with uproot (which URL form works);
  3. whether every branch the planned analysis needs exists in data and in simulation;
  4. how fast the needed branches stream, so the full skim can be sized and split.
The only event counts it prints are raw entry counts and the fraction of entries with at least
one electron and one muon above 15 GeV (to size the skim output), never a mass or angle.

Usage: python cmsprobe.py --out-dir /out [--entries 100000]
"""
import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

VERSION = "0.1.0"
API = "https://opendata.cern.ch/api/records/"

# Known records (found 2026-11-24 through the portal API).
RECORDS = {
    "data_muoneg_2016G": 30528,
    "data_muoneg_2016H": 30561,
    "mc_ttbar_2l2nu_ul16": 67801,
}

# Free-text searches for the backgrounds and the luminosity / validated-runs records.
SEARCHES = {
    "tW_top": "ST_tW_top_5f_NoFullyHadronicDecays RunIISummer20UL16NanoAODv9",
    "tW_antitop": "ST_tW_antitop_5f_NoFullyHadronicDecays RunIISummer20UL16NanoAODv9",
    "DY_M50": "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8 RunIISummer20UL16NanoAODv9",
    "DY_M10to50": "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8 RunIISummer20UL16NanoAODv9",
    "WW_2l2nu": "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8 RunIISummer20UL16NanoAODv9",
    "ttbar_hdamp_up": "TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8 RunIISummer20UL16NanoAODv9",
    "ttbar_mtop1715": "TTTo2L2Nu_mtop171p5_TuneCP5_13TeV-powheg-pythia8 RunIISummer20UL16NanoAODv9",
    "ttbar_mtop1735": "TTTo2L2Nu_mtop173p5_TuneCP5_13TeV-powheg-pythia8 RunIISummer20UL16NanoAODv9",
    "eta_t_or_A365": "HpseudoToTTTo2L2Nu_m365 RunIISummer20UL16NanoAODv9",
    "validated_runs_2016": "CMS 2016 validated runs list Collisions16 JSON",
    "luminosity_2016": "CMS luminosity information 2016",
}

COMMON_BRANCHES = [
    "run", "luminosityBlock", "event", "PV_npvsGood",
    "nElectron", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass",
    "Electron_charge", "Electron_cutBased", "Electron_deltaEtaSC", "Electron_pfRelIso03_all",
    "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge", "Muon_tightId",
    "Muon_pfRelIso04_all",
    "nJet", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_jetId", "Jet_puId",
    "Jet_btagDeepFlavB",
    "MET_pt", "MET_phi", "PuppiMET_pt", "PuppiMET_phi",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
    "HLT_IsoMu24", "HLT_IsoTkMu24", "HLT_Ele27_WPTight_Gsf",
    "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter", "Flag_eeBadScFilter",
]
MC_BRANCHES = [
    "genWeight", "Pileup_nTrueInt", "LHEScaleWeight", "PSWeight", "LHEPdfWeight",
    "nGenPart", "GenPart_pdgId", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass",
    "GenPart_statusFlags", "GenPart_genPartIdxMother", "GenPart_status",
    "nLHEPart", "LHEPart_pdgId", "LHEPart_pt", "LHEPart_eta", "LHEPart_phi", "LHEPart_mass",
    "Jet_hadronFlavour", "Jet_genJetIdx", "Electron_genPartIdx", "Muon_genPartIdx",
]


# ---------------------------------------------------------------- pure helpers (unit-tested)
def extract_root_uris(obj, out=None):
    """All strings in a JSON object that look like CMS open-data ROOT file URIs, in order, unique."""
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
    """Candidate remote URLs for one root:// URI: HTTPS through the portal, HTTPS straight to EOS,
    then the original XRootD URI."""
    prefix = "root://eospublic.cern.ch/"
    if not uri.startswith(prefix):
        raise ValueError("not an eospublic URI: %r" % uri)
    path = uri[len(prefix):]          # '/eos/opendata/...'
    path = "/" + path.lstrip("/")
    return [
        "https://opendata.cern.ch" + path,
        "https://eospublic.cern.ch:443/" + path,
        uri,
    ]


def extrapolate_hours(seconds, entries_read, total_entries, workers=1):
    """Hours to stream total_entries at the measured rate, split over workers."""
    if entries_read <= 0 or seconds <= 0:
        return None
    return total_entries * (seconds / entries_read) / 3600.0 / max(1, workers)


def missing(wanted, present):
    present = set(present)
    return [b for b in wanted if b not in present]


def record_summary(meta):
    """Compact summary of a portal record's metadata dict."""
    dist = meta.get("distribution") or {}
    return {
        "recid": meta.get("recid"),
        "title": meta.get("title"),
        "number_events": dist.get("number_events"),
        "number_files": dist.get("number_files"),
        "size": dist.get("size"),
        "doi": meta.get("doi"),
    }


# ---------------------------------------------------------------- network
def get_json(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cmsprobe/0.1 (open-data research)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(5 * (i + 1))
    raise RuntimeError("GET failed %s: %s" % (url, last))


def log(fh, line):
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def probe_file(uri, wanted, n_entries, fh):
    import numpy as np
    import uproot
    import awkward as ak
    res = {"uri": uri, "opened_with": None, "errors": []}
    for url in url_forms(uri):
        handlers = [None]
        if url.startswith("https://"):
            handlers = [uproot.source.http.HTTPSource, None]
        for h in handlers:
            try:
                t0 = time.time()
                kw = {"timeout": 120}
                if h is not None:
                    kw["handler"] = h
                f = uproot.open(url, **kw)
                tree = f["Events"]
                res["opened_with"] = {"url": url, "handler": getattr(h, "__name__", "default"),
                                      "open_s": round(time.time() - t0, 2)}
                break
            except Exception as exc:  # noqa: BLE001
                res["errors"].append("%s [%s]: %s" % (url, getattr(h, "__name__", "default"),
                                                       repr(exc)[:300]))
        if res["opened_with"]:
            break
    if not res["opened_with"]:
        log(fh, "CMSPROBE_FILE FAIL %s" % uri)
        return res
    names = list(tree.keys())
    res["n_entries"] = int(tree.num_entries)
    res["n_branches"] = len(names)
    res["missing"] = missing(wanted, names)
    use = [b for b in wanted if b in names]
    n = min(n_entries, res["n_entries"])
    src = f.file.source
    b0 = getattr(src, "num_requested_bytes", None)
    t0 = time.time()
    n_emu = 0
    done = 0
    for arr in tree.iterate(use, entry_stop=n, step_size=25000, library="ak"):
        ne = ak.sum(arr["Electron_pt"] > 15, axis=1)
        nm = ak.sum(arr["Muon_pt"] > 15, axis=1)
        n_emu += int(np.sum(np.asarray((ne >= 1) & (nm >= 1))))
        done += len(arr)
    dt = time.time() - t0
    b1 = getattr(src, "num_requested_bytes", None)
    res.update({"entries_read": done, "read_s": round(dt, 2),
                "entries_per_s": round(done / dt, 1) if dt > 0 else None,
                "bytes_read": (b1 - b0) if (b0 is not None and b1 is not None) else None,
                "frac_e15_mu15": round(n_emu / done, 5) if done else None})
    log(fh, "CMSPROBE_FILE OK %s entries=%d branches=%d missing=%d read=%d in %.1fs (%.0f/s) bytes=%s "
        "frac_e15mu15=%s" % (res["opened_with"]["url"], res["n_entries"], res["n_branches"],
                             len(res["missing"]), done, dt, res["entries_per_s"] or 0,
                             res["bytes_read"], res["frac_e15_mu15"]))
    if res["missing"]:
        log(fh, "CMSPROBE_MISSING %s" % ",".join(res["missing"]))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--entries", type=int, default=100000)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    out = {"version": VERSION, "records": {}, "searches": {}, "files": {}}
    fh = open(os.path.join(args.out_dir, "probe.txt"), "w")
    log(fh, "CMSPROBE_VERSION %s python %s" % (VERSION, sys.version.split()[0]))
    try:
        import uproot, awkward, numpy  # noqa: E401
        out["versions"] = {"uproot": uproot.__version__, "awkward": awkward.__version__,
                           "numpy": numpy.__version__}
        log(fh, "CMSPROBE_VERSIONS %s" % json.dumps(out["versions"]))
    except Exception as exc:  # noqa: BLE001
        log(fh, "CMSPROBE_VERSIONS FAIL %r" % exc)

    uris = {}
    for key, recid in RECORDS.items():
        try:
            rec = get_json(API + str(recid))
            meta = rec.get("metadata", rec)
            s = record_summary(meta)
            u = extract_root_uris(rec)
            s["n_uris_found"] = len(u)
            uris[key] = u
            out["records"][key] = s
            log(fh, "CMSPROBE_RECORD %s %s" % (key, json.dumps(s)))
        except Exception as exc:  # noqa: BLE001
            out["records"][key] = {"error": repr(exc)[:300]}
            log(fh, "CMSPROBE_RECORD %s FAIL %r" % (key, exc))

    for key, q in SEARCHES.items():
        try:
            res = get_json(API + "?" + urllib.parse.urlencode({"q": q, "size": 8}))
            hits = res.get("hits", {}).get("hits", [])
            rows = [record_summary(h.get("metadata", {})) for h in hits]
            out["searches"][key] = {"q": q, "total": res.get("hits", {}).get("total"), "hits": rows}
            log(fh, "CMSPROBE_SEARCH %s total=%s" % (key, out["searches"][key]["total"]))
            for r in rows:
                log(fh, "CMSPROBE_HIT %s %s %s events=%s files=%s size=%s" % (
                    key, r["recid"], r["title"], r["number_events"], r["number_files"], r["size"]))
        except Exception as exc:  # noqa: BLE001
            out["searches"][key] = {"q": q, "error": repr(exc)[:300]}
            log(fh, "CMSPROBE_SEARCH %s FAIL %r" % (key, exc))

    plan = [("data_muoneg_2016G", COMMON_BRANCHES), ("mc_ttbar_2l2nu_ul16", COMMON_BRANCHES + MC_BRANCHES)]
    for key, wanted in plan:
        u = uris.get(key) or []
        if not u:
            log(fh, "CMSPROBE_FILE SKIP %s no URIs" % key)
            continue
        try:
            out["files"][key] = probe_file(u[0], wanted, args.entries, fh)
        except Exception as exc:  # noqa: BLE001
            out["files"][key] = {"uri": u[0], "error": repr(exc)[:500]}
            log(fh, "CMSPROBE_FILE ERROR %s %r" % (key, exc))

    # Extrapolate the full skim: every file of the three records, at each kind's measured rate.
    est = {}
    for key, kind in (("data_muoneg_2016G", "data_muoneg_2016G"), ("data_muoneg_2016H", "data_muoneg_2016G"),
                      ("mc_ttbar_2l2nu_ul16", "mc_ttbar_2l2nu_ul16")):
        f = out["files"].get(kind, {})
        tot = (out["records"].get(key) or {}).get("number_events")
        if f.get("read_s") and f.get("entries_read") and tot:
            est[key] = {"total_entries": tot,
                        "hours_1_worker": round(extrapolate_hours(f["read_s"], f["entries_read"], tot), 2),
                        "hours_2_workers": round(extrapolate_hours(f["read_s"], f["entries_read"], tot, 2), 2)}
            if f.get("bytes_read"):
                est[key]["est_gb_transferred"] = round(f["bytes_read"] / f["entries_read"] * tot / 1e9, 1)
    out["estimate"] = est
    log(fh, "CMSPROBE_ESTIMATE %s" % json.dumps(est, sort_keys=True))
    ok = all(out["files"].get(k, {}).get("opened_with") for k, _ in plan)
    log(fh, "CMSPROBE_DONE %s" % ("READABLE" if ok else "NOT_READABLE"))
    with open(os.path.join(args.out_dir, "probe.json"), "w") as jf:
        json.dump(out, jf, indent=1, sort_keys=True)
    fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
