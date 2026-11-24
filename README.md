# cms-toponium-opendata

An open, pre-registered attempt by MV2 of Munim, Inc. to see the top-quark-pair threshold
excess ("toponium") in the **public** CMS 2016 collision data (runs G and H, NanoAOD, CERN Open
Data portal), using only public files and rules fixed before any distribution is looked at.

Made with AI agents working in MV2's name. No human expert has checked it yet. This project is
not affiliated with the CMS Collaboration. Code: Apache-2.0. Text: CC BY 4.0.

## Status (2026-11-24)

- `PREREG-CMSTOPO-v0.1.md` (**frozen**): selection, the m(llbb) x cos(phi*) histogram, the
  toponium-like signal template built by spin-reweighting the near-threshold top-pair simulation,
  the likelihood fit, four gates, four graded hypotheses and our predictions. The wow, if it works:
  the first detection of the top-pair threshold excess outside the LHC collaborations, from public
  files.
- `cmsskim.py` (0.1.0) + `test_cmsskim.py` (17 tests): streams the public NanoAOD files and fills
  the pre-registered histograms, with scale, parton-shower and jet-energy variations and the
  generator spin class of every simulated top pair. Six simulation jobs are queued
  (`jobs/cmsskim-*.json`).
- `cmsfit.py` (0.1.0) + `test_cmsfit.py` (22 tests: 12 stdlib, 10 numeric): the pre-registered binned likelihood fit
  (section 6): top-pair, toponium-like signal, tW and other templates; piecewise-linear
  nuisances; Barlow-Beeston-lite simulation statistics; five-start L-BFGS-B fits; q0 and Z;
  profile scan for the 68 % interval; Asimov expected sensitivity; gates G1-G4 and grades H1-H4.
  Implementation details that v0.1 left open are fixed in `AMENDMENT-CMSTOPO-0.1.1.md`, written
  before any histogram returned. A self-test on invented templates is queued
  (`jobs/cmsfit-selftest-v1.json`).
- Order of operations: simulation skims -> fit on simulation only (expected sensitivity, H4
  graded) -> only then the data skims and the data fit (H1-H3 graded, whatever they say).
- `SCOUT-CMSTOPO-v0.1.md`: the scouting note. `cmsprobe.py` (0.1.0): plumbing probe
  (`jobs/cmsprobe-v1.json`, queued).

## Public records used

MuonEG 2016G/H NanoAOD (30528, 30561); TTTo2L2Nu (67801); tW (64895, 64839); DY M>50 (35669);
validated runs (14220); luminosity (1059): 16.39 fb-1.

## Related

- MV2's ATLAS open-data study of the 2019 top spin-correlation excess:
  https://github.com/munimtechnologies/ttbar-threshold-opendata
- CMS, Observation of a pseudoscalar excess at the top quark pair production threshold,
  https://arxiv.org/abs/2503.22382

## License

Apache-2.0 (see LICENSE).
