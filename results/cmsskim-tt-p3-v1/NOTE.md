# cmsskim-tt-p3-v1: returned 2026-11-24 (information for gate G1; nothing graded)

Job: `jobs/cmsskim-tt-p3-v1.json` (cmsskim 0.1.0, frozen with PREREG-CMSTOPO-v0.1). Ran on 2026-11-24
21:28-22:18 UTC on 2 CPUs, peak memory 0.28 GB. The files here are copied unchanged from the job's `/out`.
`sha256.txt` shows the same cmsskim.py, test_cmsskim.py and PREREG bytes as the DY job
(`results/cmsskim-dy-v1/sha256.txt`).

- Unit tests: 17 of 17 passed (`tests.txt`).
- Record 67801 (TTTo2L2Nu powheg-pythia8, 2016 post-VFP NanoAODv9), part 3 of 4 of the sorted file list:
  12 files assigned, 12 processed, 0 failed, 0 skipped for time. 11,657,000 events were read and
  588,388 were selected (nominal e-mu selection).
- Generator-level check: has_gen = 1.0000 in every file (20,000 hard-process events checked per file,
  with 0 fallback and 0 missing). The pre-registered per-part floor is 0.95.
- Per-file weight check, |sum genWeight - genEventSumw| / genEventSumw: about 4.4e-8 in every file and
  4.4e-8 for the part total (840,566,861.2 against 840,566,823.9). The pre-registered limit is 1e-6.

For part 3, the top-pair parts of gate G1 pass: has_gen is at least 0.95 and the weight sums agree.
The G1 floor of at least 8M top-pair events is checked over all four parts once p0-p2 return.
G1 as a whole is graded by `cmsfit.py --mode asimov` once all six simulation skims are in.

Made with AI agents working in MV2's name. Apache-2.0 (code), CC BY 4.0 (text).
