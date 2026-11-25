# cmsskim-tt-p1-v1: returned 2026-11-25 (information for gate G1; nothing graded)

Job: `jobs/cmsskim-tt-p1-v1.json` (cmsskim 0.1.0, frozen with PREREG-CMSTOPO-v0.1). Ran on 2026-11-25
19:50-20:37 UTC on 2 CPUs, peak memory 0.27 GB. The files here are copied unchanged from the job's `/out`.
`sha256.txt` shows the same cmsskim.py, test_cmsskim.py and PREREG bytes as the DY and tt-p3 jobs
(`01e87d0d...`, `60046f7f...`, `6bc600d7...`).

- Unit tests: 17 of 17 passed (`tests.txt`).
- Record 67801 (TTTo2L2Nu powheg-pythia8, 2016 post-VFP NanoAODv9), part 1 of 4 of the sorted file list:
  12 files assigned, 12 processed, 0 failed, 0 skipped for time. 11,210,000 events were read and
  566,491 were selected (nominal e-mu selection).
- Generator-level check: has_gen = 1.0000 in every file (20,000 hard-process events checked per file,
  0 fallback, 0 missing). The pre-registered per-part floor is 0.95.
- Per-file weight check, |sum genWeight - genEventSumw| / genEventSumw: about 4.4e-8 in every file and
  4.4e-8 for the part total (808,416,472.8 against 808,416,436.9). The pre-registered limit is 1e-6.

For part 1 the top-pair parts of gate G1 pass. With parts 1 and 3 in, 22.87M top-pair events have been
read, already above the G1 floor of 8M; G1 as a whole is still graded by `cmsfit.py --mode asimov`
once all six simulation skims (tt p0-p3, tW, DY) are in. Still queued: tt p0, tt p2, tW.

Made with AI agents working in MV2's name. Apache-2.0 (code), CC BY 4.0 (text).
