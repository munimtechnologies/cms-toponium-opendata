# cmsskim-tt-p0-v1: returned 2026-11-26 (information for gate G1; nothing graded)

Job: `jobs/cmsskim-tt-p0-v1.json` (cmsskim 0.1.0, frozen with PREREG-CMSTOPO-v0.1). Ran on 2026-11-26
11:17-11:53 UTC on 2 CPUs, peak memory 0.28 GB. The files here are copied unchanged from the job's `/out`.
`sha256.txt` shows the same cmsskim.py, test_cmsskim.py and PREREG bytes as the DY, tt-p3 and tt-p1 jobs
(`01e87d0d...`, `60046f7f...`, `6bc600d7...`).

- Unit tests: 17 of 17 passed (`tests.txt`).
- Record 67801 (TTTo2L2Nu powheg-pythia8, 2016 post-VFP NanoAODv9), part 0 of 4 of the sorted file list:
  13 files assigned, 13 processed, 0 failed, 0 skipped for time. 9,392,000 events were read and
  473,793 were selected (nominal e-mu selection).
- Generator-level check: has_gen = 1.0000 in every file (20,000 hard-process events checked per file,
  0 fallback, 0 missing). The pre-registered per-part floor is 0.95.
- Per-file weight check, |sum genWeight - genEventSumw| / genEventSumw: about 4.4e-8 in every file and
  4.4e-8 for the part total (677,257,102.0 against 677,257,072.0). The pre-registered limit is 1e-6.

For part 0 the top-pair parts of gate G1 pass. With parts 0, 1 and 3 in, about 32.26M top-pair events
have been read (G1 floor: 8M). G1 as a whole is still graded by `cmsfit.py --mode asimov` once all six
simulation skims (tt p0-p3, tW, DY) are in. Still queued: tt p2, tW, the cmsfit self-test and cmsprobe.

(Committed in a later run than the one that wrote it: the first commit request of 2026-11-26 was not
processed, so it was filed again unchanged.)

Made with AI agents working in MV2's name. Apache-2.0 (code), CC BY 4.0 (text).
