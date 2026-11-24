# cmsskim-dy-v1: returned 2026-11-24 (information for gate G1; nothing graded)

Job: `jobs/cmsskim-dy-v1.json` (cmsskim 0.1.0, frozen with PREREG-CMSTOPO-v0.1). Ran on 2026-11-24
19:56-20:19 UTC, 2 CPUs, peak memory 0.16 GB. Files here are copied unchanged from the job's `/out`.

- Unit tests: 17 of 17 passed (`tests.txt`).
- Record 35669 (DYJetsToLL_M-50 amcatnloFXFX, 2016 post-VFP NanoAODv9): the first 6 of the 41 sorted
  files, as pre-registered. 6 of 6 processed, 0 failed, 0 skipped for time (one transient connection
  retry, then read). 9,408,214 events read, 40 selected (nominal), 42 (JES up), 40 (JES down).
- Per-file check |sum genWeight - genEventSumw| / genEventSumw: below 1e-8 in every file (largest
  about 7e-9; the pre-registered limit is 1e-6).

So the two parts of gate G1 that concern DY ("at least one file of DY" and the per-file weight check)
pass. G1 as a whole is graded by `cmsfit.py --mode asimov` once all six simulation skims are in.

Desk arithmetic, information only (the fit computes the official numbers): the 40 selected events carry a
net of 21 positive weights (aMC@NLO has negative weights), so the DY "other" template is about
21 x 25,259.67 x 16,393.4 pb^-1 x 6,077.22 pb / 1.5983e11 = about 330 expected events in 16.4 fb^-1,
with a simulation-statistics uncertainty near 30 %. That is small next to the top-pair yield, and the
pre-registered 50 % prior on "other" plus the Barlow-Beeston-lite terms cover it. Two of the 36 bins have
a small negative DY sum; the fit floors only the total expected count, which stays positive.

Made with AI agents working in MV2's name. Apache-2.0 (code), CC BY 4.0 (text).
