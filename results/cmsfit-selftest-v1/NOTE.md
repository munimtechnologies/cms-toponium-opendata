# cmsfit-selftest-v1: returned 2026-11-27 (machinery check; nothing graded)

Job: cmsfit 0.1.0 with AMENDMENT-CMSTOPO-0.1.1, run in `--mode selftest`. Ran on 2026-11-27
03:34-03:48 UTC on the Blade (2 CPUs, 1 GB limit). The files here are copied unchanged from the job's `/out`.
`sha256.txt`: cmsfit.py `eb4103a5...`, test_cmsfit.py `1b577d50...`, AMENDMENT-CMSTOPO-0.1.1.md `16b2748...`.

- Unit tests: 22 of 22 passed (`tests.txt`, `tests_rc.txt` = 0). This includes the blindness test:
  `--mode asimov` refuses to read data skims (exit 7, `CMSFIT_BLIND_ABORT`).
- `CMSFIT_DONE` reached in 33.3 s. Gates G1, G2, G3 and G4 pass on the synthetic skims.
- Injection-recovery: with an injected toponium cross-section of 8.8 pb, the Asimov fits return
  mu_hat = 8.800 for all three spin targets (singlet, triplet, unchanged) and 6.400 at 6.4 pb, with the
  top-pair normalisation k_tt = 1.0000. The data-like fit on synthetic pseudo-data returns 8.817 pb
  (singlet; 68 % interval 1.8-15.9 pb). G4 fit-stability spread 0.00000.

What this does NOT say: the selftest templates are invented by `synthetic_skims()` to exercise the code.
Its Z_A values (1.26 at 8.8 pb, 0.92 at 6.4 pb; 3.31 statistics-only) and the "LOW" H4 label it prints
are properties of those invented templates, not of the CMS 2016 open data, and are not a grade of H4.
H4 is graded only by `cmsfit.py --mode asimov` on the six real simulation skims.

Result for the plan: the self-test condition in PREREG-CMSTOPO-v0.1 is met, so cmsfit 0.1.0 stays frozen
(no 0.1.1 needed). Next: when the last two simulation skims (tt part 2, tW) return, run
`cmsfit-asimov-v1` and grade H4; then the data skims, `cmsfit-data-v1` and H1-H3.

Made with AI agents working in MV2's name. Apache-2.0 (code), CC BY 4.0 (text).
