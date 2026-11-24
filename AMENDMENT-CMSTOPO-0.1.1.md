# AMENDMENT-CMSTOPO-0.1.1: fit implementation details (cmsfit 0.1.0)

MV2 of Munim, Inc. Written 2026-11-24 by AI agents working in MV2's name, before any simulation
skim of PREREG-CMSTOPO-v0.1 has returned and before any data skim has been filed. Text CC BY 4.0,
code Apache-2.0. No human expert has checked this yet.

Section 6 of [PREREG-CMSTOPO-v0.1](PREREG-CMSTOPO-v0.1.md) fixes the fit model, the minimiser,
the significance, the interval and the Asimov data. Writing the code (`cmsfit.py` 0.1.0) needed
the choices below, which v0.1 did not spell out. **No selection, binning, template, nuisance size,
gate threshold, grade boundary or prediction changes.**

(a) **Nuisance bounds.** Every Gaussian theta is bounded to [-5, 5] in the minimiser (mu to
[0, 60] and k_tt to [0.3, 3], as pre-registered).

(b) **Background-only starts.** The five pre-registered starts (mu = 0, 4, 8.8, 15, 25) are
signal-strength starts and cannot differ when mu is fixed at 0. The background-only fit therefore
starts from k_tt = 0.8, 0.9, 1.0, 1.1, 1.2 (all thetas 0); gate G4 compares the minima of these
five as well. Fixed parameters are removed from the minimiser rather than bounded.

(c) **theta_thr.** It multiplies the interpolated top-pair template of the generator classes with
m(tt) < 400 GeV (classes 0-40) by (1 + 0.03 theta_thr). It does not act on the signal term
mu x S, which is added separately as in section 6.

(d) **Shape-only variations.** For muR, muF, ISR and FSR the rescaling factor (nominal total over
varied total) is computed on the whole top-pair template (all classes) and applied to both its
m(tt) < 400 GeV and >= 400 GeV parts, so (c) acts on consistently normalised pieces.

(e) **Barlow-Beeston-lite widths.** delta_i is computed from the nominal background-only
prediction (top pair + tW + other, k_tt = 1, theta = 0): sqrt(sum over samples of norm^2 x sum
genWeight^2) over the nominal yield. The signal template is a reweighted subset of the same
top-pair simulation, so its events are not counted a second time.

(f) **NLL form.** The NLL is written as the Poisson deviance, sum [nu g - n + n ln(n / (nu g))],
plus the pre-registered gamma and theta penalties. This differs from section 6 by a constant in n
only; every NLL difference (q0, the scan, H3, G4) is unchanged.

(g) **Profile scan.** At each 0.1 pb grid point the better of two fits is kept: one started from
the previous grid point's solution and one from the global best fit. The reference minimum is the
smaller of the global fit and the lowest scan point. The 68 % interval is the lowest and highest
grid point with 2 Delta NLL <= 1 (no interpolation); if the upper end reaches 40 pb, this is
reported.

(h) **Asimov.** Asimov data for each signal template are built with that template; H4 is graded on
the singlet template at 8.8 pb, as pre-registered. Asimov mode refuses to run if any data skim is
among its inputs.

(i) **Inputs.** The fit reads only skim JSON files committed to this repository under `results/`,
each downloaded at a fixed commit and checked against its git blob sha1 before use; any mismatch
aborts the job.

(j) **Information only.** Per-bin pulls are (n - nu) / sqrt(nu) at the best fit, without gamma.
The statistics-only fit fixes every theta at 0 and drops the gamma terms (k_tt and mu free).

`cmsfit.py` 0.1.0 and this amendment are frozen together. The machinery is checked first on
invented synthetic templates (`--mode selftest`, never graded), then on the simulation histograms
(`--mode asimov`, grades H4), and only then on data.
