# PREREG-CMSTOPO v0.1: toponium-like threshold excess in the public CMS 2016 G-H data (e-mu)

MV2 of Munim, Inc. Frozen 2026-11-24 by AI agents working in MV2's name. Code Apache-2.0, text
CC BY 4.0. Not affiliated with the CMS Collaboration. No human expert has checked this yet.

**Blindness statement.** When this file was committed, no mass, angle or data/simulation
comparison of the CMS open data had been looked at by this project (the only job filed before it,
`cmsprobe-v1`, reads plumbing only). Order of operations (section 10) keeps it that way: the
simulation skims run first, the fit code is published and its expected sensitivity graded on
simulation only, and only then is the data skim filed.

## 1. Question

CMS ([arXiv:2503.22382](https://arxiv.org/abs/2503.22382)) measured a pseudoscalar,
toponium-like excess at the top-pair threshold in 138 fb-1 of dilepton events:
sigma(eta_t) = 8.8 +- 0.5 (stat) +1.1 -1.3 (syst) pb, above five standard deviations; ATLAS and the
CMS lepton+jets channel ([arXiv:2608.19546](https://arxiv.org/abs/2608.19546)) confirmed it.
**Can it be seen by an outside group in the public CMS 2016 runs G-H (16.39 fb-1), with only
public files and every rule fixed in advance, and does its spin structure look like a spin
singlet?** The wow, if it works: the first detection of the top-pair threshold excess outside the
LHC collaborations, from public data. An honest "not seen at this luminosity" with a stated
expected sensitivity is also a result.

Statistics-only scaling of the CMS error (x sqrt(138/16.4) = 2.9, e-mu only, a cruder mass
proxy) suggests a statistical error near 3 pb on 8.8 pb; systematic terms we cannot constrain as
well as CMS will widen it. Section 8 states our predictions.

## 2. Data and simulation (CERN Open Data portal, checked 2026-11-24)

| Role | Record | Dataset (NanoAODv9, UL2016 post-VFP) | Events |
|---|---|---|---|
| data | 30528 | /MuonEG/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD | 33,854,612 |
| data | 30561 | /MuonEG/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD | 29,236,516 |
| tt (Powheg hvq + Pythia8 CP5, mt 172.5) | 67801 | TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8 | 43,546,000 |
| tW | 64895 | ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8 | 3,368,375 |
| tW-bar | 64839 | ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8 | 3,654,510 |
| other (DY) | 35669 | DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8 (first 6 of 41 sorted files) | 71,839,442 |

- Validated runs: record 14220, `Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt`.
- Luminosity: record 1059 brilcalc summaries, Run2016G 7.653261 fb-1 + Run2016H 8.740119 fb-1
  = **L = 16393.4 pb-1**, uncertainty 1.2 % (irrelevant here: the top-pair normalisation is free).
- Cross sections: sigma(tt) = 833.9 pb (the CMS value), B(W -> l nu) = 0.1086 per flavour, so
  sigma(TTTo2L2Nu) = 833.9 x 0.3258^2 = 88.51 pb; tW and tW-bar 19.47 pb each; DY M>50 6077.22 pb.
- Simulation weights: genWeight; normalisation to the sum of `genEventSumw` (Runs tree) over the
  processed files. Skims use only whole files.

## 3. Selection (code: `cmsskim.py` 0.1.0, `select_reco`)

1. Trigger: OR of the four `HLT_Mu23_..._Ele12_...` / `HLT_Mu8_..._Ele23_...` e-mu paths (with
   and without `_DZ`), those present in the file. Data from MuonEG only (no overlap removal needed).
2. Data: run/lumi in the validated-runs list. All events: `PV_npvsGood >= 1`; the noise flags
   goodVertices, globalSuperTightHalo2016Filter, HBHENoiseFilter, HBHENoiseIsoFilter,
   EcalDeadCellTriggerPrimitiveFilter, BadPFMuonFilter, BadPFMuonDzFilter, eeBadScFilter (those
   present).
3. Electrons: cutBased >= 4 (tight), pT > 20 GeV, |eta| < 2.4, |eta + deltaEtaSC| outside
   (1.4442, 1.5660). Muons: tightId, pfRelIso04_all < 0.15, pT > 20 GeV, |eta| < 2.4.
   Exactly one such electron and one such muon (third-lepton veto), opposite charge, leading
   lepton pT > 25 GeV, m(e mu) > 20 GeV.
4. Jets: jetId bit 2 (tight), pT > 30 GeV, |eta| < 2.4, Delta R > 0.4 from both leptons. At least
   two jets; at least one with DeepJet `btagDeepFlavB > 0.2489` (UL2016 post-VFP medium).
5. The two jets used: b-tagged jets first (by pT), then the hardest other jets.

## 4. Observables and binning

- m(llbb): invariant mass of electron + muon + the two jets. Edges
  **[0, 160, 200, 240, 280, 330, 400, 500, 650, inf) GeV** (9 bins).
- cos(phi*): cosine of the opening angle of the electron and muon in the rest frame of the llbb
  system, a reconstruction-free stand-in for c_hel near threshold, where the tops are slow. Edges
  **[-1, -0.5, 0, 0.5, 1]** (4 bins). Fit histogram: 9 x 4 = 36 bins.
- Information only: |Delta phi(e, mu)| in the lab, 10 bins on [0, pi].

## 5. Templates

**Top pairs.** Every simulated tt event gets a generator class from the hard-process top, antitop
and W-decay leptons (e, mu or tau; per-event fallback to last-copy tops and first-copy leptons
from a W): generator m(tt) < 350 GeV split into 40 equal bins of the generator spin angle cos(phi)
(l+ in the top rest frame, l- in the antitop rest frame, both reached through the tt rest frame;
classes 0-39), 350-400 GeV (class 40), >= 400 GeV (class 41), no generator tops (class 42).
Sums of weights are kept per class before selection (G_k) and after selection per fit bin (H_k,i),
for the nominal weight and each variation of section 6. The top-pair background T_i is the sum of
all classes, normalised to L x 88.51 pb / sum(genEventSumw).

**Toponium-like signal (no dedicated sample is public).** From the window classes k = 0..39
(generator m(tt) < 350 GeV, i.e. within about 5 GeV of 2 m_t, where CMS's eta_t at 343 GeV with a
2.8 GeV width sits after detector smearing of about 15 %), with c_k the class centre:
- D_win = -3 x sum_k G_k c_k / sum_k G_k (the window's own spin correlation, nominal weights);
- target spin distribution (1 - D_tgt c)/2, weights w_k = (1 - D_tgt c_k)/(1 - D_win c_k);
  **singlet D_tgt = -1 (primary)**, triplet-like D_tgt = +1/3 (the CMS chi_t slope), and
  unchanged w_k = 1;
- S_i per 1 pb = L x 1 pb x 0.10614564 x sum_k w_k H_k,i / sum_k w_k G_k.
The signal is added on top of the full top-pair simulation, as CMS adds eta_t to its
fixed-order sample. This reweighted-continuum template is an approximation of the CMS simplified
eta_t model (a LO MadGraph resonance at 343 GeV), named as such in every result.

**tW, other.** tW and tW-bar each normalised to L x 19.47 pb / sum(genEventSumw); "other" = DY
M>50 (6 files) normalised to L x 6077.22 pb / sum(genEventSumw) of those files; it stands in for
all small backgrounds (WW, ttV, non-prompt), hence its wide prior.

## 6. Fit (to be implemented as `cmsfit.py` 0.1.0 before any data skim is filed)

Expected count in bin i:
nu_i = gamma_i x [ k_tt x ( T_i(theta) + mu x S_i(theta_JES) ) + 1.15^theta_tW x TW_i(theta_JES)
+ 1.5^theta_oth x O_i(theta_JES) ].

- Free: mu = sigma(eta_t) in pb, mu >= 0 (bounded at 60); k_tt in [0.3, 3]. Because k_tt
  multiplies the signal too, efficiency and scale-factor differences cancel in mu.
- Gaussian unit nuisances: theta_muR, theta_muF (LHEScaleWeight), theta_ISR, theta_FSR (PSWeight)
  on the top-pair template, each **shape only** (varied template rescaled to the nominal total);
  theta_JES (jet pT and mass x (1 +- 0.03), a proxy for the full JES) on all templates with its
  rate effect; **theta_thr: all top-pair classes with generator m(tt) < 400 GeV scaled by
  (1 + 0.03 theta_thr)**, a conservative stand-in for the threshold-region effects CMS models with
  the Yukawa-coupling and electroweak-scheme nuisances; theta_tW (15 %); theta_oth (50 %).
- Interpolation: per-bin piecewise linear between the -1, 0, +1 templates; yields floored at 1e-9.
- Simulation statistics: Barlow-Beeston-lite, one gamma_i per bin with relative width
  delta_i = sqrt(sum w^2)/sum w of the total nominal prediction, profiled analytically
  (gamma = [-(nu delta^2 - 1) + sqrt((nu delta^2 - 1)^2 + 4 n delta^2)] / 2).
- NLL = sum_i [nu_i - n_i ln nu_i + (gamma_i - 1)^2 / (2 delta_i^2)] + sum theta^2 / 2.
  Minimiser: scipy L-BFGS-B from 5 starts (mu = 0, 4, 8.8, 15, 25); the best is kept.
- Significance: q0 = 2 [NLL(mu = 0) - NLL(mu_hat)] for mu_hat > 0 (else 0), Z = sqrt(q0).
- Interval: profile scan of mu in 0.1 pb steps on [0, 40]; 68 % interval where 2 Delta NLL <= 1.
- Expected sensitivity: Asimov data = nominal prediction with mu = 8.8 pb, k_tt = 1, theta = 0,
  gamma = 1 (not rounded); Z_A from the same q0. Also reported for mu = 6.4 pb (NRQCD estimate
  quoted by CMS).

## 7. Gates (all must pass, else NOT GRADABLE and the reason is published)

- **G1 simulation:** >= 8,000,000 tt events processed in total; the first-chunk generator check
  `has_gen >= 0.95` in every tt part; for every processed file |sum genWeight - genEventSumw| /
  genEventSumw < 1e-6; at least one file of each tW record and of DY.
- **G2 data:** at least 95 % of the files of each MuonEG record processed.
- **G3 normalisation:** fitted k_tt in [0.7, 1.3].
- **G4 fit health:** the 5 starts agree on the minimum NLL within 0.01 for both the background-only
  and the signal+background fit.

## 8. Hypotheses, grades and our predictions (fixed now)

| | Quantity | Grades | Our prior |
|---|---|---|---|
| H4 (graded on simulation, before data) | Z_A at 8.8 pb | HIGH >= 3; MEDIUM 1.5-3; LOW < 1.5 | 20 / 55 / 25 % |
| H1 detection | Z (singlet template) | EVIDENCE >= 3; HINT 2-3; NOT SEEN < 2 | 15 / 30 / 55 % |
| H2 size | mu_hat vs CMS 8.8 pb | CONSISTENT if \|mu_hat - 8.8\| <= 2 sqrt(s^2 + 1.3^2) (s = our half-width of the 68 % interval), else HIGHER / LOWER | 85 / 10 / 5 % |
| H3 spin | min NLL: singlet vs unchanged and vs triplet-like | SINGLET PREFERRED if 2 Delta NLL >= 1 against both; NOT SINGLET if worse by >= 1 than either; else UNRESOLVED | 30 / 10 / 60 % |

Information only (never graded): the fit without theta_thr; statistics-only mu; per-bin pulls;
fitted nuisances; |Delta phi(e, mu)| shapes; the same fit with the unchanged and triplet-like
templates as the signal.

## 9. Limits known now

e-mu only and post-VFP 2016 only (16.4 of CMS's 138 fb-1); m(llbb) and cos(phi*) are proxies for
m(tt) and c_hel; no lepton, trigger or b-tag scale factors and no pile-up reweighting (rate
effects absorbed by k_tt, shape effects not modelled); JES is a flat proxy; no top-mass, hdamp,
PDF, alternative-generator, Yukawa or electroweak terms beyond theta_thr; the signal is a
reweighted Powheg window, not a bound-state sample. A positive result would show that the
threshold excess is visible in public data under these approximations; it would not be a
measurement competitive with CMS.

## 10. Order of operations and version guard

1. Simulation skims: `cmsskim-tt-p{0,1,2,3}-v1`, `cmsskim-tw-v1`, `cmsskim-dy-v1` (filed with this
   pre-registration). Histograms committed to `results/` as returned.
2. `cmsfit.py` 0.1.0 + tests implementing section 6 exactly, published; then a fit job on the
   simulation histograms only: H4 graded and published.
3. Only then the data skims `cmsskim-data_G-v1`, `cmsskim-data_H-v1`; then the data fit; H1-H3
   graded and published whatever they say.
Frozen: `cmsskim.py` 0.1.0 and this file. Any change before the data fit is a numbered amendment
(AMENDMENT-CMSTOPO-0.1.x) committed before the data fit runs; nothing changes after it.

## 11. Prior work (searched 2026-11-24, last at 16:40 UTC)

Searches: web "toponium CMS open data 2016 NanoAOD threshold excess independent reanalysis";
web "toponium open data ttbar threshold excess independent analysis arXiv 2026"; web "CMS open
data 2016 toponium threshold excess dilepton reanalysis"; web "ATLAS open data 2015-2016 ttbar
spin correlation excess reanalysis toponium threshold arXiv"; arXiv abstracts of the papers below.
Nothing found that analyses the public collision files for the threshold excess. Closest:
- CMS, [arXiv:2503.22382](https://arxiv.org/abs/2503.22382) and
  [arXiv:2608.19546](https://arxiv.org/abs/2608.19546): the measurements, collaboration data and tools.
- Fuks, Hossain and Keaveney, [arXiv:2511.02040](https://arxiv.org/abs/2511.02040): reinterpretation
  from HEPData tables, not collision data.
- Fuks, Hagiwara, Ma and Zheng, [arXiv:2102.11281](https://arxiv.org/abs/2102.11281): the
  simplified eta_t model and spin signatures (theory).
- Aushev et al., [arXiv:2611.09797](https://arxiv.org/abs/2611.09797): top-pair reconstruction
  on CMS open data, no threshold test.
- [arXiv:2610.42143](https://arxiv.org/abs/2610.42143): strange-quark polarisation in top decays
  with the CMS 2016 G-H open data (single-lepton), not toponium.
- Baer and Barger, [arXiv:2611.22421](https://arxiv.org/abs/2611.22421): review; no Run 3 or
  outside threshold analysis exists.
- MV2's own ATLAS open-data test of the 2019 spin-correlation excess:
  https://github.com/munimtechnologies/ttbar-threshold-opendata (generator spread at threshold
  exceeded the excess there; this study adds a toponium template and a likelihood fit).
