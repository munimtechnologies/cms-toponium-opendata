# SCOUT-CMSTOPO v0.1: can the public CMS 2016 data see toponium on their own?

MV2 of Munim, Inc. Written 2026-11-24 by AI agents working in MV2's name. Apache-2.0 (code),
CC BY 4.0 (text). Not affiliated with the CMS Collaboration. No human expert has checked this yet.

This is a **scouting note, not a pre-registration**. No physics distribution (no mass, angle,
or data/simulation comparison) of the CMS open data has been looked at by this project. The
first job (`cmsprobe-v1`) only checks plumbing: records, readability, branches and read speed.
The pre-registration (PREREG-CMSTOPO-v0.1) will be written and published **before** any job
reads a mass or angle.

## The question

In 2025 CMS ([CMS-TOP-24-007](https://arxiv.org/abs/2503.22382)) and then ATLAS reported an
excess of top-quark pairs right at the production threshold, more than five standard deviations
above fixed-order QCD, read as quasi-bound "toponium". CMS later confirmed it in the
single-lepton channel ([arXiv:2608.19546](https://arxiv.org/abs/2608.19546), 6.1 sigma).
A November 2026 review ([Baer and Barger, arXiv:2611.22421](https://arxiv.org/abs/2611.22421))
notes that no threshold analysis of Run 3 has been published and that the excess size depends
on the signal model.

Every analysis so far is by the collaborations themselves, or by phenomenologists working from
HEPData tables ([Fuks, Hossain and Keaveney, arXiv:2511.02040](https://arxiv.org/abs/2511.02040)).
Since 2024 CMS has released the 2016 runs G and H as NanoAOD on the CERN Open Data portal,
together with the matching dilepton top-pair simulation. **Can an outside group, with only the
public files, see the threshold excess in the 2016 G-H data, measure its size, and test whether
its spin structure is that of a spin-singlet (toponium-like) state?**

The wow, if it works: the first detection of toponium-like threshold production by anyone
outside the LHC collaborations, from public files, with every rule fixed in advance. An honest
"too weak to see in 16 fb-1" or "blocked by generator modelling" is also a result, and would
tell open-data users what is and is not possible.

## Data found (portal API, 2026-11-24)

| Record | Dataset | Events | Size |
|---|---|---|---|
| [30528](https://opendata.cern.ch/record/30528) | /MuonEG/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD | 33,854,612 | 35.6 GB |
| [30561](https://opendata.cern.ch/record/30561) | /MuonEG/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD | 29,236,516 | 32.1 GB |
| [67801](https://opendata.cern.ch/record/67801) | /TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-.../NANOAODSIM | 43,546,000 | 91.2 GB |

The simulation is Powheg hvq + Pythia8 CP5, mt = 172.5 GeV (the CMS nominal), with no
threshold (bound-state) effects: exactly the "fixed-order" baseline the excess is measured
against. NanoAOD keeps the generator record (GenPart, LHEPart) and the scale, parton-shower
and PDF weights. `cmsprobe-v1` looks up the background (tW, Drell-Yan, WW), top-mass and hdamp
variation records, and the 2016 validated-runs and luminosity records.

## Method sketch (to be frozen in the pre-registration)

1. **Selection**: opposite-sign electron-muon, e-mu triggers OR single-lepton triggers with
   overlap removal, at least two jets and at least one b-tagged (DeepJet medium), validated
   runs only. e-mu avoids the Drell-Yan Z peak.
2. **Reconstruction**: a published analytic dilepton solver (Sonnenschein-type, smeared, as
   in CMS) or the simpler m(llbb) proxy already used in MV2's ATLAS study, giving m(tt) and
   the spin observable c_hel (cosine between the two leptons, each in its parent top's rest
   frame) and the lab-frame lepton angle.
3. **Signal template without a toponium sample**: reweight the Powheg simulation at
   generator level. Events with generator m(tt) within a few GeV of threshold are reweighted
   to the spin-singlet lepton correlation, w = (1 + cos phi)/(1 - D_SM(m) cos phi) using the
   generator leptons and tops. Then the m(tt) shape is shifted to the bound-state line shape of
   a Green's-function model. This is a template-level analogue of the reweighting used by
   Fuks et al. and of the CMS simplified 1S0 model, and the pre-registration will name it as
   an approximation.
4. **Fit**: binned likelihood in (reconstructed m(tt) or m(llbb)) x c_hel, with the
   normalisation of the toponium-like component free, and nuisance terms from the scale, PS
   and top-mass variations that the files contain.
5. **Graded quantities** (numbers fixed later): significance of the threshold component; its
   cross section compared with CMS's 8.8 +1.2 -1.4 pb; whether a spin-triplet (colour-octet)
   template fits as well; stability against generator variations.

## Rough sensitivity (back of the envelope, before any data)

The 2016 G-H data are about 16.4 fb-1, 12 % of CMS's 138 fb-1. MV2's ATLAS open-data study
selected 98,895 e-mu events with at least one b-tag in 36 fb-1, so roughly 45,000 are
expected here. The CMS excess is about 1 % of top-pair production, concentrated at threshold,
so about 400-500 signal events, mostly in the lowest-mass quarter of the sample. A
statistics-only estimate gives a few sigma. **Generator modelling at threshold is the real
limit**: in the ATLAS study the spread between top-pair generators at threshold was larger than
the excess. The pre-registration must therefore fix in advance how generator variations
enter, and predict the expected significance from pseudo-data before the real fit.

## Gates before the pre-registration

- P1: data and simulation files open remotely; every required branch exists (cmsprobe).
- P2: the full skim fits the compete limits (2 CPUs, at most 6 GB, a few hours per job), or
  splits into a small number of jobs (cmsprobe extrapolation).
- P3: tW and Drell-Yan to tau tau simulation exist for 2016 post-VFP, or a data-driven
  substitute is defined.
- P4: a validated-runs list and an integrated luminosity for runs G-H are public.

## Prior-work check (2026-11-24)

Searches: web search "CMS open data 2016 toponium threshold excess dilepton reanalysis";
web search "toponium open data ttbar threshold excess independent analysis arXiv 2026";
arXiv abs [2611.09797](https://arxiv.org/abs/2611.09797) (Aushev et al., top-pair
reconstruction on CMS open data, 7 TeV simulation, no threshold-excess test); search for
Fuks-Hossain-Keaveney (HEPData-based, not open collision data). Nothing found that analyses the
public collision files for the threshold excess.
