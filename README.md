# cms-toponium-opendata

An open, pre-registered attempt by MV2 of Munim, Inc. to see the top-quark-pair threshold
excess ("toponium") in the **public** CMS 2016 collision data (runs G and H, NanoAOD, CERN Open
Data portal), using only public files and rules fixed before any distribution is looked at.

Made with AI agents working in MV2's name. No human expert has checked it yet. This project is
not affiliated with the CMS Collaboration. Code: Apache-2.0. Text: CC BY 4.0.

## Status (2026-11-24)

- `SCOUT-CMSTOPO-v0.1.md`: the question, the public records found, a method sketch and the
  gates that must pass before the pre-registration is written. **Not a pre-registration.**
- `cmsprobe.py` (0.1.0) + `test_cmsprobe.py`: a plumbing probe. It finds the records, opens
  one data and one simulation file remotely, checks that every branch the analysis needs exists
  and times the read. It looks at no mass or angle. Job definition: `jobs/cmsprobe-v1.json`.

## Related

- MV2's ATLAS open-data study of the 2019 top spin-correlation excess:
  https://github.com/munimtechnologies/ttbar-threshold-opendata
- CMS, Observation of a pseudoscalar excess at the top quark pair production threshold,
  https://arxiv.org/abs/2503.22382

## License

Apache-2.0 (see LICENSE).
