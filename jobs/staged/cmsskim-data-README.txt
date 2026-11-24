STAGED, NOT FILED. PREREG-CMSTOPO-v0.1 section 10: the data skims are filed only AFTER H4 (expected
sensitivity, cmsfit-asimov-v1) is graded and its result committed to munimtechnologies/cms-toponium-opendata.

To file: copy each <x>.job.json unchanged to <run>/compete/cmsskim-data-<g|h>-p<0|1>/job.json.
The workdir points at the frozen cmsskim 0.1.0 folder E:/MV2-Pursuit/work/20261124T162325Z-science/cmstopo
(the same bytes as commit 803a4af536 of the repo; never edit that folder).

On return: commit out/{skim-data_*.json, skim-data_*.txt, run.txt, tests.txt, sha256.txt, versions.txt} +
job.json to results/<job name>/. Gate G2 = >= 95 % of each MuonEG record's files processed, both parts.
Exit 5 / GENFAIL is a simulation-only failure mode; for data a read failure on all URL forms -> xrootd
(pip install xrootd fsspec-xrootd) v2 + amendment, never a threshold change.
Then: cmsfit --mode data with a manifest of all 6 MC skims + 4 data skims (pinned by git blob sha1),
grade H1-H3, publish whatever they say.
