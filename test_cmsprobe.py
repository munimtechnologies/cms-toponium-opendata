# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Unit tests for the pure helpers of cmsprobe (standard library only)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmsprobe as cp  # noqa: E402

U = ("root://eospublic.cern.ch//eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/"
     "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/120000/A.root")


class TestHelpers(unittest.TestCase):
    def test_extract_nested_unique_in_order(self):
        obj = {"metadata": {"_file_indices": [
            {"key": "k1", "files": [{"uri": U}, {"uri": U.replace("A.root", "B.root")}]},
            {"key": "k2", "files": [{"uri": U}, {"uri": "root://eospublic.cern.ch//eos/opendata/x.json"}]},
        ]}, "other": ["https://example.org/C.root", 3, None]}
        self.assertEqual(cp.extract_root_uris(obj), [U, U.replace("A.root", "B.root")])

    def test_extract_empty(self):
        self.assertEqual(cp.extract_root_uris({}), [])

    def test_url_forms(self):
        forms = cp.url_forms(U)
        self.assertEqual(len(forms), 3)
        self.assertEqual(forms[0], "https://opendata.cern.ch/eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/"
                         "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/"
                         "106X_mcRun2_asymptotic_v17-v1/120000/A.root")
        self.assertTrue(forms[1].startswith("https://eospublic.cern.ch:443//eos/opendata/cms/mc/"))
        self.assertEqual(forms[2], U)

    def test_url_forms_rejects_other_hosts(self):
        with self.assertRaises(ValueError):
            cp.url_forms("root://other.cern.ch//eos/x.root")

    def test_extrapolate(self):
        # 100000 entries in 50 s -> 2000/s; 43,200,000 entries -> 21600 s = 6 h; 2 workers -> 3 h
        self.assertAlmostEqual(cp.extrapolate_hours(50, 100000, 43200000), 6.0)
        self.assertAlmostEqual(cp.extrapolate_hours(50, 100000, 43200000, workers=2), 3.0)
        self.assertIsNone(cp.extrapolate_hours(0, 100, 1000))
        self.assertIsNone(cp.extrapolate_hours(10, 0, 1000))

    def test_missing(self):
        self.assertEqual(cp.missing(["a", "b", "c"], ["c", "a"]), ["b"])

    def test_record_summary(self):
        meta = {"recid": 67801, "title": "T", "doi": "d",
                "distribution": {"number_events": 43546000, "number_files": 49, "size": 91175567890}}
        self.assertEqual(cp.record_summary(meta), {"recid": 67801, "title": "T", "number_events": 43546000,
                                                   "number_files": 49, "size": 91175567890, "doi": "d"})
        self.assertEqual(cp.record_summary({})["number_events"], None)

    def test_branch_lists_unique(self):
        allb = cp.COMMON_BRANCHES + cp.MC_BRANCHES
        self.assertEqual(len(allb), len(set(allb)))


if __name__ == "__main__":
    unittest.main()
