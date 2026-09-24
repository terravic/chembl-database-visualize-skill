# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Unit tests for chembl_api.py."""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

import chembl_api


class TestChemblApi(unittest.TestCase):
  """Test suite for ChEMBL API client and mock engine."""

  def test_unit_normalization(self):
    """Verify concentration unit conversion to nM."""
    data = {
        "activities": [
            {"standard_value": "1.0", "standard_units": "uM"},
            {"standard_value": "0.5", "standard_units": "mM"},
            {"standard_value": "250", "standard_units": "nM"},
            {"standard_value": "1000", "standard_units": "pM"},
            {"standard_value": "invalid", "standard_units": "uM"},
            {"standard_value": None, "standard_units": "nM"},
        ]
    }
    normalized = chembl_api._normalize_activity_records(data)
    acts = normalized["activities"]
    self.assertEqual(acts[0]["normalized_value_nM"], 1000.0)
    self.assertEqual(acts[1]["normalized_value_nM"], 500000.0)
    self.assertEqual(acts[2]["normalized_value_nM"], 250.0)
    self.assertEqual(acts[3]["normalized_value_nM"], 1.0)
    self.assertIsNone(acts[4]["normalized_value_nM"])
    self.assertIsNone(acts[5]["normalized_value_nM"])

  def test_url_builder(self):
    """Verify REST URL construction."""
    url1 = chembl_api._build_url("molecule", item_id="CHEMBL25")
    self.assertEqual(url1, "https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL25.json")

    url2 = chembl_api._build_url("activity", limit=10, filters=["standard_type=IC50"])
    self.assertIn("limit=10", url2)
    self.assertIn("standard_type=IC50", url2)

  def test_mock_status(self):
    """Verify status response in mock mode."""
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tf:
      tmp_path = tf.name
    try:
      parser = chembl_api.build_parser()
      args = parser.parse_args(["status", "--output", tmp_path, "--mock"])
      chembl_api.cmd_status(args)

      with open(tmp_path, "r", encoding="utf-8") as f:
        res = json.load(f)
      self.assertEqual(res.get("status"), "UP")
      self.assertTrue(res.get("mock"))
      self.assertIn("_license_notice", res)
    finally:
      if os.path.exists(tmp_path):
        os.remove(tmp_path)

  def test_mock_molecule(self):
    """Verify mock molecule fetch for CHEMBL25."""
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tf:
      tmp_path = tf.name
    try:
      parser = chembl_api.build_parser()
      args = parser.parse_args(["molecule", "--id", "CHEMBL25", "--output", tmp_path, "--mock"])
      chembl_api.cmd_generic(args)

      with open(tmp_path, "r", encoding="utf-8") as f:
        res = json.load(f)
      self.assertEqual(res.get("molecule_chembl_id"), "CHEMBL25")
      self.assertEqual(res.get("pref_name"), "ASPIRIN")
      self.assertIn("molecule_properties", res)
    finally:
      if os.path.exists(tmp_path):
        os.remove(tmp_path)


if __name__ == "__main__":
  unittest.main()
