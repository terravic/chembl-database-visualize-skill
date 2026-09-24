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

"""End-to-end integration tests for run_pipeline.py."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPELINE_SCRIPT = os.path.join(PROJECT_ROOT, "scripts", "run_pipeline.py")


class TestPipelineIntegration(unittest.TestCase):
  """End-to-end integration tests for the unified pipeline."""

  def setUp(self):
    self.temp_dir = tempfile.mkdtemp()

  def tearDown(self):
    if os.path.exists(self.temp_dir):
      shutil.rmtree(self.temp_dir)

  def test_molecule_mock_pipeline(self):
    """Test full compound profiling pipeline in mock mode."""
    cmd = [
        sys.executable,
        PIPELINE_SCRIPT,
        "--molecule", "CHEMBL25",
        "--mock",
        "-o", self.temp_dir,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")

    # Verify generated artifacts
    dashboard_path = os.path.join(self.temp_dir, "dashboard.html")
    data_path = os.path.join(self.temp_dir, "data.json")
    sdf_path = os.path.join(self.temp_dir, "structure.sdf")
    svg_path = os.path.join(self.temp_dir, "structure.svg")

    self.assertTrue(os.path.exists(dashboard_path))
    self.assertTrue(os.path.exists(data_path))
    self.assertTrue(os.path.exists(sdf_path))
    self.assertTrue(os.path.exists(svg_path))

    # Inspect dashboard contents
    with open(dashboard_path, "r", encoding="utf-8") as f:
      html_content = f.read()

    self.assertIn("tailwindcss", html_content)
    self.assertIn("3Dmol-min.js", html_content)
    self.assertIn("ASPIRIN", html_content)
    self.assertIn("CHEMBL25", html_content)
    self.assertIn("Zdrazil", html_content)
    self.assertIn("Davies", html_content)
    self.assertIn("Disclaimer:", html_content)

  def test_target_mock_pipeline(self):
    """Test target bioactivity profiling pipeline in mock mode."""
    cmd = [
        sys.executable,
        PIPELINE_SCRIPT,
        "--target", "CHEMBL203",
        "--activity_type", "IC50",
        "--mock",
        "-o", self.temp_dir,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    self.assertEqual(proc.returncode, 0, f"Target pipeline failed: {proc.stderr}")

    dashboard_path = os.path.join(self.temp_dir, "dashboard.html")
    self.assertTrue(os.path.exists(dashboard_path))

  def test_similarity_mock_pipeline(self):
    """Test chemical similarity pipeline in mock mode."""
    cmd = [
        sys.executable,
        PIPELINE_SCRIPT,
        "--smiles", "CC(=O)Oc1ccccc1C(=O)O",
        "--similarity", "85",
        "--mock",
        "-o", self.temp_dir,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    self.assertEqual(proc.returncode, 0, f"Similarity pipeline failed: {proc.stderr}")

    dashboard_path = os.path.join(self.temp_dir, "dashboard.html")
    self.assertTrue(os.path.exists(dashboard_path))

  def test_license_guard_creation(self):
    """Verify license guard verification file exists."""
    guard_file = os.path.join(PROJECT_ROOT, ".licenses", "chembl_database_visualize_LICENSE.txt")
    self.assertTrue(os.path.exists(guard_file))
    with open(guard_file, "r", encoding="utf-8") as f:
      content = f.read()
    self.assertIn("ChEMBL terms of use", content)
    self.assertIn("Timestamp:", content)


if __name__ == "__main__":
  unittest.main()
