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

"""Unit tests for visualization components."""

from __future__ import annotations

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

import generate_dashboard
import visualize_bioactivity
import visualize_compound
import visualize_similarity


class TestVisualizations(unittest.TestCase):
  """Test suite for analytical and visualization modules."""

  def test_lipinski_evaluation_compliant(self):
    """Verify Ro5 compliance scoring for drug-like molecule."""
    props = {
        "full_mwt": 180.16,
        "alogp": 1.31,
        "hbd": 1,
        "hba": 3,
        "rtb": 2,
        "psa": 63.6,
    }
    metrics = visualize_compound.evaluate_lipinski(props)
    self.assertEqual(metrics["violations"], 0)
    self.assertIn("Compliant", metrics["badge"])

  def test_lipinski_evaluation_violations(self):
    """Verify Ro5 violation counting."""
    props = {
        "full_mwt": 620.5,  # >500
        "alogp": 6.2,       # >5.0
        "hbd": 6,          # >5
        "hba": 12,         # >10
    }
    metrics = visualize_compound.evaluate_lipinski(props)
    self.assertEqual(metrics["violations"], 4)
    self.assertIn("High Violation", metrics["badge"])

  def test_radar_svg_generation(self):
    """Verify SVG radar chart markup."""
    props = {"full_mwt": 180.16, "alogp": 1.31, "hbd": 1, "hba": 3, "rtb": 2, "psa": 63.6}
    metrics = visualize_compound.evaluate_lipinski(props)
    svg = visualize_compound.generate_radar_svg(metrics)
    self.assertIn("<svg", svg)
    self.assertIn("</svg>", svg)
    self.assertIn("polygon", svg)
    self.assertIn("Ro5 Boundary", svg)

  def test_pic50_calculation(self):
    """Verify pIC50 mathematical calculation and edge cases."""
    self.assertEqual(visualize_bioactivity.calculate_pic50(1.0), 9.0)
    self.assertEqual(visualize_bioactivity.calculate_pic50(10.0), 8.0)
    self.assertEqual(visualize_bioactivity.calculate_pic50(100.0), 7.0)
    self.assertEqual(visualize_bioactivity.calculate_pic50(1000.0), 6.0)
    self.assertIsNone(visualize_bioactivity.calculate_pic50(0.0))
    self.assertIsNone(visualize_bioactivity.calculate_pic50(-5.0))
    self.assertIsNone(visualize_bioactivity.calculate_pic50(None))

  def test_bioactivity_processing(self):
    """Verify bioactivity binning and statistics."""
    acts = [
        {"standard_value": "10.0", "standard_units": "nM", "standard_type": "IC50"},
        {"standard_value": "500.0", "standard_units": "nM", "standard_type": "IC50"},
        {"standard_value": "20000.0", "standard_units": "nM", "standard_type": "IC50"},
    ]
    res = visualize_bioactivity.process_activities(acts)
    self.assertEqual(res["stats"]["count"], 3)
    self.assertEqual(res["potency_counts"]["high"], 1)
    self.assertEqual(res["potency_counts"]["moderate"], 1)
    self.assertEqual(res["potency_counts"]["weak"], 1)

  def test_similarity_gallery(self):
    """Verify similarity gallery rendering."""
    sim_data = {
        "molecules": [
            {
                "molecule_chembl_id": "CHEMBL26",
                "pref_name": "SALICYLIC ACID",
                "similarity": 91.2,
                "molecule_structures": {"canonical_smiles": "O=C(O)c1ccccc1O"},
                "molecule_properties": {"full_mwt": 138.12, "alogp": 1.48},
            }
        ]
    }
    html_out = visualize_similarity.render_similarity_gallery(sim_data)
    self.assertIn("CHEMBL26", html_out)
    self.assertIn("SALICYLIC ACID", html_out)
    self.assertIn("91.2% Match", html_out)

  def test_dashboard_compiler_physics_and_theme(self):
    """Verify interactive dashboard contains physics components and theme toggle."""
    dashboard = generate_dashboard.build_full_dashboard(
        title="Test Aspirin",
        compound_data={"molecule_chembl_id": "CHEMBL25", "pref_name": "ASPIRIN"},
        similarity_data={"molecules": [{"molecule_chembl_id": "CHEMBL26", "similarity": 91.2}]},
    )
    # Verify physics simulation elements
    self.assertIn("physics-viewport", dashboard)
    self.assertIn("stepPhysics", dashboard)
    self.assertIn("Force-Directed Similarity & SAR Network", dashboard)

    # Verify theme toggle
    self.assertIn("theme-toggle-btn", dashboard)
    self.assertIn("theme-light-label", dashboard)
    self.assertIn("theme-dark-label", dashboard)

    # Verify 3D conformer viewer
    self.assertIn("mol-3d-viewer", dashboard)
    self.assertIn("3Dmol-min.js", dashboard)


if __name__ == "__main__":
  unittest.main()
