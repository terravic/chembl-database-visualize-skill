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

"""Unit and runtime verification tests for ChEMBL dashboard components."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
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
    """Verify similarity gallery and heatmap rendering."""
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
    self.assertIn("Normalized Descriptor Heatmap Matrix", html_out)

  def test_dashboard_compiler_design_and_controls(self):
    """Verify interactive dashboard color tokens, node map, theme toggle, and zoom controls."""
    root_dir = os.path.join(os.path.dirname(__file__), "..")
    with open(os.path.join(root_dir, "sample_data", "chembl25_aspirin.json"), "r", encoding="utf-8") as f:
      compound_data = json.load(f)
    with open(os.path.join(root_dir, "sample_data", "chembl25_aspirin.sdf"), "r", encoding="utf-8") as f:
      sdf_content = f.read()
    with open(os.path.join(root_dir, "sample_data", "chembl25_aspirin.svg"), "r", encoding="utf-8") as f:
      svg_content = f.read()
    with open(os.path.join(root_dir, "sample_data", "egfr_activities.json"), "r", encoding="utf-8") as f:
      activity_data = json.load(f)
    with open(os.path.join(root_dir, "sample_data", "similarity_aspirin.json"), "r", encoding="utf-8") as f:
      similarity_data = json.load(f)

    dashboard = generate_dashboard.build_full_dashboard(
        title="ASPIRIN (CHEMBL25) - Profile & Analytics",
        compound_data=compound_data,
        sdf_content=sdf_content,
        svg_content=svg_content,
        activity_data=activity_data,
        similarity_data=similarity_data,
    )

    # Verify relationship node map elements
    self.assertIn("physics-viewport", dashboard)
    self.assertIn("stepPhysics", dashboard)
    self.assertIn("Interactive Relationship Node Map", dashboard)
    self.assertNotIn("Physics Graph", dashboard)
    self.assertNotIn("Pause Physics", dashboard)

    # Verify theme toggle
    self.assertIn("theme-toggle-btn", dashboard)
    self.assertIn("theme-light-label", dashboard)
    self.assertIn("theme-dark-label", dashboard)
    self.assertNotIn(">Light Mode<", dashboard)
    self.assertNotIn(">Dark Mode<", dashboard)

    # Verify Light & Dark theme color tokens
    for token in ("#FAFAF9", "#FFFFFF", "#E7E5E4", "#1C1917", "#787774", "#0D0E11", "#16171D", "#292A30", "#ECE8E2", "#9A9895"):
      self.assertIn(token, dashboard)

    # Verify Data Visualization Color Palette
    for pal_color in ("#0D9488", "#D97706", "#4F46E5", "#059669", "#E11D48"):
      self.assertIn(pal_color, dashboard)

    # Verify all interactive visual components have Zoom In, Zoom Out, and Reset View controls
    expected_controls = [
        "zoomView('chem2d', 1.25)",
        "zoomView('chem2d', 0.8)",
        "resetView('chem2d')",
        "zoom3D(1.25)",
        "zoom3D(0.8)",
        "reset3DView()",
        "zoomView('radar', 1.25)",
        "zoomView('radar', 0.8)",
        "resetView('radar')",
        "zoomNetwork(1.25)",
        "zoomNetwork(0.8)",
        "resetPhysicsPositions()",
        "zoomView('hist', 1.25)",
        "zoomView('hist', 0.8)",
        "resetView('hist')",
        "resetActivityFilter()",
        "resetAnalogFilters()",
    ]
    for ctrl in expected_controls:
      self.assertIn(ctrl, dashboard, f"Missing interactive control: {ctrl}")

    # Verify all 5 sections contain populated data
    for expected_data in (
        "CHEMBL25",
        "ASPIRIN",
        "CC(=O)Oc1ccccc1C(=O)O",
        "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        "Cyclooxygenase-1 inhibitor",
        "Cyclooxygenase-2 inhibitor",
        "Myocardial Infarction",
        "CHEMBL26",
        "SALICYLIC ACID",
        "CHEMBL112",
        "DIFLUNISAL",
    ):
      self.assertTrue(expected_data in dashboard, f"Missing expected data in dashboard: {expected_data}")

  def test_javascript_syntax_and_runtime_execution(self):
    """Verify embedded JS has zero syntax errors and all interactive functions execute properly."""
    root_dir = os.path.join(os.path.dirname(__file__), "..")
    dashboard_path = os.path.join(root_dir, "dashboard.html")
    self.assertTrue(os.path.exists(dashboard_path), "dashboard.html must exist in workspace root")

    with open(dashboard_path, "r", encoding="utf-8") as f:
      html_content = f.read()

    scripts = re.findall(r"<script>(.*?)</script>", html_content, re.S)
    self.assertEqual(len(scripts), 1, "Expected exactly 1 inline executable <script> block")
    js_code = scripts[0]

    # 1. Verify JS syntax via `node --check`
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False, encoding="utf-8") as tf:
      tf.write(js_code)
      js_file = tf.name

    try:
      check_res = subprocess.run(["node", "--check", js_file], capture_output=True, text=True, check=False)
      self.assertEqual(check_res.returncode, 0, f"JS syntax error: {check_res.stderr}")
    finally:
      if os.path.exists(js_file):
        os.remove(js_file)

    # 2. Execute JS inside a simulated sandboxed DOM (where localStorage throws SecurityError)
    sdf_match = re.search(r'<script id="sdf-data" type="text/plain">(.*?)</script>', html_content, re.S)
    json_match = re.search(r'<script id="dashboard-data" type="application/json">(.*?)</script>', html_content, re.S)
    self.assertIsNotNone(sdf_match)
    self.assertIsNotNone(json_match)

    harness_js = f"""
    const assert = require('assert');

    function makeEl(id, textContent = '') {{
      const attrs = {{}};
      const classes = new Set();
      const children = [];
      const childSvg = {{
        style: {{}},
        setAttribute(k, v) {{ attrs[k] = v; }},
        getAttribute(k) {{ return attrs[k]; }}
      }};
      return {{
        id,
        textContent,
        innerHTML: '',
        value: '',
        style: {{}},
        classList: {{
          toggle(c, force) {{
            if (force === undefined) {{
              if (classes.has(c)) {{ classes.delete(c); return false; }}
              classes.add(c); return true;
            }}
            if (force) classes.add(c); else classes.delete(c);
            return force;
          }},
          add(c) {{ classes.add(c); }},
          remove(c) {{ classes.delete(c); }},
          contains(c) {{ return classes.has(c); }}
        }},
        setAttribute(k, v) {{ attrs[k] = String(v); }},
        getAttribute(k) {{ return attrs[k]; }},
        appendChild(child) {{ children.push(child); this.firstChild = child; }},
        querySelector(sel) {{ return sel === 'svg' ? childSvg : null; }},
        getBoundingClientRect() {{ return {{ width: 680, height: 310, left: 0, top: 0 }}; }},
        addEventListener() {{}},
        children,
        childSvg
      }};
    }}

    const elements = {{
      'sdf-data': makeEl('sdf-data', {json.dumps(sdf_match.group(1))}),
      'dashboard-data': makeEl('dashboard-data', {json.dumps(json_match.group(1))}),
      'mol-3d-viewer': makeEl('mol-3d-viewer'),
      'physics-viewport': makeEl('physics-viewport'),
      'physics-tooltip': makeEl('physics-tooltip'),
      'physics-toggle-btn': makeEl('physics-toggle-btn'),
      'spin-toggle-btn': makeEl('spin-toggle-btn'),
      'style-btn-stick': makeEl('style-btn-stick'),
      'style-btn-ballstick': makeEl('style-btn-ballstick'),
      'style-btn-sphere': makeEl('style-btn-sphere'),
      'style-btn-surface': makeEl('style-btn-surface'),
      'toggle-fg-btn': makeEl('toggle-fg-btn'),
      'chem-svg-container': makeEl('chem-svg-container'),
      'radar-svg-container': makeEl('radar-svg-container'),
      'hist-svg-container': makeEl('hist-svg-container'),
      'activity-table-search': makeEl('activity-table-search'),
      'similarity-filter': Object.assign(makeEl('similarity-filter'), {{ value: '70' }}),
      'mwt-filter': Object.assign(makeEl('mwt-filter'), {{ value: '800' }}),
      'similarity-filter-val': makeEl('similarity-filter-val'),
      'mwt-filter-val': makeEl('mwt-filter-val'),
      'visible-analog-count': makeEl('visible-analog-count')
    }};
    elements['style-btn-stick'].setAttribute('data-style-btn', 'stick');
    elements['style-btn-ballstick'].setAttribute('data-style-btn', 'ballstick');
    elements['style-btn-sphere'].setAttribute('data-style-btn', 'sphere');
    elements['style-btn-surface'].setAttribute('data-style-btn', 'surface');

    const styleBtnList = [
      elements['style-btn-stick'],
      elements['style-btn-ballstick'],
      elements['style-btn-sphere'],
      elements['style-btn-surface']
    ];

    global.window = {{
      // Simulate sandboxed iframe where localStorage throws SecurityError
      get localStorage() {{
        throw new Error('SecurityError: Storage is disabled inside sandboxed iframe');
      }},
      matchMedia: () => ({{ matches: false }}),
      addEventListener() {{}},
      requestAnimationFrame: () => 1,
      cancelAnimationFrame() {{}}
    }};
    global.requestAnimationFrame = () => 1;
    global.cancelAnimationFrame = () => {{}};

    global.document = {{
      readyState: 'complete',
      documentElement: makeEl('html'),
      body: makeEl('body'),
      getElementById(id) {{ return elements[id] || null; }},
      createElementNS(ns, tag) {{ return makeEl(tag); }},
      createElement(tag) {{ return makeEl(tag); }},
      querySelectorAll(sel) {{ return sel === '[data-style-btn]' ? styleBtnList : []; }},
      addEventListener() {{}}
    }};

    {js_code}

    // Verify 1: Theme toggle works even when localStorage throws
    assert.strictEqual(document.documentElement.classList.contains('dark'), false);
    toggleTheme();
    assert.strictEqual(document.documentElement.classList.contains('dark'), true);
    assert.strictEqual(document.documentElement.getAttribute('data-theme'), 'dark');
    toggleTheme();
    assert.strictEqual(document.documentElement.classList.contains('dark'), false);
    assert.strictEqual(document.documentElement.getAttribute('data-theme'), 'light');

    // Verify 2: 3D Conformer ('Structural Data Diagram') rendered 21 atoms and bonds
    const svg3d = elements['mol-3d-viewer'].firstChild;
    assert(svg3d, '3D SVG must be appended to #mol-3d-viewer');
    assert(svg3d.innerHTML.includes('<circle'), '3D conformer must render atom circles');
    assert(svg3d.innerHTML.includes('<line'), '3D conformer must render bond lines');
    const circleCount3D = (svg3d.innerHTML.match(/<circle/g) || []).length;
    assert.strictEqual(circleCount3D, 21, 'Aspirin 3D conformer must render 21 atoms');

    // Verify initial active style button is 'ballstick' only
    assert.strictEqual(elements['style-btn-ballstick'].classList.contains('active'), true);
    assert.strictEqual(elements['style-btn-stick'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-sphere'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-surface'].classList.contains('active'), false);

    // Test 3D Zoom In, Zoom Out, Reset View, Style switches + Button Highlighting, and Spin
    zoom3D(1.25);
    assert.strictEqual(elements['mol-3d-viewer'].getAttribute('data-zoom'), '1.25');
    zoom3D(0.8);
    assert.strictEqual(elements['mol-3d-viewer'].getAttribute('data-zoom'), '1');

    setViewerStyle('stick');
    assert.strictEqual(elements['mol-3d-viewer'].getAttribute('data-style'), 'stick');
    assert.strictEqual(elements['style-btn-stick'].classList.contains('active'), true);
    assert.strictEqual(elements['style-btn-ballstick'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-sphere'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-surface'].classList.contains('active'), false);

    setViewerStyle('surface');
    assert.strictEqual(elements['mol-3d-viewer'].getAttribute('data-style'), 'surface');
    assert(svg3d.innerHTML.includes('fill-opacity="0.16"'), 'Surface style must render VDW surface');
    assert.strictEqual(elements['style-btn-surface'].classList.contains('active'), true);
    assert.strictEqual(elements['style-btn-stick'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-ballstick'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-sphere'].classList.contains('active'), false);

    setViewerStyle('sphere');
    assert.strictEqual(elements['mol-3d-viewer'].getAttribute('data-style'), 'sphere');
    assert.strictEqual(elements['style-btn-sphere'].classList.contains('active'), true);
    assert.strictEqual(elements['style-btn-surface'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-ballstick'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-stick'].classList.contains('active'), false);

    toggleSpin3D();
    assert.strictEqual(elements['spin-toggle-btn'].textContent, 'Pause Spin');
    reset3DView();
    assert.strictEqual(elements['mol-3d-viewer'].getAttribute('data-zoom'), '1');
    assert.strictEqual(elements['mol-3d-viewer'].getAttribute('data-style'), 'ballstick');
    assert.strictEqual(elements['spin-toggle-btn'].textContent, 'Auto-Spin');
    assert.strictEqual(elements['style-btn-ballstick'].classList.contains('active'), true);
    assert.strictEqual(elements['style-btn-stick'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-sphere'].classList.contains('active'), false);
    assert.strictEqual(elements['style-btn-surface'].classList.contains('active'), false);

    // Verify 3: Interactive Relationship Node Map rendered all nodes & edges
    const netHtml = elements['physics-viewport'].innerHTML;
    assert(netHtml.includes('id="physics-graph-group"'), 'Node map must render physics-graph-group');
    assert(netHtml.includes('CHEMBL25'), 'Node map must include query compound CHEMBL25');
    assert(netHtml.includes('CHEMBL26'), 'Node map must include analog CHEMBL26');
    assert(netHtml.includes('CHEMBL112'), 'Node map must include analog CHEMBL112');
    const nodeCircles = (netHtml.match(/<circle/g) || []).length;
    assert.strictEqual(nodeCircles, 5, 'Node map must render 5 compound nodes');

    // Test Node Map Zoom In, Zoom Out, Reset View
    zoomNetwork(1.25);
    assert.strictEqual(elements['physics-viewport'].getAttribute('data-zoom'), '1.25');
    assert(elements['physics-viewport'].innerHTML.includes('scale(1.250)'));
    zoomNetwork(0.8);
    assert.strictEqual(elements['physics-viewport'].getAttribute('data-zoom'), '1');
    resetPhysicsPositions();
    assert.strictEqual(elements['physics-viewport'].getAttribute('data-zoom'), '1');

    // Verify 4: 2D Structure, Radar, and Histogram Zoom In / Zoom Out / Reset View
    for (const target of ['chem2d', 'radar', 'hist']) {{
      const cid = viewContainers[target];
      zoomView(target, 1.25);
      assert.strictEqual(elements[cid].getAttribute('data-zoom'), '1.25');
      assert(elements[cid].childSvg.style.transform.includes('scale(1.250)'));
      zoomView(target, 0.8);
      assert.strictEqual(elements[cid].getAttribute('data-zoom'), '1');
      resetView(target);
      assert.strictEqual(elements[cid].getAttribute('data-zoom'), '1');
    }}

    clearInterval(spinTimer3D);
    console.log('ALL_RUNTIME_ASSERTIONS_PASSED');
    """

    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False, encoding="utf-8") as hf:
      hf.write(harness_js)
      harness_path = hf.name

    try:
      run_res = subprocess.run(["node", harness_path], capture_output=True, text=True, check=False)
      self.assertEqual(run_res.returncode, 0, f"Runtime verification failed:\nSTDOUT: {run_res.stdout}\nSTDERR: {run_res.stderr}")
      self.assertIn("ALL_RUNTIME_ASSERTIONS_PASSED", run_res.stdout)
    finally:
      if os.path.exists(harness_path):
        os.remove(harness_path)


if __name__ == "__main__":
  unittest.main()
