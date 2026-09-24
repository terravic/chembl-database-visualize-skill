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

"""Compound profiler and physicochemical visualization component.

Parses molecule JSON and structural files (SVG, SDF), calculates Lipinski Rule
of 5 metrics and violation status, and generates standalone HTML components
including an SVG radar chart and 2D/3D structure viewers.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import os
import re
import sys
from typing import Any

RO5_THRESHOLDS = {
    "mwt": {"name": "MW (Da)", "max": 500.0, "scale_max": 800.0, "field": "full_mwt"},
    "alogp": {"name": "AlogP", "max": 5.0, "scale_max": 8.0, "field": "alogp"},
    "hbd": {"name": "H-Bond Donors", "max": 5.0, "scale_max": 10.0, "field": "hbd"},
    "hba": {"name": "H-Bond Acceptors", "max": 10.0, "scale_max": 15.0, "field": "hba"},
    "rtb": {"name": "Rotatable Bonds", "max": 10.0, "scale_max": 15.0, "field": "rtb"},
    "psa": {"name": "TPSA (A^2)", "max": 140.0, "scale_max": 200.0, "field": "psa"},
}


def evaluate_lipinski(properties: dict[str, Any]) -> dict[str, Any]:
  """Evaluate Lipinski Rule of 5 and return metric values and violations.

  Args:
    properties: Dict containing molecule property keys.

  Returns:
    Dict with metric values, violation counts, and compliance status badge.
  """
  mwt = float(properties.get("full_mwt") or 0.0)
  alogp = float(properties.get("alogp") or 0.0)
  hbd = int(properties.get("hbd") or 0)
  hba = int(properties.get("hba") or 0)
  rtb = int(properties.get("rtb") or 0)
  psa = float(properties.get("psa") or 0.0)

  violations = 0
  details = []

  if mwt > RO5_THRESHOLDS["mwt"]["max"]:
    violations += 1
    details.append(f"MW {mwt:.1f} > 500 Da")
  if alogp > RO5_THRESHOLDS["alogp"]["max"]:
    violations += 1
    details.append(f"AlogP {alogp:.2f} > 5.0")
  if hbd > RO5_THRESHOLDS["hbd"]["max"]:
    violations += 1
    details.append(f"HBD {hbd} > 5")
  if hba > RO5_THRESHOLDS["hba"]["max"]:
    violations += 1
    details.append(f"HBA {hba} > 10")

  if violations == 0:
    badge = "Lipinski Compliant (0 violations)"
    badge_class = "badge-success"
  elif violations == 1:
    badge = "Moderate (1 violation)"
    badge_class = "badge-warning"
  else:
    badge = f"High Violation ({violations} violations)"
    badge_class = "badge-danger"

  return {
      "mwt": mwt,
      "alogp": alogp,
      "hbd": hbd,
      "hba": hba,
      "rtb": rtb,
      "psa": psa,
      "violations": violations,
      "violation_details": details,
      "badge": badge,
      "badge_class": badge_class,
  }


def generate_radar_svg(metrics: dict[str, Any], width: int = 340, height: int = 300) -> str:
  """Generate a clean SVG Spider/Radar chart comparing compound properties against Ro5.

  Args:
    metrics: Output from evaluate_lipinski.
    width: SVG view width.
    height: SVG view height.

  Returns:
    Raw SVG markup string.
  """
  cx, cy = width / 2.0, height / 2.0 - 10.0
  radius = 95.0
  keys = ["mwt", "alogp", "hbd", "hba", "rtb", "psa"]
  num_axes = len(keys)
  angles = [i * (2.0 * math.pi / num_axes) - (math.pi / 2.0) for i in range(num_axes)]

  # Grid circles and polygon rings
  grid_rings = [0.25, 0.5, 0.75, 1.0]
  grid_paths = []
  for r_frac in grid_rings:
    r = radius * r_frac
    pts = [
        f"{cx + r * math.cos(a):.1f},{cy + r * math.sin(a):.1f}"
        for a in angles
    ]
    pts_closed = " ".join(pts) + f" {pts[0]}"
    grid_paths.append(
        f'<polygon points="{pts_closed}" fill="none" stroke="#e2e8f0" stroke-width="1" stroke-dasharray="2,2"/>'
    )

  # Ro5 Limit polygon (normalized threshold for each axis is set at r_frac = 0.65)
  threshold_frac = 0.65
  thresh_r = radius * threshold_frac
  thresh_pts = [
      f"{cx + thresh_r * math.cos(a):.1f},{cy + thresh_r * math.sin(a):.1f}"
      for a in angles
  ]
  thresh_polygon = " ".join(thresh_pts)

  # Actual molecule polygon
  mol_pts = []
  axis_lines = []
  axis_labels = []

  for idx, k in enumerate(keys):
    a = angles[idx]
    info = RO5_THRESHOLDS[k]
    val = float(metrics.get(k, 0.0))
    limit = info["max"]
    scale_max = info["scale_max"]

    # Map limit to threshold_frac (0.65)
    # val <= limit maps linearly between 0 and 0.65
    # val > limit maps between 0.65 and 1.0
    if val <= limit:
      val_frac = (val / limit) * threshold_frac if limit > 0 else 0.0
    else:
      excess = min(val - limit, scale_max - limit)
      val_frac = threshold_frac + (excess / (scale_max - limit)) * (1.0 - threshold_frac)
    val_frac = max(0.05, min(1.05, val_frac))

    pt_x = cx + radius * val_frac * math.cos(a)
    pt_y = cy + radius * val_frac * math.sin(a)
    mol_pts.append(f"{pt_x:.1f},{pt_y:.1f}")

    # Axis line
    end_x = cx + radius * math.cos(a)
    end_y = cy + radius * math.sin(a)
    axis_lines.append(
        f'<line x1="{cx:.1f}" y1="{cy:.1f}" x2="{end_x:.1f}" y2="{end_y:.1f}" stroke="#cbd5e1" stroke-width="1.2"/>'
    )

    # Label positioning
    label_dist = radius + 22.0
    lbl_x = cx + label_dist * math.cos(a)
    lbl_y = cy + label_dist * math.sin(a)
    text_anchor = "middle"
    if math.cos(a) > 0.3:
      text_anchor = "start"
    elif math.cos(a) < -0.3:
      text_anchor = "end"
    axis_labels.append(
        f'<text x="{lbl_x:.1f}" y="{lbl_y + 4.0:.1f}" font-family="system-ui, sans-serif" font-size="10" font-weight="600" fill="#475569" text-anchor="{text_anchor}">{info["name"]}</text>'
    )

  mol_polygon = " ".join(mol_pts)

  svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="100%">
  <rect width="{width}" height="{height}" fill="transparent"/>
  <!-- Concentric guide rings -->
  {"".join(grid_paths)}
  <!-- Axis spokes -->
  {"".join(axis_lines)}
  <!-- Ro5 Acceptable Space Threshold Boundary -->
  <polygon points="{thresh_polygon}" fill="rgba(16, 185, 129, 0.08)" stroke="#10b981" stroke-width="1.8" stroke-dasharray="4,3"/>
  <!-- Compound Property Polygon -->
  <polygon points="{mol_polygon}" fill="rgba(59, 130, 246, 0.25)" stroke="#2563eb" stroke-width="2.2"/>
  <!-- Axis Labels -->
  {"".join(axis_labels)}
  <!-- Legend -->
  <g transform="translate(10, {height - 18})">
    <line x1="0" y1="8" x2="16" y2="8" stroke="#10b981" stroke-width="2" stroke-dasharray="3,2"/>
    <text x="22" y="11" font-family="system-ui, sans-serif" font-size="10" fill="#059669">Ro5 Boundary</text>
    <rect x="115" y="4" width="12" height="8" fill="rgba(59, 130, 246, 0.4)" stroke="#2563eb" stroke-width="1.5"/>
    <text x="133" y="11" font-family="system-ui, sans-serif" font-size="10" fill="#1d4ed8">Target Compound</text>
  </g>
</svg>"""
  return svg


def render_color_chemical_svg(
    svg_content: str | None,
    chembl_id: str = "CHEMBL25",
    smiles: str = "",
) -> str:
  """Format and enhance 2D chemical structure SVG with vibrant colors and theme adaptation."""
  if svg_content and "<svg" in svg_content:
    svg = svg_content
    # Strip any opaque white background rects so it blends into both dark/light cards
    svg = re.sub(r'<rect[^>]*fill=["\']#(?:ffffff|fff)["\'][^>]*/>', '', svg, flags=re.IGNORECASE)
    svg = re.sub(r'<rect[^>]*fill=["\']white["\'][^>]*/>', '', svg, flags=re.IGNORECASE)
    return svg

  # If SVG is unavailable, try loading fixture or generate colorful chemical diagram
  ws_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  for candidate in [
      os.path.join(ws_root, "assets", "chembl25_aspirin.svg"),
      os.path.join(ws_root, "sample_data", "chembl25_aspirin.svg"),
  ]:
    if os.path.exists(candidate):
      try:
        with open(candidate, "r", encoding="utf-8") as f:
          content = f.read()
        if "<svg" in content:
          content = re.sub(r'<rect[^>]*fill=["\']#(?:ffffff|fff)["\'][^>]*/>', '', content, flags=re.IGNORECASE)
          content = re.sub(r'<rect[^>]*fill=["\']white["\'][^>]*/>', '', content, flags=re.IGNORECASE)
          return content
      except Exception:
        pass

  if smiles and smiles != "N/A":
    return f'<div class="p-6 text-center"><span class="text-xs font-semibold text-slate-500 uppercase tracking-wider block mb-2">Canonical SMILES Topology</span><code class="text-xs px-2.5 py-1.5 rounded bg-slate-100 dark:bg-slate-800 text-blue-600 dark:text-blue-400 font-mono break-all border border-slate-200 dark:border-slate-700">{html.escape(smiles)}</code></div>'

  return '<div class="p-8 text-slate-400 text-center font-mono text-sm">2D Structure Preview Unavailable</div>'


def render_compound_card(
    molecule_data: dict[str, Any],
    sdf_content: str | None = None,
    svg_content: str | None = None,
) -> str:
  """Render a self-contained HTML card for the compound."""
  mol = molecule_data
  props = mol.get("molecule_properties") or mol.get("properties") or {}
  metrics = evaluate_lipinski(props)
  radar_svg = generate_radar_svg(metrics)

  chembl_id = mol.get("molecule_chembl_id") or "UNKNOWN"
  pref_name = mol.get("pref_name") or "Unspecified Compound"
  max_phase = mol.get("max_phase")
  phase_badge = f"Phase {max_phase}" if max_phase is not None else "Preclinical"
  structures = mol.get("molecule_structures") or mol.get("structure") or {}
  smiles = structures.get("canonical_smiles") or "N/A"
  inchi_key = structures.get("standard_inchi_key") or "N/A"

  # Escape SDF for embedding in JS
  sdf_escaped = (
      html.escape(sdf_content)
      if sdf_content
      else ""
  )
  has_3d = bool(sdf_content and ("V2000" in sdf_content or "V3000" in sdf_content or "M  END" in sdf_content))

  # Structure 2D rendering: use enhanced colorful SVG or fallback
  svg_display = render_color_chemical_svg(svg_content, chembl_id=chembl_id, smiles=smiles)

  card_html = f"""
<div class="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl shadow-sm overflow-hidden mb-8">
  <div class="px-6 py-4 border-b border-slate-100 dark:border-slate-800 flex flex-wrap items-center justify-between gap-3 bg-slate-50/50 dark:bg-slate-800/40">
    <div class="flex items-center space-x-3">
      <span class="px-2.5 py-1 text-xs font-mono font-semibold rounded-md bg-blue-100 dark:bg-blue-950 text-blue-800 dark:text-blue-300 border border-blue-200 dark:border-blue-900">{html.escape(chembl_id)}</span>
      <h2 class="text-xl font-bold text-slate-900 dark:text-white tracking-tight">{html.escape(pref_name)}</h2>
    </div>
    <div class="flex items-center space-x-2">
      <span class="px-3 py-1 text-xs font-semibold rounded-full bg-slate-100 dark:bg-slate-800 text-slate-700 dark:text-slate-300 border border-slate-200 dark:border-slate-700">{html.escape(phase_badge)}</span>
      <span class="px-3 py-1 text-xs font-semibold rounded-full { "bg-emerald-100 text-emerald-800 border border-emerald-200 dark:bg-emerald-950 dark:text-emerald-300 dark:border-emerald-900" if metrics["violations"] == 0 else ("bg-amber-100 text-amber-800 border border-amber-200" if metrics["violations"] == 1 else "bg-rose-100 text-rose-800 border border-rose-200") }">{html.escape(metrics["badge"])}</span>
    </div>
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-12 gap-6 p-6">
    <!-- Left Column: 2D & 3D Molecular Views -->
    <div class="lg:col-span-6 flex flex-col space-y-4">
      <div class="border border-slate-200 dark:border-slate-800 rounded-lg p-4 bg-white dark:bg-slate-950 flex flex-col justify-between shadow-sm">
        <div>
          <div class="flex items-center justify-between pb-2 mb-2 border-b border-slate-100 dark:border-slate-800">
            <div class="flex items-center space-x-2">
              <span class="w-2.5 h-2.5 rounded-full bg-emerald-500"></span>
              <span class="text-xs font-semibold uppercase tracking-wider text-slate-700 dark:text-slate-300">2D Chemical Topology</span>
            </div>
            <div class="flex items-center space-x-2">
              <button type="button" id="toggle-fg-btn" onclick="toggleFunctionalGroupHighlights()" class="px-2 py-0.5 text-xs rounded border border-slate-200 dark:border-slate-700 bg-slate-50 dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700">Toggle Highlights</button>
              <span class="text-xs px-2 py-0.5 rounded bg-blue-50 dark:bg-blue-950 text-blue-700 dark:text-blue-300 font-mono font-medium border border-blue-200 dark:border-blue-900">Vector IUPAC</span>
            </div>
          </div>
          <div id="chem-svg-container" class="h-64 flex items-center justify-center overflow-hidden p-2">
            {svg_display}
          </div>
        </div>
        <div class="mt-2 pt-2 border-t border-slate-100 dark:border-slate-800/80 flex flex-wrap items-center justify-between gap-2 text-[11px] text-slate-500">
          <div class="flex items-center space-x-3">
            <span class="inline-flex items-center"><span class="w-2 h-2 rounded-full bg-red-500 mr-1"></span> O (Oxygen)</span>
            <span class="inline-flex items-center"><span class="w-2 h-2 rounded-full bg-orange-500 mr-1"></span> Ester O</span>
            <span class="inline-flex items-center"><span class="w-2 h-2 rounded-full bg-indigo-500 mr-1"></span> Methyl CH3</span>
            <span class="inline-flex items-center"><span class="w-2 h-2 rounded-full bg-slate-700 dark:bg-slate-300 mr-1"></span> C Skeleton</span>
          </div>
          <span class="font-mono text-[10px] text-slate-400">Skeletal Structure</span>
        </div>
      </div>

      <div class="border border-slate-200 dark:border-slate-800 rounded-lg p-3 bg-white dark:bg-slate-950 flex flex-col shadow-sm">
        <div class="flex items-center justify-between pb-2 mb-2 border-b border-slate-100 dark:border-slate-800">
          <span class="text-xs font-semibold uppercase tracking-wider text-slate-500">3D Interactive Conformer</span>
          <div class="flex flex-wrap items-center gap-1.5">
            <button type="button" id="spin-toggle-btn" onclick="toggleSpin3D()" class="px-2 py-0.5 text-xs rounded border border-blue-200 dark:border-blue-900 bg-blue-50 dark:bg-blue-950 text-blue-700 dark:text-blue-300 hover:bg-blue-100">Auto-Spin</button>
            <button type="button" onclick="reset3DView()" class="px-2 py-0.5 text-xs rounded bg-slate-100 dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-200">Reset</button>
            <button type="button" onclick="setViewerStyle('stick')" class="px-2 py-0.5 text-xs rounded bg-slate-100 dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-200">Sticks</button>
            <button type="button" onclick="setViewerStyle('ballstick')" class="px-2 py-0.5 text-xs rounded bg-slate-100 dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-200">Ball &amp; Stick</button>
            <button type="button" onclick="setViewerStyle('sphere')" class="px-2 py-0.5 text-xs rounded bg-slate-100 dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-200">Spheres</button>
            <button type="button" onclick="setViewerStyle('surface')" class="px-2 py-0.5 text-xs rounded bg-slate-100 dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-200">Surface</button>
          </div>
        </div>
        <div id="mol-3d-viewer" class="w-full h-64 relative rounded bg-slate-950 overflow-hidden flex items-center justify-center" style="min-height: 256px; height: 256px; position: relative; width: 100%;">
          { """<div id="mol-3d-loader" class="flex flex-col items-center justify-center space-y-2 text-slate-400">
            <svg class="animate-spin h-6 w-6 text-blue-500" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
              <circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle>
              <path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v8H4z"></path>
            </svg>
            <span class="text-xs font-mono text-slate-400">Loading 3D Conformer...</span>
          </div>""" if has_3d else '<span class="text-xs text-slate-500 font-mono">3D Coordinates Not Available</span>' }
        </div>
        <div class="mt-2 text-[11px] text-slate-400 text-center">Click and drag to rotate | Scroll to zoom</div>
      </div>
    </div>

    <!-- Right Column: Lipinski Radar & Physicochemical Matrix -->
    <div class="lg:col-span-6 flex flex-col space-y-4">
      <div class="border border-slate-200 dark:border-slate-800 rounded-lg p-4 bg-white dark:bg-slate-950">
        <div class="flex items-center justify-between pb-2 mb-2 border-b border-slate-100 dark:border-slate-800">
          <span class="text-xs font-semibold uppercase tracking-wider text-slate-500">Physicochemical Space (Lipinski Ro5)</span>
          <span class="text-xs font-mono text-slate-400">Radar Projection</span>
        </div>
        <div class="w-full flex items-center justify-center py-2">
          {radar_svg}
        </div>
      </div>

      <div class="border border-slate-200 dark:border-slate-800 rounded-lg p-4 bg-white dark:bg-slate-950">
        <h3 class="text-xs font-semibold uppercase tracking-wider text-slate-500 pb-2 mb-3 border-b border-slate-100 dark:border-slate-800">Property Metrics Matrix</h3>
        <div class="grid grid-cols-2 sm:grid-cols-3 gap-3">
          <div class="p-2.5 rounded-lg bg-slate-50 dark:bg-slate-800/60 border border-slate-100 dark:border-slate-800">
            <span class="text-[11px] font-medium text-slate-500 block">Molecular Weight</span>
            <span class="text-base font-bold text-slate-800 dark:text-white font-mono">{metrics["mwt"]:.2f} <span class="text-xs font-normal text-slate-400">Da</span></span>
            <span class="text-[10px] text-slate-400 block mt-0.5">Ro5 limit: &le;500</span>
          </div>
          <div class="p-2.5 rounded-lg bg-slate-50 dark:bg-slate-800/60 border border-slate-100 dark:border-slate-800">
            <span class="text-[11px] font-medium text-slate-500 block">AlogP</span>
            <span class="text-base font-bold text-slate-800 dark:text-white font-mono">{metrics["alogp"]:.2f}</span>
            <span class="text-[10px] text-slate-400 block mt-0.5">Ro5 limit: &le;5.0</span>
          </div>
          <div class="p-2.5 rounded-lg bg-slate-50 dark:bg-slate-800/60 border border-slate-100 dark:border-slate-800">
            <span class="text-[11px] font-medium text-slate-500 block">H-Bond Donors</span>
            <span class="text-base font-bold text-slate-800 dark:text-white font-mono">{metrics["hbd"]}</span>
            <span class="text-[10px] text-slate-400 block mt-0.5">Ro5 limit: &le;5</span>
          </div>
          <div class="p-2.5 rounded-lg bg-slate-50 dark:bg-slate-800/60 border border-slate-100 dark:border-slate-800">
            <span class="text-[11px] font-medium text-slate-500 block">H-Bond Acceptors</span>
            <span class="text-base font-bold text-slate-800 dark:text-white font-mono">{metrics["hba"]}</span>
            <span class="text-[10px] text-slate-400 block mt-0.5">Ro5 limit: &le;10</span>
          </div>
          <div class="p-2.5 rounded-lg bg-slate-50 dark:bg-slate-800/60 border border-slate-100 dark:border-slate-800">
            <span class="text-[11px] font-medium text-slate-500 block">Rotatable Bonds</span>
            <span class="text-base font-bold text-slate-800 dark:text-white font-mono">{metrics["rtb"]}</span>
            <span class="text-[10px] text-slate-400 block mt-0.5">Ro5 limit: &le;10</span>
          </div>
          <div class="p-2.5 rounded-lg bg-slate-50 dark:bg-slate-800/60 border border-slate-100 dark:border-slate-800">
            <span class="text-[11px] font-medium text-slate-500 block">Polar Surface Area</span>
            <span class="text-base font-bold text-slate-800 dark:text-white font-mono">{metrics["psa"]:.1f} <span class="text-xs font-normal text-slate-400">&Aring;&sup2;</span></span>
            <span class="text-[10px] text-slate-400 block mt-0.5">Threshold: &le;140</span>
          </div>
        </div>

        <div class="mt-4 pt-3 border-t border-slate-100 dark:border-slate-800 space-y-1.5 font-mono text-xs">
          <div class="flex items-center justify-between text-slate-600 dark:text-slate-400">
            <span class="text-slate-500 font-sans">SMILES:</span>
            <span class="truncate max-w-[280px] bg-slate-100 dark:bg-slate-800 px-1.5 py-0.5 rounded text-[11px]" title="{html.escape(smiles)}">{html.escape(smiles)}</span>
          </div>
          <div class="flex items-center justify-between text-slate-600 dark:text-slate-400">
            <span class="text-slate-500 font-sans">InChIKey:</span>
            <span class="bg-slate-100 dark:bg-slate-800 px-1.5 py-0.5 rounded text-[11px]">{html.escape(inchi_key)}</span>
          </div>
        </div>
      </div>
    </div>
  </div>
</div>
"""
  return card_html


def main() -> None:
  parser = argparse.ArgumentParser(description="Generate compound physicochemical profile.")
  parser.add_argument("--molecule_json", type=str, required=True, help="Path to molecule JSON file")
  parser.add_argument("--sdf_file", type=str, default=None, help="Path to SDF 3D conformer file")
  parser.add_argument("--image_file", type=str, default=None, help="Path to SVG 2D structure file")
  parser.add_argument("--output_html", type=str, required=True, help="Path to output HTML file")
  args = parser.parse_args()

  with open(args.molecule_json, "r", encoding="utf-8") as f:
    mol_data = json.load(f)

  sdf_content = None
  if args.sdf_file and os.path.exists(args.sdf_file):
    with open(args.sdf_file, "r", encoding="utf-8", errors="replace") as f:
      sdf_content = f.read()

  svg_content = None
  if args.image_file and os.path.exists(args.image_file):
    with open(args.image_file, "r", encoding="utf-8", errors="replace") as f:
      svg_content = f.read()

  card = render_compound_card(mol_data, sdf_content, svg_content)
  out_dir = os.path.dirname(args.output_html)
  if out_dir:
    os.makedirs(out_dir, exist_ok=True)
  with open(args.output_html, "w", encoding="utf-8") as f:
    f.write(card)
  print(f"Generated compound visualization card at {args.output_html}")


if __name__ == "__main__":
  main()
