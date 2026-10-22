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
including an SVG radar chart and 2D/3D structure viewers strictly adhering to
DESIGN.md guidelines.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import os
import re
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
  """Evaluate Lipinski Rule of 5 and return metric values and violations."""
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
    badge_class = "status-normal"
  elif violations == 1:
    badge = "Moderate (1 violation)"
    badge_class = "status-warning"
  else:
    badge = f"High Violation ({violations} violations)"
    badge_class = "status-critical"

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


def generate_radar_svg(metrics: dict[str, Any], width: int = 360, height: int = 300) -> str:
  """Generate an SVG Radar chart comparing compound properties against Ro5."""
  cx, cy = width / 2.0, height / 2.0 - 8.0
  radius = 95.0
  keys = ["mwt", "alogp", "hbd", "hba", "rtb", "psa"]
  num_axes = len(keys)
  angles = [i * (2.0 * math.pi / num_axes) - (math.pi / 2.0) for i in range(num_axes)]

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
        f'<polygon points="{pts_closed}" fill="none" stroke="var(--chart-grid, #E7E5E4)" stroke-width="1" stroke-dasharray="2,2"/>'
    )

  threshold_frac = 0.65
  thresh_r = radius * threshold_frac
  thresh_pts = [
      f"{cx + thresh_r * math.cos(a):.1f},{cy + thresh_r * math.sin(a):.1f}"
      for a in angles
  ]
  thresh_polygon = " ".join(thresh_pts)

  mol_pts = []
  mol_dots = []
  axis_lines = []
  axis_labels = []

  has_violations = metrics.get("violations", 0) > 0
  poly_stroke = "#E11D48" if has_violations else "#0D9488"
  poly_fill = "rgba(225, 29, 72, 0.20)" if has_violations else "rgba(13, 148, 136, 0.22)"

  for idx, k in enumerate(keys):
    a = angles[idx]
    info = RO5_THRESHOLDS[k]
    val = float(metrics.get(k, 0.0))
    limit = info["max"]
    scale_max = info["scale_max"]

    if val <= limit:
      val_frac = (val / limit) * threshold_frac if limit > 0 else 0.0
    else:
      excess = min(val - limit, scale_max - limit)
      val_frac = threshold_frac + (excess / (scale_max - limit)) * (1.0 - threshold_frac)
    val_frac = max(0.05, min(1.05, val_frac))

    pt_x = cx + radius * val_frac * math.cos(a)
    pt_y = cy + radius * val_frac * math.sin(a)
    mol_pts.append(f"{pt_x:.1f},{pt_y:.1f}")
    dot_color = "#E11D48" if val > limit else "#0D9488"
    mol_dots.append(
        f'<circle cx="{pt_x:.1f}" cy="{pt_y:.1f}" r="4" fill="{dot_color}"><title>{info["name"]}: {val} (Limit: {limit})</title></circle>'
    )

    end_x = cx + radius * math.cos(a)
    end_y = cy + radius * math.sin(a)
    axis_lines.append(
        f'<line x1="{cx:.1f}" y1="{cy:.1f}" x2="{end_x:.1f}" y2="{end_y:.1f}" stroke="var(--chart-grid, #E7E5E4)" stroke-width="1"/>'
    )

    label_dist = radius + 22.0
    lbl_x = cx + label_dist * math.cos(a)
    lbl_y = cy + label_dist * math.sin(a)
    text_anchor = "middle"
    if math.cos(a) > 0.3:
      text_anchor = "start"
    elif math.cos(a) < -0.3:
      text_anchor = "end"
    axis_labels.append(
        f'<text x="{lbl_x:.1f}" y="{lbl_y + 4.0:.1f}" font-family="IBM Plex Mono, monospace" font-size="10" font-weight="500" fill="var(--text-secondary, #787774)" text-anchor="{text_anchor}">{info["name"]}</text>'
    )

  mol_polygon = " ".join(mol_pts)

  svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="100%">
  <rect width="{width}" height="{height}" fill="transparent"/>
  {"".join(grid_paths)}
  {"".join(axis_lines)}
  <polygon points="{thresh_polygon}" fill="rgba(5, 150, 105, 0.08)" stroke="#059669" stroke-width="1.6" stroke-dasharray="4,3"/>
  <polygon points="{mol_polygon}" fill="{poly_fill}" stroke="{poly_stroke}" stroke-width="2.2"/>
  {"".join(mol_dots)}
  {"".join(axis_labels)}
  <g transform="translate(14, {height - 16})">
    <line x1="0" y1="6" x2="16" y2="6" stroke="#059669" stroke-width="2" stroke-dasharray="3,2"/>
    <text x="22" y="9" font-family="IBM Plex Mono, monospace" font-size="10" fill="#059669">Ro5 Boundary</text>
    <rect x="125" y="2" width="12" height="8" fill="{poly_fill}" stroke="{poly_stroke}" stroke-width="1.5"/>
    <text x="143" y="9" font-family="IBM Plex Mono, monospace" font-size="10" fill="{poly_stroke}">Target Compound</text>
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
    svg = re.sub(r'<rect[^>]*fill=["\']#(?:ffffff|fff)["\'][^>]*/>', '', svg, flags=re.IGNORECASE)
    svg = re.sub(r'<rect[^>]*fill=["\']white["\'][^>]*/>', '', svg, flags=re.IGNORECASE)
    svg = re.sub(r'Roboto,?\s*', '', svg, flags=re.IGNORECASE)
    svg = re.sub(r'Inter,?\s*', '', svg, flags=re.IGNORECASE)
    return svg

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
    return f'<div class="empty-state"><span class="section-label">Canonical SMILES</span><code class="mono-code">{html.escape(smiles)}</code></div>'

  return '<div class="empty-state">2D Structure Preview Unavailable</div>'


def render_compound_card(
    molecule_data: dict[str, Any],
    sdf_content: str | None = None,
    svg_content: str | None = None,
) -> str:
  """Render a flat Bento grid section for the compound with Zoom/Reset controls."""
  mol = molecule_data
  props = mol.get("molecule_properties") or mol.get("properties") or {}
  metrics = evaluate_lipinski(props)
  radar_svg = generate_radar_svg(metrics)

  chembl_id = mol.get("molecule_chembl_id") or "UNKNOWN"
  pref_name = mol.get("pref_name") or "Unspecified Compound"
  max_phase = mol.get("max_phase")
  phase_text = f"Phase {max_phase}" if max_phase is not None else "Preclinical"
  formula = props.get("full_molformula") or "N/A"
  structures = mol.get("molecule_structures") or mol.get("structure") or {}
  smiles = structures.get("canonical_smiles") or "N/A"
  inchi_key = structures.get("standard_inchi_key") or "N/A"

  has_3d = bool(sdf_content and ("V2000" in sdf_content or "V3000" in sdf_content or "M  END" in sdf_content))
  svg_display = render_color_chemical_svg(svg_content, chembl_id=chembl_id, smiles=smiles)

  # Per DESIGN.md Section 1: Only display status indicator when in an error/violation state
  violation_indicator = ""
  if metrics["violations"] > 0:
    alert_class = "alert-critical" if metrics["violations"] > 1 else "alert-warning"
    violation_indicator = f'<span class="status-alert {alert_class}">{html.escape(metrics["badge"])}</span>'

  card_html = f"""
<section class="bento-panel" id="compound-overview-section">
  <div class="bento-header">
    <div class="header-title-group">
      <h2 class="panel-title">{html.escape(pref_name)}</h2>
      <span class="mono-meta">{html.escape(chembl_id)} &middot; {html.escape(formula)} &middot; {html.escape(phase_text)}</span>
    </div>
    <div class="header-actions">
      {violation_indicator}
    </div>
  </div>

  <div class="bento-grid-2">
    <!-- Left Cell: 2D Chemical Structure & 3D Conformer -->
    <div class="bento-cell bento-split-vertical">
      <div class="bento-subcell">
        <div class="subcell-toolbar">
          <h3 class="cell-title">2D Chemical Structure</h3>
          <div class="control-group">
            <button type="button" id="chem2d-zoom-in" onclick="zoomView('chem2d', 1.25)" class="btn-control" title="Zoom In">Zoom In</button>
            <button type="button" id="chem2d-zoom-out" onclick="zoomView('chem2d', 0.8)" class="btn-control" title="Zoom Out">Zoom Out</button>
            <button type="button" id="chem2d-reset" onclick="resetView('chem2d')" class="btn-control" title="Reset View">Reset View</button>
            <button type="button" id="toggle-fg-btn" onclick="toggleFunctionalGroupHighlights()" class="btn-control">Toggle Highlights</button>
          </div>
        </div>
        <div id="chem-svg-container" class="viz-canvas-2d interactive-viewport" data-zoom-target="chem2d">
          <div id="chem2d-stage" class="zoom-stage">
            {svg_display}
          </div>
        </div>
        <div class="viz-legend">
          <span class="legend-item"><span class="legend-dot" style="background:#E11D48"></span>Carbonyl / Hydroxyl O</span>
          <span class="legend-item"><span class="legend-dot" style="background:#D97706"></span>Ester O</span>
          <span class="legend-item"><span class="legend-dot" style="background:#4F46E5"></span>Methyl CH3</span>
          <span class="legend-item"><span class="legend-dot" style="background:#0D9488"></span>Aromatic Ring</span>
        </div>
      </div>

      <div class="bento-subcell">
        <div class="subcell-toolbar">
          <h3 class="cell-title">Structural Data Diagram</h3>
          <div class="control-group">
            <button type="button" id="mol3d-zoom-in" onclick="zoom3D(1.25)" class="btn-control" title="Zoom In">Zoom In</button>
            <button type="button" id="mol3d-zoom-out" onclick="zoom3D(0.8)" class="btn-control" title="Zoom Out">Zoom Out</button>
            <button type="button" id="mol3d-reset" onclick="reset3DView()" class="btn-control" title="Reset View">Reset View</button>
            <button type="button" id="spin-toggle-btn" onclick="toggleSpin3D()" class="btn-control">Auto-Spin</button>
            <button type="button" id="style-btn-stick" data-style-btn="stick" aria-pressed="false" onclick="setViewerStyle('stick')" class="btn-control">Sticks</button>
            <button type="button" id="style-btn-ballstick" data-style-btn="ballstick" aria-pressed="true" onclick="setViewerStyle('ballstick')" class="btn-control active">Ball &amp; Stick</button>
            <button type="button" id="style-btn-sphere" data-style-btn="sphere" aria-pressed="false" onclick="setViewerStyle('sphere')" class="btn-control">Spheres</button>
            <button type="button" id="style-btn-surface" data-style-btn="surface" aria-pressed="false" onclick="setViewerStyle('surface')" class="btn-control">Surface</button>
          </div>
        </div>
        <div id="mol-3d-viewer" class="viz-canvas-3d" data-has-3d="{str(has_3d).lower()}">
        </div>
        <div class="viz-legend">
          <span class="mono-meta">Drag to rotate 3D conformer &middot; Scroll or use buttons to zoom</span>
          <div class="viz-legend-inline">
            <span class="legend-item"><span class="legend-dot" style="background:#0D9488"></span>C</span>
            <span class="legend-item"><span class="legend-dot" style="background:#E11D48"></span>O</span>
            <span class="legend-item"><span class="legend-dot" style="background:#787774"></span>H</span>
          </div>
        </div>
      </div>
    </div>

    <!-- Right Cell: Lipinski Radar & Flat Property Matrix -->
    <div class="bento-cell bento-split-vertical">
      <div class="bento-subcell">
        <div class="subcell-toolbar">
          <h3 class="cell-title">Physicochemical Profile (Lipinski Ro5)</h3>
          <div class="control-group">
            <span class="mono-meta">{metrics["violations"]} Violations</span>
            <button type="button" id="radar-zoom-in" onclick="zoomView('radar', 1.25)" class="btn-control" title="Zoom In">Zoom In</button>
            <button type="button" id="radar-zoom-out" onclick="zoomView('radar', 0.8)" class="btn-control" title="Zoom Out">Zoom Out</button>
            <button type="button" id="radar-reset" onclick="resetView('radar')" class="btn-control" title="Reset View">Reset View</button>
          </div>
        </div>
        <div id="radar-svg-container" class="viz-canvas-radar interactive-viewport" data-zoom-target="radar">
          <div id="radar-stage" class="zoom-stage">
            {radar_svg}
          </div>
        </div>
      </div>

      <div class="bento-subcell no-padding">
        <div class="subcell-toolbar padded-toolbar">
          <h3 class="cell-title">Physicochemical Descriptors</h3>
        </div>
        <div class="descriptor-matrix">
          <div class="descriptor-cell">
            <span class="descriptor-label">Molecular Weight</span>
            <span class="descriptor-value" style="color: {'#E11D48' if metrics['mwt'] > 500 else 'var(--text-primary)'}">{metrics["mwt"]:.2f} <small>Da</small></span>
            <span class="descriptor-limit">Threshold &le; 500</span>
          </div>
          <div class="descriptor-cell">
            <span class="descriptor-label">Lipophilicity (AlogP)</span>
            <span class="descriptor-value" style="color: {'#E11D48' if metrics['alogp'] > 5.0 else 'var(--text-primary)'}">{metrics["alogp"]:.2f}</span>
            <span class="descriptor-limit">Threshold &le; 5.0</span>
          </div>
          <div class="descriptor-cell">
            <span class="descriptor-label">H-Bond Donors</span>
            <span class="descriptor-value" style="color: {'#E11D48' if metrics['hbd'] > 5 else 'var(--text-primary)'}">{metrics["hbd"]}</span>
            <span class="descriptor-limit">Threshold &le; 5</span>
          </div>
          <div class="descriptor-cell">
            <span class="descriptor-label">H-Bond Acceptors</span>
            <span class="descriptor-value" style="color: {'#E11D48' if metrics['hba'] > 10 else 'var(--text-primary)'}">{metrics["hba"]}</span>
            <span class="descriptor-limit">Threshold &le; 10</span>
          </div>
          <div class="descriptor-cell">
            <span class="descriptor-label">Rotatable Bonds</span>
            <span class="descriptor-value">{metrics["rtb"]}</span>
            <span class="descriptor-limit">Threshold &le; 10</span>
          </div>
          <div class="descriptor-cell">
            <span class="descriptor-label">Polar Surface Area</span>
            <span class="descriptor-value">{metrics["psa"]:.1f} <small>&Aring;&sup2;</small></span>
            <span class="descriptor-limit">Threshold &le; 140</span>
          </div>
        </div>
        <div class="identifier-strip">
          <div class="identifier-row">
            <span class="identifier-key">SMILES</span>
            <span class="identifier-val" title="{html.escape(smiles)}">{html.escape(smiles)}</span>
          </div>
          <div class="identifier-row">
            <span class="identifier-key">InChIKey</span>
            <span class="identifier-val">{html.escape(inchi_key)}</span>
          </div>
        </div>
      </div>
    </div>
  </div>
</section>
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
