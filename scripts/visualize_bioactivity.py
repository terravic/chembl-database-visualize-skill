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

"""SAR and bioactivity distribution visualizer component.

Transforms normalized bioactivity records into interactive potency histograms,
target affinity comparisons, and filterable monospace data tables strictly
adhering to DESIGN.md guidelines.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import os
from typing import Any


def calculate_pic50(normalized_val_nm: float | None) -> float | None:
  """Calculate pIC50 = 9 - log10(nM)."""
  if normalized_val_nm is None:
    return None
  try:
    val = float(normalized_val_nm)
    if val <= 0:
      return None
    return round(9.0 - math.log10(val), 2)
  except (ValueError, TypeError):
    return None


def process_activities(raw_activities: list[dict[str, Any]]) -> dict[str, Any]:
  """Process activity records, compute pIC50, and calculate distributions."""
  processed = []
  pic50_values = []
  potency_counts = {"high": 0, "moderate": 0, "weak": 0}
  target_affinities: dict[str, list[float]] = {}

  for rec in raw_activities:
    std_val = rec.get("normalized_value_nM")
    if std_val is None:
      try:
        raw_std = rec.get("standard_value")
        raw_units = (rec.get("standard_units") or "").strip().lower()
        if raw_std is not None and raw_units in ("nm", ""):
          std_val = float(raw_std)
      except (ValueError, TypeError):
        std_val = None

    pic50 = rec.get("pchembl_value")
    if pic50 is not None:
      try:
        pic50 = float(pic50)
      except (ValueError, TypeError):
        pic50 = None
    if pic50 is None:
      pic50 = calculate_pic50(std_val)

    if std_val is not None:
      if std_val < 100.0:
        tier = "High (<100 nM)"
        potency_counts["high"] += 1
      elif std_val <= 10000.0:
        tier = "Moderate (100-10,000 nM)"
        potency_counts["moderate"] += 1
      else:
        tier = "Low (>10,000 nM)"
        potency_counts["weak"] += 1
    else:
      tier = "Unclassified"

    item = {
        "activity_id": rec.get("activity_id") or "N/A",
        "assay_chembl_id": rec.get("assay_chembl_id") or "N/A",
        "target_chembl_id": rec.get("target_chembl_id") or "N/A",
        "target_pref_name": rec.get("target_pref_name") or rec.get("target_name") or "Unspecified Target",
        "target_organism": rec.get("target_organism") or "Unknown",
        "standard_type": rec.get("standard_type") or "N/A",
        "standard_relation": rec.get("standard_relation") or "=",
        "normalized_value_nM": std_val,
        "pIC50": pic50,
        "potency_tier": tier,
        "is_inactive_anomaly": bool(std_val is not None and std_val > 10000.0),
        "document_chembl_id": rec.get("document_chembl_id") or "N/A",
    }
    processed.append(item)

    if pic50 is not None:
      pic50_values.append(pic50)
      t_name = item["target_pref_name"]
      target_affinities.setdefault(t_name, []).append(pic50)

  if pic50_values:
    pic50_values.sort()
    n = len(pic50_values)
    median_pic50 = (
        pic50_values[n // 2]
        if n % 2 != 0
        else (pic50_values[n // 2 - 1] + pic50_values[n // 2]) / 2.0
    )
    stats = {
        "count": n,
        "min": min(pic50_values),
        "max": max(pic50_values),
        "median": round(median_pic50, 2),
    }
  else:
    stats = {"count": 0, "min": 0.0, "max": 0.0, "median": 0.0}

  return {
      "records": processed,
      "stats": stats,
      "potency_counts": potency_counts,
      "pic50_values": pic50_values,
      "target_affinities": target_affinities,
  }


def generate_histogram_svg(pic50_values: list[float], width: int = 420, height: int = 230) -> str:
  """Generate an SVG histogram of pIC50 values using the DESIGN.md palette."""
  if not pic50_values:
    return f"""<svg viewBox="0 0 {width} {height}" width="100%" height="100%">
      <rect width="{width}" height="{height}" fill="transparent"/>
      <text x="{width/2}" y="{height/2}" text-anchor="middle" font-family="IBM Plex Mono, monospace" font-size="11" fill="var(--text-secondary, #787774)">No numerical pIC50 records available</text>
    </svg>"""

  bin_min, bin_max = 3.0, 11.0
  num_bins = 8
  bin_width = (bin_max - bin_min) / num_bins
  bins = [0] * num_bins

  for v in pic50_values:
    idx = int((v - bin_min) / bin_width)
    if idx < 0:
      idx = 0
    elif idx >= num_bins:
      idx = num_bins - 1
    bins[idx] += 1

  max_count = max(bins) if max(bins) > 0 else 1
  margin_left = 36
  margin_bottom = 34
  margin_top = 22
  margin_right = 16

  plot_w = width - margin_left - margin_right
  plot_h = height - margin_top - margin_bottom
  bar_w = plot_w / num_bins

  grid_lines = []
  for frac in (0.25, 0.5, 0.75, 1.0):
    gy = margin_top + plot_h * (1.0 - frac)
    grid_lines.append(
        f'<line x1="{margin_left}" y1="{gy:.1f}" x2="{width - margin_right}" y2="{gy:.1f}" stroke="var(--chart-grid, #E7E5E4)" stroke-width="1" stroke-dasharray="2,2"/>'
    )

  svg_bars = []
  svg_labels = []

  for i, count in enumerate(bins):
    x = margin_left + i * bar_w + 4
    bw = max(2, bar_w - 8)
    bh = (count / max_count) * (plot_h - 10)
    y = margin_top + plot_h - bh

    mid_bin = bin_min + (i + 0.5) * bin_width
    if mid_bin >= 8.0:
      fill_color = "#059669"
    elif mid_bin >= 7.0:
      fill_color = "#0D9488"
    elif mid_bin >= 5.0:
      fill_color = "#D97706"
    else:
      fill_color = "#E11D48"

    svg_bars.append(
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw:.1f}" height="{bh:.1f}" fill="{fill_color}" rx="3"/>'
    )
    if count > 0:
      svg_bars.append(
          f'<text x="{x + bw/2:.1f}" y="{y - 5:.1f}" font-family="IBM Plex Mono, monospace" font-size="10" font-weight="600" fill="var(--text-primary, #1C1917)" text-anchor="middle">{count}</text>'
      )

    lbl_x = margin_left + i * bar_w + bar_w / 2
    lbl_text = f"{bin_min + i * bin_width:.0f}"
    svg_labels.append(
        f'<text x="{lbl_x:.1f}" y="{height - 16}" font-family="IBM Plex Mono, monospace" font-size="10" fill="var(--text-secondary, #787774)" text-anchor="middle">{lbl_text}</text>'
    )

  svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="100%">
  <rect width="{width}" height="{height}" fill="transparent"/>
  {"".join(grid_lines)}
  <line x1="{margin_left}" y1="{margin_top + plot_h}" x2="{width - margin_right}" y2="{margin_top + plot_h}" stroke="var(--chart-grid, #E7E5E4)" stroke-width="1.5"/>
  <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_h}" stroke="var(--chart-grid, #E7E5E4)" stroke-width="1.5"/>
  {"".join(svg_bars)}
  {"".join(svg_labels)}
  <text x="{width/2}" y="{height - 2}" font-family="IBM Plex Mono, monospace" font-size="10" fill="var(--text-secondary, #787774)" text-anchor="middle">-log10 Affinity (pIC50 / pKi)</text>
</svg>"""
  return svg


def render_bioactivity_card(activity_data: dict[str, Any]) -> str:
  """Render flat Bento bioactivity distribution and monospace table section."""
  raw_list = activity_data.get("activities", [])
  processed = process_activities(raw_list)
  records = processed["records"]
  stats = processed["stats"]
  counts = processed["potency_counts"]
  hist_svg = generate_histogram_svg(processed["pic50_values"])

  table_rows = []
  for r in records[:50]:
    val_str = f"{r['normalized_value_nM']:.2f}" if r["normalized_value_nM"] is not None else "N/A"
    pic50_str = f"{r['pIC50']:.2f}" if r["pIC50"] is not None else "-"
    if r["is_inactive_anomaly"]:
      anomaly_cell = '<span class="status-alert alert-critical">&gt;10 &mu;M Threshold</span>'
      val_color = "color:#E11D48;font-weight:600;"
    else:
      anomaly_cell = f'<span class="mono-secondary">{html.escape(r["assay_chembl_id"])}</span>'
      val_color = "color:#059669;font-weight:600;" if (r["pIC50"] and r["pIC50"] >= 7.0) else "color:#D97706;"

    table_rows.append(
        f"""<tr class="table-row">
          <td class="cell-mono">{html.escape(str(r['activity_id']))}</td>
          <td class="cell-mono cell-truncate" title="{html.escape(r['target_pref_name'])}">{html.escape(r['target_pref_name'])}</td>
          <td class="cell-mono mono-secondary">{html.escape(r['target_organism'])}</td>
          <td class="cell-mono" style="color:#4F46E5;font-weight:600;">{html.escape(r['standard_type'])}</td>
          <td class="cell-mono">{html.escape(r['standard_relation'])} {val_str}</td>
          <td class="cell-mono" style="{val_color}">{pic50_str}</td>
          <td class="cell-mono">{anomaly_cell}</td>
        </tr>"""
    )

  rows_html = (
      "".join(table_rows)
      if table_rows
      else '<tr><td colspan="7" class="empty-table-cell">No bioactivity measurements found</td></tr>'
  )

  card_html = f"""
<section class="bento-panel" id="bioactivity-section">
  <div class="bento-header">
    <div class="header-title-group">
      <h2 class="panel-title">Bioactivity &amp; Potency Distribution</h2>
      <span class="mono-meta">Concentration normalized to nM</span>
    </div>
    <div class="control-group">
      <input type="text" id="activity-table-search" placeholder="Filter target or assay..." oninput="filterActivityTable()" class="input-filter"/>
      <button type="button" onclick="resetActivityFilter()" class="btn-control">Reset</button>
    </div>
  </div>

  <!-- Flat Metric Summary Strip -->
  <div class="metric-strip-4">
    <div class="metric-cell">
      <span class="descriptor-label">Assays Screened</span>
      <span class="metric-big" style="color:#0D9488">{stats['count']}</span>
    </div>
    <div class="metric-cell">
      <span class="descriptor-label">Median pIC50</span>
      <span class="metric-big">{stats['median']}</span>
    </div>
    <div class="metric-cell">
      <span class="descriptor-label">Peak Potency (pIC50)</span>
      <span class="metric-big" style="color:#059669">{stats['max']}</span>
    </div>
    <div class="metric-cell">
      <span class="descriptor-label">Sub-100 nM Hits</span>
      <span class="metric-big" style="color:#4F46E5">{counts['high']}</span>
    </div>
  </div>

  <div class="bento-grid-5-7">
    <!-- Left Cell: Potency Histogram -->
    <div class="bento-cell">
      <div class="subcell-toolbar padded-toolbar">
        <h3 class="cell-title">Potency Spread (pIC50)</h3>
        <div class="control-group">
          <button type="button" id="hist-zoom-in" onclick="zoomView('hist', 1.25)" class="btn-control" title="Zoom In">Zoom In</button>
          <button type="button" id="hist-zoom-out" onclick="zoomView('hist', 0.8)" class="btn-control" title="Zoom Out">Zoom Out</button>
          <button type="button" id="hist-reset" onclick="resetView('hist')" class="btn-control" title="Reset View">Reset View</button>
        </div>
      </div>
      <div id="hist-svg-container" class="viz-canvas-hist interactive-viewport" data-zoom-target="hist">
        <div id="hist-stage" class="zoom-stage">
          {hist_svg}
        </div>
      </div>
      <div class="viz-legend" style="padding: 8px 16px 12px; margin-top: 0;">
        <div class="viz-legend-inline">
          <span class="legend-item"><span class="legend-dot" style="background:#059669"></span>&ge;8.0</span>
          <span class="legend-item"><span class="legend-dot" style="background:#0D9488"></span>7.0–7.9</span>
          <span class="legend-item"><span class="legend-dot" style="background:#D97706"></span>5.0–6.9</span>
          <span class="legend-item"><span class="legend-dot" style="background:#E11D48"></span>&lt;5.0</span>
        </div>
      </div>
    </div>

    <!-- Right Cell: Monospace Bioactivity Table -->
    <div class="bento-cell no-padding">
      <div class="table-scroll-container">
        <table id="activity-table" class="data-table">
          <thead>
            <tr>
              <th>ID</th>
              <th>Target</th>
              <th>Organism</th>
              <th>Type</th>
              <th>Value (nM)</th>
              <th>pIC50</th>
              <th>Assay / Alert</th>
            </tr>
          </thead>
          <tbody>
            {rows_html}
          </tbody>
        </table>
      </div>
    </div>
  </div>
</section>
"""
  return card_html


def main() -> None:
  parser = argparse.ArgumentParser(description="Generate bioactivity visualizer card.")
  parser.add_argument("--activity_json", type=str, required=True, help="Path to activity JSON")
  parser.add_argument("--output_html", type=str, required=True, help="Path to output HTML")
  args = parser.parse_args()

  with open(args.activity_json, "r", encoding="utf-8") as f:
    data = json.load(f)

  html_out = render_bioactivity_card(data)
  out_dir = os.path.dirname(args.output_html)
  if out_dir:
    os.makedirs(out_dir, exist_ok=True)
  with open(args.output_html, "w", encoding="utf-8") as f:
    f.write(html_out)
  print(f"Generated bioactivity visualization card at {args.output_html}")


if __name__ == "__main__":
  main()
