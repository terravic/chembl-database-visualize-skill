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
target affinity comparisons, and filterable data tables.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import os
import sys
from typing import Any


def calculate_pic50(normalized_val_nm: float | None) -> float | None:
  """Calculate pIC50 = 9 - log10(nM).

  Args:
    normalized_val_nm: Value in nM.

  Returns:
    Calculated pIC50 rounded to 2 decimals, or None if invalid.
  """
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
      # Try standard_value if normalized_value_nM is missing
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

    # Potency tier
    if std_val is not None:
      if std_val < 100.0:
        tier = "High Potency (<100 nM)"
        potency_counts["high"] += 1
      elif std_val <= 10000.0:
        tier = "Moderate (100 - 10,000 nM)"
        potency_counts["moderate"] += 1
      else:
        tier = "Weak / Inactive (>10,000 nM)"
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
        "document_chembl_id": rec.get("document_chembl_id") or "N/A",
    }
    processed.append(item)

    if pic50 is not None:
      pic50_values.append(pic50)
      t_name = item["target_pref_name"]
      target_affinities.setdefault(t_name, []).append(pic50)

  # Stats
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


def generate_histogram_svg(pic50_values: list[float], width: int = 400, height: int = 220) -> str:
  """Generate an SVG histogram of pIC50 values."""
  if not pic50_values:
    return f"""<svg viewBox="0 0 {width} {height}" class="w-full h-full flex items-center justify-center">
      <rect width="{width}" height="{height}" fill="transparent"/>
      <text x="{width/2}" y="{height/2}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#94a3b8">No numerical pIC50 records available</text>
    </svg>"""

  # Define bins: 3.0 to 11.0 in steps of 1.0 (8 bins)
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
  margin_left = 35
  margin_bottom = 30
  margin_top = 20
  margin_right = 15

  plot_w = width - margin_left - margin_right
  plot_h = height - margin_top - margin_bottom
  bar_w = plot_w / num_bins

  svg_bars = []
  svg_labels = []

  for i, count in enumerate(bins):
    x = margin_left + i * bar_w + 3
    bw = max(2, bar_w - 6)
    bh = (count / max_count) * (plot_h - 10)
    y = margin_top + plot_h - bh

    # Color code by potency: pIC50 >= 7.0 is emerald/blue, 5-7 is amber, <5 is slate
    mid_bin = bin_min + (i + 0.5) * bin_width
    fill_color = "#3b82f6" if mid_bin >= 7.0 else ("#f59e0b" if mid_bin >= 5.0 else "#94a3b8")

    svg_bars.append(
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw:.1f}" height="{bh:.1f}" fill="{fill_color}" rx="3"/>'
    )
    if count > 0:
      svg_bars.append(
          f'<text x="{x + bw/2:.1f}" y="{y - 4:.1f}" font-family="system-ui, sans-serif" font-size="10" font-weight="600" fill="#475569" text-anchor="middle">{count}</text>'
      )

    lbl_x = margin_left + i * bar_w + bar_w / 2
    lbl_text = f"{bin_min + i * bin_width:.0f}"
    svg_labels.append(
        f'<text x="{lbl_x:.1f}" y="{height - 12}" font-family="system-ui, sans-serif" font-size="10" fill="#64748b" text-anchor="middle">{lbl_text}</text>'
    )

  svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="100%">
  <rect width="{width}" height="{height}" fill="transparent"/>
  <!-- Axes lines -->
  <line x1="{margin_left}" y1="{margin_top + plot_h}" x2="{width - margin_right}" y2="{margin_top + plot_h}" stroke="#cbd5e1" stroke-width="1.5"/>
  <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_h}" stroke="#cbd5e1" stroke-width="1.5"/>
  <!-- Bars -->
  {"".join(svg_bars)}
  <!-- X Labels -->
  {"".join(svg_labels)}
  <text x="{width/2}" y="{height - 2}" font-family="system-ui, sans-serif" font-size="10" font-weight="500" fill="#64748b" text-anchor="middle">-log10 Affinity (pIC50 / pKi)</text>
</svg>"""
  return svg


def render_bioactivity_card(activity_data: dict[str, Any]) -> str:
  """Render complete bioactivity visualizer section."""
  raw_list = activity_data.get("activities", [])
  processed = process_activities(raw_list)
  records = processed["records"]
  stats = processed["stats"]
  counts = processed["potency_counts"]
  hist_svg = generate_histogram_svg(processed["pic50_values"])

  # Generate table rows
  table_rows = []
  for r in records[:50]:
    val_str = f"{r['normalized_value_nM']:.2f} nM" if r["normalized_value_nM"] is not None else "N/A"
    pic50_str = f"{r['pIC50']:.2f}" if r["pIC50"] is not None else "-"
    badge_bg = (
        "bg-emerald-100 text-emerald-800 dark:bg-emerald-950 dark:text-emerald-300"
        if "High" in r["potency_tier"]
        else ("bg-amber-100 text-amber-800 dark:bg-amber-950 dark:text-amber-300" if "Moderate" in r["potency_tier"] else "bg-slate-100 text-slate-700 dark:bg-slate-800 dark:text-slate-300")
    )
    table_rows.append(
        f"""<tr class="border-b border-slate-100 dark:border-slate-800/80 hover:bg-slate-50/50 dark:hover:bg-slate-800/30 text-xs">
          <td class="py-2.5 px-3 font-mono font-medium text-slate-700 dark:text-slate-300">{html.escape(str(r['activity_id']))}</td>
          <td class="py-2.5 px-3 text-slate-900 dark:text-slate-100 font-medium max-w-[200px] truncate" title="{html.escape(r['target_pref_name'])}">{html.escape(r['target_pref_name'])}</td>
          <td class="py-2.5 px-3 text-slate-500 italic">{html.escape(r['target_organism'])}</td>
          <td class="py-2.5 px-3 font-semibold text-blue-600 dark:text-blue-400">{html.escape(r['standard_type'])}</td>
          <td class="py-2.5 px-3 font-mono text-slate-800 dark:text-slate-200">{html.escape(r['standard_relation'])} {val_str}</td>
          <td class="py-2.5 px-3 font-mono font-semibold text-slate-900 dark:text-white">{pic50_str}</td>
          <td class="py-2.5 px-3"><span class="px-2 py-0.5 rounded-full text-[10px] font-semibold {badge_bg}">{html.escape(r['potency_tier'])}</span></td>
          <td class="py-2.5 px-3 font-mono text-slate-400">{html.escape(r['assay_chembl_id'])}</td>
        </tr>"""
    )

  rows_html = (
      "".join(table_rows)
      if table_rows
      else '<tr><td colspan="8" class="py-6 text-center text-slate-400">No bioactivity measurements found</td></tr>'
  )

  card_html = f"""
<div class="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl shadow-sm overflow-hidden mb-8">
  <div class="px-6 py-4 border-b border-slate-100 dark:border-slate-800 flex flex-wrap items-center justify-between gap-3 bg-slate-50/50 dark:bg-slate-800/40">
    <div class="flex items-center space-x-3">
      <span class="px-2.5 py-1 text-xs font-mono font-semibold rounded-md bg-purple-100 dark:bg-purple-950 text-purple-800 dark:text-purple-300 border border-purple-200 dark:border-purple-900">SAR Profile</span>
      <h2 class="text-xl font-bold text-slate-900 dark:text-white tracking-tight">Bioactivity & Potency Distribution</h2>
    </div>
    <div class="flex items-center space-x-2 text-xs">
      <span class="px-2.5 py-1 rounded-md bg-slate-100 dark:bg-slate-800 text-slate-700 dark:text-slate-300 font-medium">Assays: {len(records)}</span>
      <span class="px-2.5 py-1 rounded-md bg-emerald-50 dark:bg-emerald-950/60 text-emerald-700 dark:text-emerald-400 font-medium font-mono">High Potency: {counts['high']}</span>
    </div>
  </div>

  <div class="p-6 grid grid-cols-1 lg:grid-cols-12 gap-6">
    <!-- Left Column: Summary Metrics & Histogram -->
    <div class="lg:col-span-5 flex flex-col space-y-4">
      <div class="grid grid-cols-3 gap-2">
        <div class="p-3 rounded-lg bg-slate-50 dark:bg-slate-800/60 border border-slate-100 dark:border-slate-800 text-center">
          <span class="text-[11px] text-slate-500 font-medium block">Median pIC50</span>
          <span class="text-lg font-bold font-mono text-slate-800 dark:text-white">{stats['median']}</span>
        </div>
        <div class="p-3 rounded-lg bg-slate-50 dark:bg-slate-800/60 border border-slate-100 dark:border-slate-800 text-center">
          <span class="text-[11px] text-slate-500 font-medium block">Max Potency</span>
          <span class="text-lg font-bold font-mono text-emerald-600 dark:text-emerald-400">{stats['max']}</span>
        </div>
        <div class="p-3 rounded-lg bg-slate-50 dark:bg-slate-800/60 border border-slate-100 dark:border-slate-800 text-center">
          <span class="text-[11px] text-slate-500 font-medium block">Total Screened</span>
          <span class="text-lg font-bold font-mono text-blue-600 dark:text-blue-400">{stats['count']}</span>
        </div>
      </div>

      <div class="border border-slate-200 dark:border-slate-800 rounded-lg p-3 bg-white dark:bg-slate-950">
        <div class="flex items-center justify-between pb-2 mb-2 border-b border-slate-100 dark:border-slate-800">
          <span class="text-xs font-semibold uppercase tracking-wider text-slate-500">Potency Spread (pIC50)</span>
          <span class="text-[11px] text-slate-400 font-mono">Count by Bin</span>
        </div>
        <div class="h-56 flex items-center justify-center">
          {hist_svg}
        </div>
      </div>
    </div>

    <!-- Right Column: Interactive Searchable Table -->
    <div class="lg:col-span-7 flex flex-col">
      <div class="border border-slate-200 dark:border-slate-800 rounded-lg bg-white dark:bg-slate-950 overflow-hidden flex flex-col h-full">
        <div class="p-3 border-b border-slate-100 dark:border-slate-800 flex items-center justify-between gap-3 bg-slate-50/50 dark:bg-slate-900/40">
          <span class="text-xs font-semibold uppercase tracking-wider text-slate-500">Normalized Assay Data (nM)</span>
          <input type="text" id="activity-table-search" placeholder="Filter target, assay..." oninput="filterActivityTable()" class="px-2.5 py-1 text-xs rounded border border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-900 text-slate-800 dark:text-slate-200 focus:outline-none focus:ring-1 focus:ring-blue-500"/>
        </div>
        <div class="overflow-x-auto max-h-[320px] overflow-y-auto">
          <table id="activity-table" class="w-full text-left border-collapse">
            <thead class="bg-slate-50 dark:bg-slate-900 sticky top-0 border-b border-slate-200 dark:border-slate-800 text-[11px] uppercase tracking-wider text-slate-500">
              <tr>
                <th class="py-2 px-3">ID</th>
                <th class="py-2 px-3">Target</th>
                <th class="py-2 px-3">Organism</th>
                <th class="py-2 px-3">Type</th>
                <th class="py-2 px-3">Standard Val (nM)</th>
                <th class="py-2 px-3">pIC50</th>
                <th class="py-2 px-3">Classification</th>
                <th class="py-2 px-3">Assay</th>
              </tr>
            </thead>
            <tbody>
              {rows_html}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  </div>
</div>
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
