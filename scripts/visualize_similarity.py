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

"""Chemical similarity, concordance heatmap matrix, and analog gallery component.

Transforms ChEMBL similarity and substructure search results into a high-contrast
Vivid Indigo / Vibrant Teal heatmap matrix and flat Bento analog cells with
real-time client-side filtering sliders strictly adhering to DESIGN.md.
"""

from __future__ import annotations

import argparse
import html
import json
import os
from typing import Any


def _heatmap_style(norm_val: float, palette: str = "teal") -> str:
  """Compute high-contrast gradient cell background & text color per DESIGN.md Section 5.

  Transitions from a soft desaturated tint into saturated Vivid Indigo (#4F46E5)
  or Vibrant Teal (#0D9488).
  """
  clamped = max(0.05, min(1.0, norm_val))
  alpha = round(0.10 + clamped * 0.78, 2)
  if palette == "indigo":
    bg = f"rgba(79, 70, 229, {alpha})"
  else:
    bg = f"rgba(13, 148, 136, {alpha})"
  text_col = "#FFFFFF" if alpha >= 0.52 else "var(--text-primary)"
  return f"background-color:{bg};color:{text_col};"


def render_similarity_gallery(similarity_data: dict[str, Any]) -> str:
  """Render high-contrast similarity heatmap matrix and flat Bento analog grid."""
  molecules = similarity_data.get("molecules", [])
  total_count = len(molecules)

  cards = []
  heatmap_rows = []

  for idx, mol in enumerate(molecules):
    cid = mol.get("molecule_chembl_id") or f"ANALOG_{idx+1}"
    name = mol.get("pref_name") or "Unspecified Analog"
    sim = float(mol.get("similarity") or 0.0)
    sim_pct = f"{sim:.1f}%" if sim > 0 else "N/A"

    structures = mol.get("molecule_structures") or {}
    smiles = structures.get("canonical_smiles") or ""
    inchi_key = structures.get("standard_inchi_key") or ""

    props = mol.get("molecule_properties") or {}
    mwt = float(props.get("full_mwt") or 0.0)
    alogp = float(props.get("alogp") or 0.0)
    psa = float(props.get("psa") or 0.0)
    hbd = int(props.get("hbd") or 0)
    hba = int(props.get("hba") or 0)
    max_phase = mol.get("max_phase")
    phase_str = f"Phase {max_phase}" if max_phase is not None else "Preclinical"

    # DESIGN.md palette bar colors: Functional Emerald (#059669), Vibrant Teal (#0D9488), Energetic Amber (#D97706)
    if sim >= 90.0:
      bar_color = "#059669"
    elif sim >= 84.0:
      bar_color = "#0D9488"
    else:
      bar_color = "#D97706"

    # Heatmap row normalized intensities
    sim_norm = max(0.0, min(1.0, (sim - 60.0) / 40.0))
    mwt_norm = max(0.0, min(1.0, mwt / 500.0))
    alogp_norm = max(0.0, min(1.0, alogp / 5.0))
    psa_norm = max(0.0, min(1.0, psa / 140.0))
    hb_norm = max(0.0, min(1.0, (hbd + hba) / 10.0))

    heatmap_rows.append(
        f"""<tr class="table-row analog-heatmap-row" data-similarity="{sim:.1f}" data-mwt="{mwt:.1f}">
          <td class="cell-mono" style="font-weight:600">{html.escape(cid)}</td>
          <td class="cell-mono">{html.escape(name)}</td>
          <td class="cell-mono heatmap-cell" style="{_heatmap_style(sim_norm, 'teal')}">{sim:.1f}%</td>
          <td class="cell-mono heatmap-cell" style="{_heatmap_style(mwt_norm, 'indigo')}">{mwt:.1f}</td>
          <td class="cell-mono heatmap-cell" style="{_heatmap_style(alogp_norm, 'indigo')}">{alogp:.2f}</td>
          <td class="cell-mono heatmap-cell" style="{_heatmap_style(psa_norm, 'teal')}">{psa:.1f}</td>
          <td class="cell-mono heatmap-cell" style="{_heatmap_style(hb_norm, 'indigo')}">{hbd} / {hba}</td>
        </tr>"""
    )

    cards.append(
        f"""
        <div class="analog-card" data-similarity="{sim:.1f}" data-mwt="{mwt:.1f}">
          <div class="analog-card-top">
            <div class="analog-id-row">
              <span class="mono-id">{html.escape(cid)}</span>
              <span class="mono-match" style="color:{bar_color}">{sim_pct} Match</span>
            </div>
            <h4 class="analog-name" title="{html.escape(name)}">{html.escape(name)}</h4>
            <div class="sim-bar-track">
              <div class="sim-bar-fill" style="width:{min(100.0, sim):.1f}%;background-color:{bar_color}"></div>
            </div>
            <div class="analog-smiles" title="{html.escape(smiles)}">
              {html.escape(smiles) if smiles else "Structure coordinates unavailable"}
            </div>
          </div>

          <div class="analog-metrics-row">
            <div class="analog-metric-pair">
              <span class="descriptor-label">MW</span>
              <span class="mono-val">{mwt:.1f} Da</span>
            </div>
            <div class="analog-metric-pair">
              <span class="descriptor-label">AlogP</span>
              <span class="mono-val">{alogp:.2f}</span>
            </div>
            <div class="analog-metric-pair">
              <span class="descriptor-label">Stage</span>
              <span class="mono-val">{html.escape(phase_str)}</span>
            </div>
          </div>
          <div class="analog-footer-row">
            <span class="mono-secondary" title="{html.escape(inchi_key)}">{html.escape(inchi_key)}</span>
          </div>
        </div>
        """
    )

  cards_html = (
      "".join(cards)
      if cards
      else '<div class="empty-state">No chemical analogs found for the specified query.</div>'
  )
  heatmap_html = (
      "".join(heatmap_rows)
      if heatmap_rows
      else '<tr><td colspan="7" class="empty-table-cell">No analog records available</td></tr>'
  )

  gallery_html = f"""
<section class="bento-panel" id="similarity-section">
  <div class="bento-header">
    <div class="header-title-group">
      <h2 class="panel-title">Chemical Similarity &amp; Descriptor Matrix</h2>
      <span class="mono-meta">Showing <span id="visible-analog-count">{total_count}</span> of {total_count} analogs</span>
    </div>
    <div class="filter-controls-bar">
      <div class="slider-group">
        <label for="similarity-filter" class="slider-label">Min Similarity</label>
        <input type="range" id="similarity-filter" min="50" max="100" value="70" step="1" oninput="updateAnalogFilters()" class="range-slider"/>
        <span id="similarity-filter-val" class="mono-slider-val">70%</span>
      </div>
      <div class="slider-group">
        <label for="mwt-filter" class="slider-label">Max MW</label>
        <input type="range" id="mwt-filter" min="100" max="800" value="800" step="25" oninput="updateAnalogFilters()" class="range-slider"/>
        <span id="mwt-filter-val" class="mono-slider-val">800 Da</span>
      </div>
      <button type="button" onclick="resetAnalogFilters()" class="btn-control">Reset</button>
    </div>
  </div>

  <!-- High-Contrast Physicochemical & Similarity Heatmap Matrix -->
  <div class="bento-subcell no-padding">
    <div class="subcell-toolbar padded-toolbar">
      <h3 class="cell-title">Normalized Descriptor Heatmap Matrix</h3>
      <div class="viz-legend-inline">
        <span class="legend-item"><span class="legend-dot" style="background:#0D9488"></span>Similarity / TPSA Density</span>
        <span class="legend-item"><span class="legend-dot" style="background:#4F46E5"></span>MW / Lipophilicity / H-Bond Scale</span>
      </div>
    </div>
    <div class="table-scroll-container">
      <table class="data-table">
        <thead>
          <tr>
            <th>ChEMBL ID</th>
            <th>Compound</th>
            <th>Tanimoto Similarity</th>
            <th>MW (Da)</th>
            <th>AlogP</th>
            <th>TPSA (&Aring;&sup2;)</th>
            <th>HBD / HBA</th>
          </tr>
        </thead>
        <tbody>
          {heatmap_html}
        </tbody>
      </table>
    </div>
  </div>

  <!-- Flat Hairline Analog Grid -->
  <div id="analog-grid" class="analog-bento-grid">
    {cards_html}
  </div>
</section>
"""
  return gallery_html


def main() -> None:
  parser = argparse.ArgumentParser(description="Generate chemical similarity gallery.")
  parser.add_argument("--similarity_json", type=str, required=True, help="Path to similarity JSON")
  parser.add_argument("--output_html", type=str, required=True, help="Path to output HTML")
  args = parser.parse_args()

  with open(args.similarity_json, "r", encoding="utf-8") as f:
    data = json.load(f)

  gallery = render_similarity_gallery(data)
  out_dir = os.path.dirname(args.output_html)
  if out_dir:
    os.makedirs(out_dir, exist_ok=True)
  with open(args.output_html, "w", encoding="utf-8") as f:
    f.write(gallery)
  print(f"Generated similarity gallery card at {args.output_html}")


if __name__ == "__main__":
  main()
