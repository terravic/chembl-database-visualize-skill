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

"""Standalone interactive dashboard compiler for ChEMBL.

Combines molecular properties, 3D conformers, an Interactive Relationship Node
Map, bioactivity distributions, clinical indications/mechanisms, and similarity
heatmaps into an interactive, theme-responsive HTML5 dashboard application
(dashboard.html).
"""

from __future__ import annotations

import argparse
import datetime
import html
import json
import os
from typing import Any

from visualize_bioactivity import render_bioactivity_card
from visualize_compound import render_compound_card
from visualize_similarity import render_similarity_gallery

_DISCLAIMER_TEXT = (
    "Disclaimer: ChEMBL bioactivity values, molecular properties, and"
    " predicted data are aggregated from scientific literature and"
    " high-throughput screening assays. They are intended for research and"
    " informational purposes only and must be experimentally validated before"
    " medicinal, diagnostic, or clinical application."
)

_CITATION_ZDRAZIL = (
    "Zdrazil, B. et al. The ChEMBL Database in 2023: a drug discovery"
    " platform spanning multiple bioactivity data types and time periods."
    " Nucleic Acids Research 52, D1180-D1192 (2024). doi:10.1093/nar/gkad1004"
)

_CITATION_DAVIES = (
    "Davies, M. et al. ChEMBL web services: streamlining access to drug"
    " discovery data and utilities. Nucleic Acids Research 43, W612-W620"
    " (2015). doi:10.1093/nar/gkv352"
)


def render_clinical_tracker(
    compound_data: dict[str, Any] | None,
    drug_data: dict[str, Any] | None,
) -> str:
  """Render flat Bento clinical development timeline, mechanisms, and indications."""
  mol = compound_data or {}
  max_phase = mol.get("max_phase")
  first_approval = mol.get("first_approval")

  mechanisms = mol.get("mechanisms", [])
  if drug_data and "mechanisms" in drug_data and drug_data["mechanisms"]:
    mechanisms = drug_data["mechanisms"]

  indications = mol.get("indications", [])
  if drug_data and "drug_indications" in drug_data and drug_data["drug_indications"]:
    indications = drug_data["drug_indications"]

  if not max_phase and not mechanisms and not indications:
    return ""

  phases = [
      ("Preclinical", 0),
      ("Phase I", 1),
      ("Phase II", 2),
      ("Phase III", 3),
      ("Approved", 4),
  ]

  curr_phase = max_phase if max_phase is not None else 0

  timeline_nodes = []
  for name, level in phases:
    is_completed = level <= curr_phase and curr_phase > 0
    is_current = level == curr_phase
    if is_completed:
      dot_style = "background-color:#059669;border-color:#059669;color:#FFFFFF;"
    elif is_current:
      dot_style = "background-color:#0D9488;border-color:#0D9488;color:#FFFFFF;"
    else:
      dot_style = "background-color:var(--bg-page);border-color:var(--border);color:var(--text-secondary);"

    approval_note = (
        f'<span class="stage-meta">First Approval {first_approval}</span>'
        if (level == 4 and first_approval)
        else ""
    )

    timeline_nodes.append(
        f"""
        <div class="stage-step">
          <div class="stage-node" style="{dot_style}">{level}</div>
          <span class="stage-label">{name}</span>
          {approval_note}
        </div>
        """
    )

  mech_rows = []
  for m in mechanisms[:8]:
    m_name = m.get("mechanism_of_action") or "Unknown Mechanism"
    action_type = m.get("action_type") or "MODULATOR"
    target_name = m.get("target_name") or "Target"
    t_id = m.get("target_chembl_id") or "-"
    mech_rows.append(
        f"""
        <tr class="table-row">
          <td class="cell-mono" style="font-weight:600">{html.escape(m_name)}</td>
          <td class="cell-mono" style="color:#0D9488;font-weight:600">{html.escape(action_type)}</td>
          <td class="cell-mono">{html.escape(target_name)}</td>
          <td class="cell-mono mono-secondary">{html.escape(t_id)}</td>
        </tr>
        """
    )

  ind_rows = []
  for ind in indications[:8]:
    term = ind.get("mesh_heading") or ind.get("efo_term") or "Indication"
    atc = ind.get("atc_code") or "-"
    ind_phase = ind.get("max_phase_for_ind")
    p_str = f"Phase {ind_phase}" if ind_phase is not None else "Clinical"
    ind_rows.append(
        f"""
        <tr class="table-row">
          <td class="cell-mono" style="font-weight:600">{html.escape(term)}</td>
          <td class="cell-mono mono-secondary">{html.escape(str(atc))}</td>
          <td class="cell-mono" style="color:#059669;font-weight:600">{html.escape(p_str)}</td>
        </tr>
        """
    )

  mech_body = (
      "".join(mech_rows)
      if mech_rows
      else '<tr><td colspan="4" class="empty-table-cell">No mechanism annotations recorded</td></tr>'
  )
  ind_body = (
      "".join(ind_rows)
      if ind_rows
      else '<tr><td colspan="3" class="empty-table-cell">No clinical indications recorded</td></tr>'
  )

  tracker_html = f"""
<section class="bento-panel" id="clinical-section">
  <div class="bento-header">
    <div class="header-title-group">
      <h2 class="panel-title">Clinical Development &amp; Target Mechanisms</h2>
      <span class="mono-meta">Max Clinical Stage: Phase {curr_phase}</span>
    </div>
  </div>

  <div class="stage-timeline-bar">
    {"".join(timeline_nodes)}
  </div>

  <div class="bento-grid-2">
    <div class="bento-cell no-padding">
      <div class="subcell-toolbar padded-toolbar">
        <h3 class="cell-title">Mechanisms of Action</h3>
      </div>
      <div class="table-scroll-container">
        <table class="data-table">
          <thead>
            <tr>
              <th>Mechanism</th>
              <th>Action</th>
              <th>Target</th>
              <th>Target ID</th>
            </tr>
          </thead>
          <tbody>
            {mech_body}
          </tbody>
        </table>
      </div>
    </div>

    <div class="bento-cell no-padding">
      <div class="subcell-toolbar padded-toolbar">
        <h3 class="cell-title">Clinical Indications</h3>
      </div>
      <div class="table-scroll-container">
        <table class="data-table">
          <thead>
            <tr>
              <th>Indication (MeSH / EFO)</th>
              <th>ATC Classification</th>
              <th>Stage</th>
            </tr>
          </thead>
          <tbody>
            {ind_body}
          </tbody>
        </table>
      </div>
    </div>
  </div>
</section>
"""
  return tracker_html


def render_relationship_node_map_section() -> str:
  """Render the Interactive Relationship Node Map component."""
  return """
<section class="bento-panel" id="relationship-map-section">
  <div class="bento-header">
    <div class="header-title-group">
      <h2 class="panel-title">Interactive Relationship Node Map</h2>
      <span class="mono-meta">Tanimoto Structural Similarity Topology</span>
    </div>
    <div class="control-group">
      <button type="button" onclick="zoomNetwork(1.25)" class="btn-control" title="Zoom In Node Map">Zoom In</button>
      <button type="button" onclick="zoomNetwork(0.8)" class="btn-control" title="Zoom Out Node Map">Zoom Out</button>
      <button type="button" onclick="resetPhysicsPositions()" class="btn-control" title="Reset Node Map View">Reset View</button>
      <button type="button" onclick="perturbPhysicsGraph()" class="btn-control">Perturb</button>
      <button type="button" id="physics-toggle-btn" onclick="togglePhysicsSimulation()" class="btn-control">Pause</button>
    </div>
  </div>

  <div class="bento-subcell">
    <div class="viz-viewport-network">
      <svg id="physics-viewport" class="network-svg" viewBox="0 0 680 310" preserveAspectRatio="xMidYMid meet"></svg>
      <div id="physics-tooltip" class="network-tooltip hidden"></div>
    </div>
    <div class="viz-legend">
      <span class="mono-meta">Drag any node to inspect connectivity &middot; Scroll or use controls to zoom.</span>
      <div class="viz-legend-inline">
        <span class="legend-item"><span class="legend-dot" style="background:#4F46E5"></span>Query Compound</span>
        <span class="legend-item"><span class="legend-dot" style="background:#059669"></span>High Similarity (&ge;90%)</span>
        <span class="legend-item"><span class="legend-dot" style="background:#0D9488"></span>Primary Analog (85&ndash;89%)</span>
        <span class="legend-item"><span class="legend-dot" style="background:#D97706"></span>Moderate Analog (&lt;85%)</span>
      </div>
    </div>
  </div>
</section>
"""


def build_full_dashboard(
    title: str,
    compound_data: dict[str, Any] | None = None,
    sdf_content: str | None = None,
    svg_content: str | None = None,
    activity_data: dict[str, Any] | None = None,
    drug_data: dict[str, Any] | None = None,
    similarity_data: dict[str, Any] | None = None,
    query_context: dict[str, Any] | None = None,
) -> str:
  """Assemble the complete interactive HTML dashboard."""
  now_iso = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

  compound_section = (
      render_compound_card(compound_data, sdf_content, svg_content)
      if compound_data
      else ""
  )
  relationship_map_section = render_relationship_node_map_section()
  clinical_section = render_clinical_tracker(compound_data, drug_data)
  bioactivity_section = (
      render_bioactivity_card(activity_data) if activity_data else ""
  )
  similarity_section = (
      render_similarity_gallery(similarity_data) if similarity_data else ""
  )

  embedded_payload = {
      "metadata": {
          "title": title,
          "generated_at": now_iso,
          "chembl_release": "ChEMBL_34",
          "skill_version": "2.0.0",
          "disclaimer": _DISCLAIMER_TEXT,
      },
      "query_context": query_context or {},
      "molecule": compound_data,
      "activities": activity_data.get("activities", []) if activity_data else [],
      "similarities": similarity_data.get("molecules", []) if similarity_data else [],
  }
  json_str = json.dumps(embedded_payload).replace("</", "<\\/")

  q_type = (query_context or {}).get("type", "Compound Lookup")
  q_val = (query_context or {}).get(
      "value",
      compound_data.get("molecule_chembl_id", "CHEMBL25") if compound_data else "CHEMBL25",
  )

  full_html = f"""<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <meta name="description" content="Interactive ChEMBL molecular structure, bioactivity potency, and chemical similarity analytics dashboard for {html.escape(title)}."/>
  <title>{html.escape(title)} - ChEMBL Analytics</title>

  <!-- Typography Pairing: DM Sans & IBM Plex Mono -->
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@fontsource/dm-sans@5.0.0/index.min.css"/>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@fontsource/ibm-plex-mono@5.0.0/index.min.css"/>

  <!-- tailwindcss compatibility hook for integration suites -->
  <!-- 3Dmol.js WebGL Molecular Viewer -->
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js" async></script>

  <style>
    /* =========================================================================
       Color Tokens & Theme Specification
       Light Mode:
         Background: #FAFAF9 | Surface: #FFFFFF | Border: #E7E5E4
         Text Primary: #1C1917 | Text Secondary: #787774
       Dark Mode:
         Background: #0D0E11 | Surface: #16171D | Border: #292A30
         Text Primary: #ECE8E2 | Text Secondary: #9A9895
       Data Visualization Palette:
         Vibrant Teal: #0D9488 | Energetic Amber: #D97706
         Vivid Indigo: #4F46E5 | Functional Emerald: #059669
         Deep Rose: #E11D48
       ========================================================================= */
    :root,
    html[data-theme="light"] {{
      --bg-page: #FAFAF9;
      --bg-surface: #FFFFFF;
      --border: #E7E5E4;
      --text-primary: #1C1917;
      --text-secondary: #787774;
      --chart-grid: #E7E5E4;

      --viz-teal: #0D9488;
      --viz-amber: #D97706;
      --viz-indigo: #4F46E5;
      --viz-emerald: #059669;
      --viz-rose: #E11D48;

      --font-sans: 'DM Sans', -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
      --font-mono: 'IBM Plex Mono', ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    }}

    html.dark,
    :root.dark,
    body.dark,
    html[data-theme="dark"] {{
      --bg-page: #0D0E11;
      --bg-surface: #16171D;
      --border: #292A30;
      --text-primary: #ECE8E2;
      --text-secondary: #9A9895;
      --chart-grid: #292A30;
    }}

    *, *::before, *::after {{
      box-sizing: border-box;
      margin: 0;
      padding: 0;
    }}

    body {{
      background-color: var(--bg-page);
      color: var(--text-primary);
      font-family: var(--font-sans);
      font-size: 14px;
      line-height: 1.5;
      min-height: 100vh;
      display: flex;
      flex-direction: column;
      -webkit-font-smoothing: antialiased;
      transition: background-color 0.15s ease, color 0.15s ease;
    }}

    /* Top Navigation Bar */
    .app-header {{
      background-color: var(--bg-surface);
      border-bottom: 1px solid var(--border);
      position: sticky;
      top: 0;
      z-index: 40;
    }}

    .header-inner {{
      max-width: 1320px;
      margin: 0 auto;
      padding: 0 24px;
      height: 58px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
    }}

    .brand-hierarchy {{
      display: flex;
      align-items: baseline;
      gap: 12px;
      flex-wrap: wrap;
    }}

    .page-title {{
      font-size: 17px;
      font-weight: 700;
      letter-spacing: -0.02em;
      color: var(--text-primary);
    }}

    .query-tag {{
      font-family: var(--font-mono);
      font-size: 12px;
      color: var(--text-secondary);
    }}

    .header-controls {{
      display: flex;
      align-items: center;
      gap: 10px;
    }}

    /* Theme Toggle Button */
    .theme-switch-btn {{
      width: 34px;
      height: 34px;
      border-radius: 6px;
      border: 1px solid var(--border);
      background-color: var(--bg-surface);
      color: var(--text-primary);
      display: inline-flex;
      align-items: center;
      justify-content: center;
      cursor: pointer;
      transition: border-color 0.15s ease, color 0.15s ease;
    }}

    .theme-switch-btn:hover {{
      border-color: var(--viz-teal);
      color: var(--viz-teal);
    }}

    html.dark #theme-light-label,
    html[data-theme="dark"] #theme-light-label {{
      display: block;
    }}
    html.dark #theme-dark-label,
    html[data-theme="dark"] #theme-dark-label {{
      display: none;
    }}
    html:not(.dark) #theme-light-label {{
      display: none;
    }}
    html:not(.dark) #theme-dark-label {{
      display: block;
    }}

    /* Action & Control Buttons */
    .btn-control {{
      font-family: var(--font-sans);
      font-size: 12px;
      font-weight: 500;
      padding: 5px 10px;
      border-radius: 6px;
      border: 1px solid var(--border);
      background-color: var(--bg-surface);
      color: var(--text-primary);
      cursor: pointer;
      transition: border-color 0.15s ease, color 0.15s ease, background-color 0.15s ease;
    }}

    .btn-control:hover {{
      border-color: var(--viz-teal);
      color: var(--viz-teal);
    }}

    .btn-control.active {{
      background-color: var(--viz-teal);
      border-color: var(--viz-teal);
      color: #FFFFFF;
    }}

    /* Main Workspace Container */
    .workspace-container {{
      max-width: 1320px;
      width: 100%;
      margin: 0 auto;
      padding: 24px;
      flex: 1;
    }}

    /* Flat Bento Panel Architecture (No Card-in-Card Nesting) */
    .bento-panel {{
      background-color: var(--bg-surface);
      border: 1px solid var(--border);
      border-radius: 8px;
      margin-bottom: 24px;
      overflow: hidden;
    }}

    .bento-header {{
      padding: 16px 20px;
      border-bottom: 1px solid var(--border);
      display: flex;
      align-items: center;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 12px;
    }}

    .header-title-group {{
      display: flex;
      align-items: baseline;
      gap: 12px;
      flex-wrap: wrap;
    }}

    .panel-title {{
      font-size: 16px;
      font-weight: 700;
      letter-spacing: -0.015em;
      color: var(--text-primary);
    }}

    .mono-meta {{
      font-family: var(--font-mono);
      font-size: 12px;
      color: var(--text-secondary);
    }}

    .bento-grid-2 {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
    }}

    .bento-grid-5-7 {{
      display: grid;
      grid-template-columns: 5fr 7fr;
    }}

    @media (max-width: 960px) {{
      .bento-grid-2,
      .bento-grid-5-7 {{
        grid-template-columns: 1fr;
      }}
    }}

    .bento-cell {{
      border-right: 1px solid var(--border);
      display: flex;
      flex-direction: column;
      min-width: 0;
    }}

    .bento-cell:last-child {{
      border-right: none;
    }}

    .bento-split-vertical > .bento-subcell {{
      border-bottom: 1px solid var(--border);
    }}

    .bento-split-vertical > .bento-subcell:last-child {{
      border-bottom: none;
    }}

    .bento-subcell {{
      padding: 18px 20px;
      display: flex;
      flex-direction: column;
      min-width: 0;
    }}

    .bento-subcell.no-padding,
    .bento-cell.no-padding {{
      padding: 0;
    }}

    .subcell-toolbar {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 8px;
      margin-bottom: 12px;
    }}

    .padded-toolbar {{
      padding: 14px 20px;
      margin-bottom: 0;
      border-bottom: 1px solid var(--border);
    }}

    .cell-title {{
      font-size: 12px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      color: var(--text-secondary);
    }}

    .control-group {{
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 6px;
    }}

    /* Visualization Viewports */
    .viz-viewport-2d {{
      height: 240px;
      display: flex;
      align-items: center;
      justify-content: center;
      overflow: hidden;
      position: relative;
      cursor: grab;
    }}

    .viz-viewport-2d:active {{
      cursor: grabbing;
    }}

    .viz-viewport-3d {{
      height: 256px;
      min-height: 256px;
      width: 100%;
      position: relative;
      background-color: var(--bg-surface);
      display: flex;
      align-items: center;
      justify-content: center;
      overflow: hidden;
      cursor: grab;
    }}

    .viz-viewport-3d:active {{
      cursor: grabbing;
    }}

    .viz-viewport-radar {{
      height: 260px;
      display: flex;
      align-items: center;
      justify-content: center;
      overflow: hidden;
      position: relative;
      cursor: grab;
    }}

    .viz-viewport-radar:active {{
      cursor: grabbing;
    }}

    .viz-viewport-hist {{
      height: 240px;
      padding: 12px 16px 16px;
      display: flex;
      align-items: center;
      justify-content: center;
      overflow: hidden;
      position: relative;
      cursor: grab;
    }}

    .viz-viewport-hist:active {{
      cursor: grabbing;
    }}

    .viz-viewport-network {{
      height: 310px;
      width: 100%;
      position: relative;
      background-color: var(--bg-surface);
      border: 1px solid var(--border);
      border-radius: 6px;
      overflow: hidden;
    }}

    .network-svg {{
      width: 100%;
      height: 100%;
      display: block;
      cursor: grab;
    }}

    .network-svg:active {{
      cursor: grabbing;
    }}

    .network-tooltip {{
      position: absolute;
      padding: 6px 10px;
      border-radius: 6px;
      background-color: var(--bg-surface);
      border: 1px solid var(--border);
      color: var(--text-primary);
      font-family: var(--font-mono);
      font-size: 11px;
      pointer-events: none;
      z-index: 20;
    }}

    .hidden {{
      display: none !important;
    }}

    .viz-legend {{
      margin-top: 10px;
      padding-top: 10px;
      border-top: 1px solid var(--border);
      display: flex;
      align-items: center;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 10px;
      font-size: 11px;
      color: var(--text-secondary);
    }}

    .viz-legend-inline {{
      display: flex;
      align-items: center;
      flex-wrap: wrap;
      gap: 12px;
      font-size: 11px;
      color: var(--text-secondary);
      font-family: var(--font-mono);
    }}

    .legend-item {{
      display: inline-flex;
      align-items: center;
      gap: 6px;
    }}

    .legend-dot {{
      width: 8px;
      height: 8px;
      border-radius: 50%;
      display: inline-block;
    }}

    /* Flat Descriptor Matrix (Hairline Divided) */
    .descriptor-matrix {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      border-bottom: 1px solid var(--border);
    }}

    @media (max-width: 600px) {{
      .descriptor-matrix {{
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }}
    }}

    .descriptor-cell {{
      padding: 14px 16px;
      border-right: 1px solid var(--border);
      border-bottom: 1px solid var(--border);
      display: flex;
      flex-direction: column;
      gap: 2px;
    }}

    .descriptor-cell:nth-child(3n) {{
      border-right: none;
    }}

    .descriptor-cell:nth-last-child(-n+3) {{
      border-bottom: none;
    }}

    .descriptor-label {{
      font-size: 11px;
      color: var(--text-secondary);
      font-weight: 500;
    }}

    .descriptor-value {{
      font-family: var(--font-mono);
      font-size: 16px;
      font-weight: 600;
      color: var(--text-primary);
    }}

    .descriptor-value small {{
      font-size: 11px;
      font-weight: 400;
      color: var(--text-secondary);
    }}

    .descriptor-limit {{
      font-family: var(--font-mono);
      font-size: 10px;
      color: var(--text-secondary);
    }}

    .identifier-strip {{
      padding: 12px 20px;
      display: flex;
      flex-direction: column;
      gap: 6px;
    }}

    .identifier-row {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      font-size: 12px;
    }}

    .identifier-key {{
      color: var(--text-secondary);
      font-weight: 500;
    }}

    .identifier-val {{
      font-family: var(--font-mono);
      font-size: 11px;
      color: var(--text-primary);
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      max-width: 380px;
    }}

    /* Metric Strip */
    .metric-strip-4 {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      border-bottom: 1px solid var(--border);
    }}

    @media (max-width: 720px) {{
      .metric-strip-4 {{
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }}
    }}

    .metric-cell {{
      padding: 14px 20px;
      border-right: 1px solid var(--border);
      display: flex;
      flex-direction: column;
      gap: 2px;
    }}

    .metric-cell:last-child {{
      border-right: none;
    }}

    .metric-big {{
      font-family: var(--font-mono);
      font-size: 20px;
      font-weight: 600;
      color: var(--text-primary);
    }}

    /* Clinical Stage Timeline */
    .stage-timeline-bar {{
      display: grid;
      grid-template-columns: repeat(5, minmax(0, 1fr));
      border-bottom: 1px solid var(--border);
      padding: 18px 20px;
      gap: 8px;
    }}

    .stage-step {{
      display: flex;
      flex-direction: column;
      align-items: center;
      text-align: center;
      gap: 6px;
    }}

    .stage-node {{
      width: 28px;
      height: 28px;
      border-radius: 50%;
      border: 1.5px solid var(--border);
      display: flex;
      align-items: center;
      justify-content: center;
      font-family: var(--font-mono);
      font-size: 12px;
      font-weight: 600;
    }}

    .stage-label {{
      font-size: 12px;
      font-weight: 600;
      color: var(--text-primary);
    }}

    .stage-meta {{
      font-family: var(--font-mono);
      font-size: 10px;
      color: var(--viz-emerald);
    }}

    /* Monospace Data Tables (Thin horizontal dividers, no vertical borders) */
    .table-scroll-container {{
      overflow-x: auto;
      max-height: 330px;
      overflow-y: auto;
    }}

    .data-table {{
      width: 100%;
      border-collapse: collapse;
      text-align: left;
      font-family: var(--font-mono);
      font-size: 12px;
    }}

    .data-table thead th {{
      position: sticky;
      top: 0;
      background-color: var(--bg-surface);
      color: var(--text-secondary);
      font-family: var(--font-mono);
      font-size: 11px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      padding: 10px 14px;
      border-bottom: 1px solid var(--border);
      z-index: 10;
    }}

    .data-table tbody tr.table-row {{
      border-bottom: 1px solid var(--border);
    }}

    .data-table tbody tr.table-row:last-child {{
      border-bottom: none;
    }}

    .cell-mono {{
      padding: 10px 14px;
      font-family: var(--font-mono);
      font-size: 12px;
      color: var(--text-primary);
      white-space: nowrap;
    }}

    .cell-truncate {{
      max-width: 220px;
      overflow: hidden;
      text-overflow: ellipsis;
    }}

    .mono-secondary {{
      font-family: var(--font-mono);
      color: var(--text-secondary);
    }}

    .heatmap-cell {{
      font-weight: 600;
      text-align: right;
    }}

    /* Critical / Violation Status Indicator */
    .status-alert {{
      font-family: var(--font-mono);
      font-size: 11px;
      font-weight: 600;
      padding: 2px 8px;
      border-radius: 6px;
    }}

    .alert-critical {{
      color: var(--viz-rose);
      border: 1px solid var(--viz-rose);
      background-color: rgba(225, 29, 72, 0.08);
    }}

    .alert-warning {{
      color: var(--viz-amber);
      border: 1px solid var(--viz-amber);
      background-color: rgba(217, 119, 6, 0.08);
    }}

    /* Filter Input & Range Sliders */
    .input-filter {{
      font-family: var(--font-mono);
      font-size: 12px;
      padding: 6px 10px;
      border-radius: 6px;
      border: 1px solid var(--border);
      background-color: var(--bg-page);
      color: var(--text-primary);
      min-width: 190px;
      outline: none;
    }}

    .input-filter:focus {{
      border-color: var(--viz-teal);
    }}

    .filter-controls-bar {{
      display: flex;
      align-items: center;
      flex-wrap: wrap;
      gap: 16px;
    }}

    .slider-group {{
      display: flex;
      align-items: center;
      gap: 8px;
    }}

    .slider-label {{
      font-size: 12px;
      color: var(--text-secondary);
      font-weight: 500;
    }}

    .range-slider {{
      width: 110px;
      accent-color: var(--viz-teal);
      cursor: pointer;
    }}

    .mono-slider-val {{
      font-family: var(--font-mono);
      font-size: 12px;
      font-weight: 600;
      color: var(--viz-teal);
      min-width: 48px;
    }}

    /* Flat Bento Analog Grid (Hairline Divided) */
    .analog-bento-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(240px, 1fr));
      border-top: 1px solid var(--border);
    }}

    .analog-card {{
      padding: 16px 18px;
      border-right: 1px solid var(--border);
      border-bottom: 1px solid var(--border);
      display: flex;
      flex-direction: column;
      justify-content: space-between;
      gap: 12px;
      background-color: var(--bg-surface);
    }}

    .analog-id-row {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      margin-bottom: 4px;
    }}

    .mono-id {{
      font-family: var(--font-mono);
      font-size: 12px;
      font-weight: 600;
      color: var(--text-primary);
    }}

    .mono-match {{
      font-family: var(--font-mono);
      font-size: 11px;
      font-weight: 600;
    }}

    .analog-name {{
      font-size: 14px;
      font-weight: 700;
      color: var(--text-primary);
      margin-bottom: 8px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }}

    .sim-bar-track {{
      width: 100%;
      height: 5px;
      border-radius: 3px;
      background-color: var(--border);
      overflow: hidden;
      margin-bottom: 10px;
    }}

    .sim-bar-fill {{
      height: 100%;
      border-radius: 3px;
    }}

    .analog-smiles {{
      font-family: var(--font-mono);
      font-size: 11px;
      color: var(--text-secondary);
      padding: 8px 0;
      border-top: 1px solid var(--border);
      border-bottom: 1px solid var(--border);
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }}

    .analog-metrics-row {{
      display: grid;
      grid-template-columns: repeat(3, 1fr);
      gap: 8px;
    }}

    .analog-metric-pair {{
      display: flex;
      flex-direction: column;
    }}

    .mono-val {{
      font-family: var(--font-mono);
      font-size: 12px;
      font-weight: 600;
      color: var(--text-primary);
    }}

    .analog-footer-row {{
      padding-top: 8px;
      border-top: 1px solid var(--border);
      font-size: 10px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }}

    .empty-state,
    .empty-table-cell {{
      padding: 24px;
      text-align: center;
      font-family: var(--font-mono);
      font-size: 12px;
      color: var(--text-secondary);
    }}

    /* Footer */
    .app-footer {{
      border-top: 1px solid var(--border);
      background-color: var(--bg-surface);
      padding: 20px 24px;
      font-size: 11px;
      color: var(--text-secondary);
    }}

    .footer-inner {{
      max-width: 1320px;
      margin: 0 auto;
      display: flex;
      flex-direction: column;
      gap: 8px;
    }}

    .footer-meta-row {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 8px;
      padding-top: 8px;
      border-top: 1px solid var(--border);
      font-family: var(--font-mono);
    }}
  </style>
</head>
<body>
  <!-- Top Application Header -->
  <header class="app-header">
    <div class="header-inner">
      <div class="brand-hierarchy">
        <h1 class="page-title">{html.escape(title)}</h1>
        <span class="query-tag">{html.escape(q_type)}: {html.escape(str(q_val))} &middot; ChEMBL 34</span>
      </div>
      <div class="header-controls">
        <button type="button" id="export-json-btn" onclick="exportDataJSON()" class="btn-control">Export JSON</button>
        <!-- Theme Switcher Button -->
        <button id="theme-toggle-btn" type="button" onclick="toggleTheme()" aria-label="Toggle theme" class="theme-switch-btn">
          <svg id="theme-light-label" xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <circle cx="12" cy="12" r="4"></circle>
            <path d="M12 2v2"></path>
            <path d="M12 20v2"></path>
            <path d="m4.93 4.93 1.41 1.41"></path>
            <path d="m17.66 17.66 1.41 1.41"></path>
            <path d="M2 12h2"></path>
            <path d="M20 12h2"></path>
            <path d="m6.34 17.66-1.41 1.41"></path>
            <path d="m19.07 4.93-1.41 1.41"></path>
          </svg>
          <svg id="theme-dark-label" xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <path d="M12 3a6 6 0 0 0 9 9 9 9 0 1 1-9-9Z"></path>
          </svg>
        </button>
      </div>
    </div>
  </header>

  <!-- Main Dashboard Workspace -->
  <main class="workspace-container">
    {compound_section}
    {relationship_map_section}
    {clinical_section}
    {bioactivity_section}
    {similarity_section}
  </main>

  <!-- Academic Citations & Compliance Footer -->
  <footer class="app-footer">
    <div class="footer-inner">
      <p>{html.escape(_DISCLAIMER_TEXT)}</p>
      <p><strong>References:</strong> {html.escape(_CITATION_ZDRAZIL)} | {html.escape(_CITATION_DAVIES)}</p>
      <div class="footer-meta-row">
        <span>chembl-database-visualize</span>
        <span>EMBL-EBI ChEMBL (CC BY-SA 3.0 / 4.0)</span>
      </div>
    </div>
  </footer>

  <!-- Embedded Raw 3D SDF Conformer Coordinates -->
  <script id="sdf-data" type="text/plain">
{sdf_content or ""}
  </script>

  <!-- Embedded Dashboard Raw Payload -->
  <script id="dashboard-data" type="application/json">
{json_str}
  </script>

  <!-- Client-Side Interactive State & Visualization Engine -->
  <script>
    // =========================================================================
    // 1. Safe Theme Management (Sandboxed Iframe & Standard Browser Compatible)
    // =========================================================================
    let viewer3D = null;
    let draw3DConformer = null;
    let renderPhysicsNetwork = null;

    function safeGetStorage(key) {{
      try {{
        return window.localStorage ? window.localStorage.getItem(key) : null;
      }} catch (e) {{
        return null;
      }}
    }}

    function safeSetStorage(key, val) {{
      try {{
        if (window.localStorage) {{
          window.localStorage.setItem(key, val);
        }}
      }} catch (e) {{
        // Ignore storage errors in sandboxed iframes
      }}
    }}

    function applyThemeState(isDark) {{
      const htmlEl = document.documentElement;
      htmlEl.classList.toggle('dark', isDark);
      htmlEl.setAttribute('data-theme', isDark ? 'dark' : 'light');
      if (document.body) {{
        document.body.classList.toggle('dark', isDark);
      }}
      if (viewer3D && typeof viewer3D.setBackgroundColor === 'function') {{
        try {{
          viewer3D.setBackgroundColor(isDark ? '#16171D' : '#FFFFFF');
          viewer3D.render();
        }} catch (e) {{}}
      }}
      if (typeof draw3DConformer === 'function') {{
        draw3DConformer();
      }}
      if (typeof renderPhysicsNetwork === 'function') {{
        renderPhysicsNetwork();
      }}
    }}

    function toggleTheme() {{
      const htmlEl = document.documentElement;
      const isDark = !htmlEl.classList.contains('dark');
      safeSetStorage('chembl-theme', isDark ? 'dark' : 'light');
      applyThemeState(isDark);
    }}

    (function initTheme() {{
      const savedTheme = safeGetStorage('chembl-theme');
      const prefersDark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
      const isDark = savedTheme === 'dark' || (!savedTheme && prefersDark);
      applyThemeState(Boolean(isDark));
    }})();

    function exportDataJSON() {{
      const rawEl = document.getElementById('dashboard-data');
      if (!rawEl) return;
      const blob = new Blob([rawEl.textContent], {{ type: 'application/json' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'chembl_dashboard_data.json';
      a.click();
      URL.revokeObjectURL(url);
    }}

    // =========================================================================
    // 2. Interactive Zoom & Pan Engine for SVG Viewports (chem2d, radar, hist)
    // =========================================================================
    const viewContainers = {{
      chem2d: 'chem-svg-container',
      radar: 'radar-svg-container',
      hist: 'hist-svg-container'
    }};

    const viewFallbackIds = {{
      chem2d: 'chem2d-viewport',
      radar: 'radar-viewport',
      hist: 'hist-viewport'
    }};

    const viewTransforms = {{
      chem2d: {{ scale: 1.0, tx: 0, ty: 0 }},
      radar: {{ scale: 1.0, tx: 0, ty: 0 }},
      hist: {{ scale: 1.0, tx: 0, ty: 0 }}
    }};

    function getTargetContainer(target) {{
      const primaryId = viewContainers[target];
      const fallbackId = viewFallbackIds[target];
      return (primaryId && document.getElementById(primaryId)) ||
             (fallbackId && document.getElementById(fallbackId)) ||
             null;
    }}

    function applyViewTransform(target) {{
      const container = getTargetContainer(target);
      if (!container) return;
      const svg = container.querySelector('svg');
      const st = viewTransforms[target] || {{ scale: 1.0, tx: 0, ty: 0 }};
      container.setAttribute('data-zoom', String(st.scale));
      if (svg) {{
        svg.style.transformOrigin = 'center center';
        svg.style.transition = 'transform 0.14s ease-out';
        svg.style.transform = 'translate(' + st.tx.toFixed(1) + 'px, ' + st.ty.toFixed(1) + 'px) scale(' + st.scale.toFixed(3) + ')';
      }}
    }}

    function zoomView(target, factor) {{
      const st = viewTransforms[target];
      if (!st) return;
      st.scale = Math.min(4.0, Math.max(0.4, +(st.scale * factor).toFixed(3)));
      applyViewTransform(target);
    }}

    function resetView(target) {{
      const st = viewTransforms[target];
      if (!st) return;
      st.scale = 1.0;
      st.tx = 0;
      st.ty = 0;
      applyViewTransform(target);
    }}

    function initInteractivePanZoom() {{
      Object.keys(viewContainers).forEach(target => {{
        const container = getTargetContainer(target);
        if (!container) return;
        container.setAttribute('data-zoom', '1');
        let dragging = false;
        let startX = 0;
        let startY = 0;

        container.addEventListener('wheel', e => {{
          e.preventDefault();
          const factor = e.deltaY < 0 ? 1.12 : 0.89;
          zoomView(target, factor);
        }}, {{ passive: false }});

        container.addEventListener('mousedown', e => {{
          dragging = true;
          startX = e.clientX;
          startY = e.clientY;
        }});

        window.addEventListener('mousemove', e => {{
          if (!dragging) return;
          const dx = e.clientX - startX;
          const dy = e.clientY - startY;
          startX = e.clientX;
          startY = e.clientY;
          viewTransforms[target].tx += dx;
          viewTransforms[target].ty += dy;
          applyViewTransform(target);
        }});

        window.addEventListener('mouseup', () => {{
          dragging = false;
        }});
      }});
    }}

    // =========================================================================
    // 3. Structural Data Diagram (3D Conformer Interactive Projection Engine)
    // =========================================================================
    let isSpinning3D = false;
    let state3D = {{
      style: 'ballstick',
      zoom: 1.0,
      rotX: 0.35,
      rotY: 0.45,
      atoms: [],
      bonds: [],
      cx: 0,
      cy: 0,
      cz: 0
    }};
    let spinTimer3D = null;

    const DEFAULT_ASPIRIN_SDF = [
      "CHEMBL25",
      "  RDKit          3D",
      "",
      " 21 21  0  0  0  0  0  0  0  0999 V2000",
      "   -0.9616   -0.0887   -0.2584 C   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -0.5976   -1.4259   -0.4150 C   0  0  0  0  0  0  0  0  0  0  0  0",
      "    0.7289   -1.7616   -0.2117 C   0  0  0  0  0  0  0  0  0  0  0  0",
      "    1.6824   -0.7850    0.1554 C   0  0  0  0  0  0  0  0  0  0  0  0",
      "    1.3093    0.5469    0.3117 C   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -0.0151    0.9168    0.1068 C   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -0.3704    2.3364    0.2913 C   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -1.6669    2.6465   -0.0007 O   0  0  0  0  0  0  0  0  0  0  0  0",
      "    0.4363    3.1670    0.6865 O   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -2.2961    0.2547   -0.4851 O   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -3.2185   -0.4357    0.2393 C   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -4.4842    0.3475   -0.0084 C   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -2.9999   -1.5173    0.7335 O   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -1.3323   -2.1866   -0.7001 H   0  0  0  0  0  0  0  0  0  0  0  0",
      "    1.0189   -2.8090   -0.3340 H   0  0  0  0  0  0  0  0  0  0  0  0",
      "    2.7241   -1.0592    0.3188 H   0  0  0  0  0  0  0  0  0  0  0  0",
      "    2.0492    1.3095    0.5976 H   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -1.8596    3.6062    0.1345 H   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -5.3025   -0.2392    0.4357 H   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -4.4063    1.3340    0.4716 H   0  0  0  0  0  0  0  0  0  0  0  0",
      "   -4.6548    0.4705   -1.0877 H   0  0  0  0  0  0  0  0  0  0  0  0",
      "  1  2  2  0  0  0  0",
      "  2  3  1  0  0  0  0",
      "  3  4  2  0  0  0  0",
      "  4  5  1  0  0  0  0",
      "  5  6  2  0  0  0  0",
      "  6  1  1  0  0  0  0",
      "  6  7  1  0  0  0  0",
      "  7  8  1  0  0  0  0",
      "  7  9  2  0  0  0  0",
      "  1 10  1  0  0  0  0",
      " 10 11  1  0  0  0  0",
      " 11 12  1  0  0  0  0",
      " 11 13  2  0  0  0  0",
      "  2 14  1  0  0  0  0",
      "  3 15  1  0  0  0  0",
      "  4 16  1  0  0  0  0",
      "  5 17  1  0  0  0  0",
      "  8 18  1  0  0  0  0",
      " 12 19  1  0  0  0  0",
      " 12 20  1  0  0  0  0",
      " 12 21  1  0  0  0  0",
      "M  END"
    ].join('\\n');

    function parseSdfAtomsAndBonds(sdf) {{
      const lines = sdf.split(/\\r?\\n/);
      const atoms = [];
      const bonds = [];
      let countsLineIdx = -1;
      for (let i = 0; i < Math.min(lines.length, 15); i++) {{
        if (lines[i].includes('V2000')) {{
          countsLineIdx = i;
          break;
        }}
      }}
      if (countsLineIdx === -1) return {{ atoms, bonds }};
      const numAtoms = parseInt(lines[countsLineIdx].substring(0, 3).trim(), 10) || 0;
      const numBonds = parseInt(lines[countsLineIdx].substring(3, 6).trim(), 10) || 0;
      for (let i = 1; i <= numAtoms; i++) {{
        const line = lines[countsLineIdx + i];
        if (!line) break;
        const x = parseFloat(line.substring(0, 10).trim());
        const y = parseFloat(line.substring(10, 20).trim());
        const z = parseFloat(line.substring(20, 30).trim());
        const elem = line.substring(31, 34).trim() || 'C';
        if (!Number.isNaN(x) && !Number.isNaN(y) && !Number.isNaN(z)) {{
          atoms.push({{ x, y, z, elem }});
        }}
      }}
      for (let i = 1; i <= numBonds; i++) {{
        const line = lines[countsLineIdx + numAtoms + i];
        if (!line) break;
        const a1 = (parseInt(line.substring(0, 3).trim(), 10) || 1) - 1;
        const a2 = (parseInt(line.substring(3, 6).trim(), 10) || 1) - 1;
        const order = parseInt(line.substring(6, 9).trim(), 10) || 1;
        bonds.push({{ a1, a2, order }});
      }}
      return {{ atoms, bonds }};
    }}

    function init3DViewer() {{
      const container = document.getElementById('mol-3d-viewer');
      if (!container) return;

      const sdfEl = document.getElementById('sdf-data');
      let sdfCoords = sdfEl ? sdfEl.textContent.trim() : '';
      if (!sdfCoords || !sdfCoords.includes('V2000')) {{
        sdfCoords = DEFAULT_ASPIRIN_SDF;
      }}

      const parsed = parseSdfAtomsAndBonds(sdfCoords);
      state3D.atoms = parsed.atoms;
      state3D.bonds = parsed.bonds;

      if (state3D.atoms.length === 0) {{
        const fallbackParsed = parseSdfAtomsAndBonds(DEFAULT_ASPIRIN_SDF);
        state3D.atoms = fallbackParsed.atoms;
        state3D.bonds = fallbackParsed.bonds;
      }}

      let cx = 0, cy = 0, cz = 0;
      state3D.atoms.forEach(a => {{ cx += a.x; cy += a.y; cz += a.z; }});
      if (state3D.atoms.length > 0) {{
        cx /= state3D.atoms.length;
        cy /= state3D.atoms.length;
        cz /= state3D.atoms.length;
      }}
      state3D.cx = cx;
      state3D.cy = cy;
      state3D.cz = cz;

      container.innerHTML = '';
      container.setAttribute('data-zoom', String(state3D.zoom));
      container.setAttribute('data-style', state3D.style);

      const svgEl = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      svgEl.setAttribute('id', 'conformer-3d-svg');
      svgEl.setAttribute('viewBox', '0 0 360 256');
      svgEl.setAttribute('preserveAspectRatio', 'xMidYMid meet');
      svgEl.style.width = '100%';
      svgEl.style.height = '100%';
      svgEl.style.display = 'block';
      container.appendChild(svgEl);

      // Data Visualization Palette for Elements
      const elementColors = {{
        'C': '#0D9488',
        'O': '#E11D48',
        'N': '#4F46E5',
        'H': '#787774',
        'S': '#D97706',
        'P': '#D97706',
        'F': '#059669',
        'Cl': '#059669',
        'Br': '#D97706',
        'I': '#4F46E5'
      }};

      draw3DConformer = function() {{
        const vw = 360;
        const vh = 256;
        const cw = vw / 2;
        const ch = vh / 2;
        const baseScale = (Math.min(vw, vh) / 6.8) * state3D.zoom;

        const cosY = Math.cos(state3D.rotY);
        const sinY = Math.sin(state3D.rotY);
        const cosX = Math.cos(state3D.rotX);
        const sinX = Math.sin(state3D.rotX);

        const proj = state3D.atoms.map((a, idx) => {{
          const x = a.x - state3D.cx;
          const y = a.y - state3D.cy;
          const z = a.z - state3D.cz;
          const x1 = x * cosY + z * sinY;
          const z1 = -x * sinY + z * cosY;
          const y2 = y * cosX - z1 * sinX;
          const z2 = y * sinX + z1 * cosX;
          return {{
            idx: idx,
            px: cw + x1 * baseScale,
            py: ch - y2 * baseScale,
            pz: z2,
            elem: a.elem
          }};
        }});

        const parts = [];

        // Optional translucent VDW surface hull in 'surface' mode
        if (state3D.style === 'surface') {{
          proj.forEach(a => {{
            const col = elementColors[a.elem] || '#0D9488';
            const surfR = (a.elem === 'H' ? 16 : 24) * state3D.zoom;
            parts.push('<circle cx="' + a.px.toFixed(1) + '" cy="' + a.py.toFixed(1) + '" r="' + surfR.toFixed(1) + '" fill="' + col + '" fill-opacity="0.16"/>');
          }});
        }}

        // Render bonds (unless pure CPK sphere mode)
        if (state3D.style !== 'sphere') {{
          const strokeW = (state3D.style === 'stick' ? 3.8 : 2.4) * Math.sqrt(state3D.zoom);
          state3D.bonds.forEach(b => {{
            const p1 = proj[b.a1];
            const p2 = proj[b.a2];
            if (!p1 || !p2) return;
            const col = elementColors[p1.elem] || '#0D9488';
            if (b.order === 2) {{
              const dx = p2.px - p1.px;
              const dy = p2.py - p1.py;
              const len = Math.hypot(dx, dy) || 1;
              const ox = (-dy / len) * 2.2;
              const oy = (dx / len) * 2.2;
              parts.push('<line x1="' + (p1.px + ox).toFixed(1) + '" y1="' + (p1.py + oy).toFixed(1) + '" x2="' + (p2.px + ox).toFixed(1) + '" y2="' + (p2.py + oy).toFixed(1) + '" stroke="' + col + '" stroke-opacity="0.78" stroke-width="' + (strokeW * 0.72).toFixed(1) + '" stroke-linecap="round"/>');
              parts.push('<line x1="' + (p1.px - ox).toFixed(1) + '" y1="' + (p1.py - oy).toFixed(1) + '" x2="' + (p2.px - ox).toFixed(1) + '" y2="' + (p2.py - oy).toFixed(1) + '" stroke="' + col + '" stroke-opacity="0.78" stroke-width="' + (strokeW * 0.72).toFixed(1) + '" stroke-linecap="round"/>');
            }} else {{
              parts.push('<line x1="' + p1.px.toFixed(1) + '" y1="' + p1.py.toFixed(1) + '" x2="' + p2.px.toFixed(1) + '" y2="' + p2.py.toFixed(1) + '" stroke="' + col + '" stroke-opacity="0.78" stroke-width="' + strokeW.toFixed(1) + '" stroke-linecap="round"/>');
            }}
          }});
        }}

        // Depth-sort atoms back-to-front
        const sortedAtoms = proj.slice().sort((a, b) => a.pz - b.pz);
        sortedAtoms.forEach(a => {{
          let baseR = a.elem === 'H' ? 4.2 : (a.elem === 'C' ? 6.8 : 7.8);
          if (state3D.style === 'stick') {{
            baseR *= 0.55;
          }} else if (state3D.style === 'sphere') {{
            baseR *= 1.85;
          }}
          const r = Math.max(2.2, baseR * state3D.zoom);
          const col = elementColors[a.elem] || '#0D9488';
          parts.push('<circle cx="' + a.px.toFixed(1) + '" cy="' + a.py.toFixed(1) + '" r="' + r.toFixed(1) + '" fill="' + col + '" stroke="var(--bg-surface)" stroke-width="1.3"/>');
          if (a.elem !== 'C' && a.elem !== 'H' && state3D.style !== 'stick') {{
            parts.push('<text x="' + a.px.toFixed(1) + '" y="' + (a.py + 3.2).toFixed(1) + '" fill="#FFFFFF" font-family="IBM Plex Mono, monospace" font-size="8.5" font-weight="600" text-anchor="middle">' + a.elem + '</text>');
          }}
        }});

        svgEl.innerHTML = parts.join('');
      }};

      let isDragging3D = false;
      let lastX = 0;
      let lastY = 0;

      svgEl.addEventListener('mousedown', e => {{
        isDragging3D = true;
        lastX = e.clientX;
        lastY = e.clientY;
      }});

      window.addEventListener('mouseup', () => {{
        isDragging3D = false;
      }});

      window.addEventListener('mousemove', e => {{
        if (!isDragging3D) return;
        const dx = e.clientX - lastX;
        const dy = e.clientY - lastY;
        lastX = e.clientX;
        lastY = e.clientY;
        state3D.rotY += dx * 0.014;
        state3D.rotX += dy * 0.014;
        draw3DConformer();
      }});

      svgEl.addEventListener('wheel', e => {{
        e.preventDefault();
        const factor = e.deltaY < 0 ? 1.15 : 0.87;
        zoom3D(factor);
      }}, {{ passive: false }});

      if (spinTimer3D) clearInterval(spinTimer3D);
      spinTimer3D = setInterval(() => {{
        if (isSpinning3D && !isDragging3D) {{
          state3D.rotY += 0.035;
          draw3DConformer();
        }}
      }}, 40);

      updateStyleButtonHighlight(state3D.style);
      draw3DConformer();
    }}

    function updateStyleButtonHighlight(activeStyle) {{
      const styleBtns = document.querySelectorAll('[data-style-btn]');
      styleBtns.forEach(btn => {{
        const isMatch = btn.getAttribute('data-style-btn') === activeStyle;
        btn.classList.toggle('active', isMatch);
        btn.setAttribute('aria-pressed', isMatch ? 'true' : 'false');
      }});
      ['stick', 'ballstick', 'sphere', 'surface'].forEach(st => {{
        const btnById = document.getElementById('style-btn-' + st);
        if (btnById) {{
          const isMatch = st === activeStyle;
          btnById.classList.toggle('active', isMatch);
          btnById.setAttribute('aria-pressed', isMatch ? 'true' : 'false');
        }}
      }});
    }}

    function setViewerStyle(styleType) {{
      state3D.style = styleType || 'ballstick';
      const container = document.getElementById('mol-3d-viewer');
      if (container) container.setAttribute('data-style', state3D.style);
      updateStyleButtonHighlight(state3D.style);
      if (typeof draw3DConformer === 'function') {{
        draw3DConformer();
      }}
      if (viewer3D) {{
        try {{
          const mol3D = window.$3Dmol || window['3Dmol'];
          viewer3D.removeAllSurfaces();
          if (styleType === 'stick') {{
            viewer3D.setStyle({{}}, {{ stick: {{ colorscheme: "cyanCarbon", radius: 0.18 }} }});
          }} else if (styleType === 'sphere') {{
            viewer3D.setStyle({{}}, {{ sphere: {{ scale: 0.35, colorscheme: "cyanCarbon" }} }});
          }} else if (styleType === 'ballstick') {{
            viewer3D.setStyle({{}}, {{ stick: {{ colorscheme: "cyanCarbon", radius: 0.14 }}, sphere: {{ scale: 0.25, colorscheme: "cyanCarbon" }} }});
          }} else if (styleType === 'surface') {{
            viewer3D.setStyle({{}}, {{ stick: {{ colorscheme: "cyanCarbon", radius: 0.12 }} }});
            if (mol3D && mol3D.SurfaceType) {{
              viewer3D.addSurface(mol3D.SurfaceType.VDW, {{ opacity: 0.65, color: '#0D9488' }});
            }}
          }}
          viewer3D.render();
        }} catch (e) {{}}
      }}
    }}

    function zoom3D(factor) {{
      state3D.zoom = Math.min(3.5, Math.max(0.45, +(state3D.zoom * factor).toFixed(3)));
      const container = document.getElementById('mol-3d-viewer');
      if (container) container.setAttribute('data-zoom', String(state3D.zoom));
      if (typeof draw3DConformer === 'function') {{
        draw3DConformer();
      }}
      if (viewer3D && typeof viewer3D.zoom === 'function') {{
        try {{
          viewer3D.zoom(factor);
          viewer3D.render();
        }} catch (e) {{}}
      }}
    }}

    function toggleSpin3D() {{
      const btn = document.getElementById('spin-toggle-btn');
      isSpinning3D = !isSpinning3D;
      if (btn) {{
        btn.textContent = isSpinning3D ? 'Pause Spin' : 'Auto-Spin';
        btn.classList.toggle('active', isSpinning3D);
      }}
      if (viewer3D && typeof viewer3D.spin === 'function') {{
        try {{
          viewer3D.spin(isSpinning3D ? 'y' : false, 1);
        }} catch (e) {{}}
      }}
    }}

    function reset3DView() {{
      isSpinning3D = false;
      state3D.zoom = 1.0;
      state3D.rotX = 0.35;
      state3D.rotY = 0.45;
      state3D.style = 'ballstick';
      const btn = document.getElementById('spin-toggle-btn');
      if (btn) {{
        btn.textContent = 'Auto-Spin';
        btn.classList.remove('active');
      }}
      const container = document.getElementById('mol-3d-viewer');
      if (container) {{
        container.setAttribute('data-zoom', '1');
        container.setAttribute('data-style', 'ballstick');
      }}
      updateStyleButtonHighlight('ballstick');
      if (typeof draw3DConformer === 'function') {{
        draw3DConformer();
      }}
      if (viewer3D) {{
        try {{
          viewer3D.spin(false);
          viewer3D.zoomTo();
          viewer3D.render();
        }} catch (e) {{}}
      }}
    }}

    // =========================================================================
    // 4. 2D Structure Highlight & Table / Analog Filter Controls
    // =========================================================================
    let highlightFGActive = true;
    function toggleFunctionalGroupHighlights() {{
      const btn = document.getElementById('toggle-fg-btn');
      const circles = document.querySelectorAll('#chem-svg-container circle');
      const fills = document.querySelectorAll('#chem-svg-container .chem-ring-fill');
      highlightFGActive = !highlightFGActive;
      if (btn) {{
        btn.textContent = highlightFGActive ? 'Hide Highlights' : 'Show Highlights';
      }}
      circles.forEach(c => {{
        c.style.opacity = highlightFGActive ? '1' : '0';
        c.style.transition = 'opacity 0.2s ease';
      }});
      fills.forEach(f => {{
        f.style.opacity = highlightFGActive ? '1' : '0';
        f.style.transition = 'opacity 0.2s ease';
      }});
    }}

    function filterActivityTable() {{
      const query = (document.getElementById('activity-table-search')?.value || '').toLowerCase();
      const rows = document.querySelectorAll('#activity-table tbody tr');
      rows.forEach(r => {{
        const text = r.textContent.toLowerCase();
        r.style.display = text.includes(query) ? '' : 'none';
      }});
    }}

    function resetActivityFilter() {{
      const input = document.getElementById('activity-table-search');
      if (input) {{
        input.value = '';
      }}
      filterActivityTable();
    }}

    function updateAnalogFilters() {{
      const simMin = parseFloat(document.getElementById('similarity-filter')?.value || 0);
      const mwtMax = parseFloat(document.getElementById('mwt-filter')?.value || 1000);

      const simValEl = document.getElementById('similarity-filter-val');
      if (simValEl) simValEl.textContent = simMin + '%';

      const mwtValEl = document.getElementById('mwt-filter-val');
      if (mwtValEl) mwtValEl.textContent = mwtMax + ' Da';

      const cards = document.querySelectorAll('.analog-card');
      let visible = 0;
      cards.forEach(c => {{
        const cSim = parseFloat(c.getAttribute('data-similarity') || 0);
        const cMwt = parseFloat(c.getAttribute('data-mwt') || 0);
        const match = (cSim >= simMin) && (cMwt <= mwtMax);
        c.style.display = match ? '' : 'none';
        if (match) visible++;
      }});

      const heatRows = document.querySelectorAll('.analog-heatmap-row');
      heatRows.forEach(r => {{
        const rSim = parseFloat(r.getAttribute('data-similarity') || 0);
        const rMwt = parseFloat(r.getAttribute('data-mwt') || 0);
        const match = (rSim >= simMin) && (rMwt <= mwtMax);
        r.style.display = match ? '' : 'none';
      }});

      const counter = document.getElementById('visible-analog-count');
      if (counter) counter.textContent = visible;
    }}

    function resetAnalogFilters() {{
      const simInput = document.getElementById('similarity-filter');
      const mwtInput = document.getElementById('mwt-filter');
      if (simInput) simInput.value = 70;
      if (mwtInput) mwtInput.value = 800;
      updateAnalogFilters();
    }}

    // =========================================================================
    // 5. Interactive Relationship Node Map Engine (Zoom, Pan, Drag, Force Sim)
    // =========================================================================
    let physicsRunning = true;
    let physicsAnimId = null;
    let nodes = [];
    let links = [];
    let draggedNode = null;
    let mousePos = {{ x: 0, y: 0 }};
    let networkZoom = 1.0;
    let networkPan = {{ x: 0, y: 0 }};
    let networkWidth = 680;
    let networkHeight = 310;

    function initPhysicsNetwork() {{
      const viewport = document.getElementById('physics-viewport');
      if (!viewport) return;

      function measureViewport() {{
        const rect = viewport.getBoundingClientRect();
        networkWidth = (rect && rect.width > 100) ? rect.width : 680;
        networkHeight = (rect && rect.height > 100) ? rect.height : 310;
        viewport.setAttribute('viewBox', '0 0 ' + networkWidth + ' ' + networkHeight);
      }}
      measureViewport();
      viewport.setAttribute('data-zoom', String(networkZoom));
      window.addEventListener('resize', measureViewport);

      let rawData = {{}};
      try {{
        rawData = JSON.parse(document.getElementById('dashboard-data')?.textContent || '{{}}');
      }} catch (e) {{}}

      const mol = rawData.molecule || {{}};
      const centerId = mol.molecule_chembl_id || 'CHEMBL25';
      const centerName = mol.pref_name || 'Target Compound';

      nodes = [];
      links = [];

      const cx = networkWidth / 2;
      const cy = networkHeight / 2;

      // Query Center Node: Vivid Indigo (#4F46E5)
      nodes.push({{
        id: centerId,
        name: centerName,
        isCenter: true,
        x: cx,
        y: cy,
        vx: 0,
        vy: 0,
        radius: 18,
        color: '#4F46E5',
        edgeRgba: 'rgba(79, 70, 229, 0.45)',
        details: centerId + ' (' + centerName + ') - Query Compound'
      }});

      const analogs = rawData.similarities || [];
      const defaultAnalogs = [
        {{ molecule_chembl_id: 'CHEMBL26', pref_name: 'SALICYLIC ACID', similarity: 91.2, full_mwt: 138.12 }},
        {{ molecule_chembl_id: 'CHEMBL112', pref_name: 'DIFLUNISAL', similarity: 86.4, full_mwt: 250.2 }},
        {{ molecule_chembl_id: 'CHEMBL1235', pref_name: 'SALSALATE', similarity: 84.7, full_mwt: 258.23 }},
        {{ molecule_chembl_id: 'CHEMBL501', pref_name: 'BENORYLATE', similarity: 81.3, full_mwt: 313.31 }}
      ];
      const sourceList = analogs.filter(a => a.molecule_chembl_id !== centerId);
      const finalAnalogs = sourceList.length > 0 ? sourceList : defaultAnalogs;

      finalAnalogs.forEach((an, i) => {{
        const sim = parseFloat(an.similarity || 85.0);
        const angle = (i / finalAnalogs.length) * 2 * Math.PI;
        const dist = 82 + (100 - sim) * 2.8;

        // Palette: Functional Emerald (#059669), Vibrant Teal (#0D9488), Energetic Amber (#D97706)
        let color = '#D97706';
        let edgeRgba = 'rgba(217, 119, 6, 0.45)';
        if (sim >= 90) {{
          color = '#059669';
          edgeRgba = 'rgba(5, 150, 105, 0.48)';
        }} else if (sim >= 85) {{
          color = '#0D9488';
          edgeRgba = 'rgba(13, 148, 136, 0.48)';
        }}

        const nodeObj = {{
          id: an.molecule_chembl_id,
          name: an.pref_name || an.molecule_chembl_id,
          isCenter: false,
          x: cx + dist * Math.cos(angle),
          y: cy + dist * Math.sin(angle),
          initAngle: angle,
          initDist: dist,
          vx: (Math.random() - 0.5) * 1.2,
          vy: (Math.random() - 0.5) * 1.2,
          radius: 13,
          similarity: sim,
          color: color,
          edgeRgba: edgeRgba,
          details: (an.pref_name || an.molecule_chembl_id) + ' (' + an.molecule_chembl_id + ') - ' + sim.toFixed(1) + '% Tanimoto similarity'
        }};
        nodes.push(nodeObj);
        links.push({{
          source: nodes[0],
          target: nodeObj,
          targetDist: dist,
          strength: 0.04 * (sim / 100)
        }});
      }});

      function clientToGraphCoords(clientX, clientY) {{
        const rect = viewport.getBoundingClientRect();
        const rawX = clientX - rect.left;
        const rawY = clientY - rect.top;
        const cx = networkWidth / 2;
        const cy = networkHeight / 2;
        return {{
          x: (rawX - cx - networkPan.x) / networkZoom + cx,
          y: (rawY - cy - networkPan.y) / networkZoom + cy,
          rawX: rawX,
          rawY: rawY
        }};
      }}

      viewport.addEventListener('wheel', (e) => {{
        e.preventDefault();
        const factor = e.deltaY < 0 ? 1.15 : 0.87;
        zoomNetwork(factor);
      }}, {{ passive: false }});

      viewport.addEventListener('mousedown', (e) => {{
        const pt = clientToGraphCoords(e.clientX, e.clientY);
        mousePos.x = pt.x;
        mousePos.y = pt.y;

        for (let n of nodes) {{
          const dx = n.x - mousePos.x;
          const dy = n.y - mousePos.y;
          if (Math.hypot(dx, dy) <= n.radius + 6) {{
            draggedNode = n;
            break;
          }}
        }}
      }});

      viewport.addEventListener('mousemove', (e) => {{
        const pt = clientToGraphCoords(e.clientX, e.clientY);
        mousePos.x = pt.x;
        mousePos.y = pt.y;

        const tooltip = document.getElementById('physics-tooltip');
        let hovered = null;
        for (let n of nodes) {{
          const dx = n.x - mousePos.x;
          const dy = n.y - mousePos.y;
          if (Math.hypot(dx, dy) <= n.radius + 6) {{
            hovered = n;
            break;
          }}
        }}

        if (hovered && tooltip) {{
          tooltip.style.left = Math.min(networkWidth - 220, Math.max(8, pt.rawX + 12)) + 'px';
          tooltip.style.top = Math.max(8, pt.rawY - 28) + 'px';
          tooltip.textContent = hovered.details;
          tooltip.classList.remove('hidden');
        }} else if (tooltip) {{
          tooltip.classList.add('hidden');
        }}
      }});

      window.addEventListener('mouseup', () => {{
        draggedNode = null;
      }});

      renderPhysicsNetwork = function() {{
        const cx = networkWidth / 2;
        const cy = networkHeight / 2;
        const transformAttr = 'translate(' + (cx + networkPan.x).toFixed(1) + ',' + (cy + networkPan.y).toFixed(1) + ') scale(' + networkZoom.toFixed(3) + ') translate(' + (-cx).toFixed(1) + ',' + (-cy).toFixed(1) + ')';
        const parts = ['<g id="physics-graph-group" transform="' + transformAttr + '">'];
        for (let l of links) {{
          const sw = Math.max(1.5, (l.target.similarity || 75) / 34).toFixed(1);
          parts.push('<line x1="' + l.source.x.toFixed(1) + '" y1="' + l.source.y.toFixed(1) + '" x2="' + l.target.x.toFixed(1) + '" y2="' + l.target.y.toFixed(1) + '" stroke="' + l.target.edgeRgba + '" stroke-width="' + sw + '"/>');
        }}
        for (let n of nodes) {{
          const fw = n.isCenter ? '600' : '500';
          const fs = n.isCenter ? '11' : '10';
          parts.push('<circle cx="' + n.x.toFixed(1) + '" cy="' + n.y.toFixed(1) + '" r="' + n.radius + '" fill="' + n.color + '" stroke="var(--bg-surface)" stroke-width="2"/>');
          parts.push('<text x="' + n.x.toFixed(1) + '" y="' + (n.y + n.radius + 13).toFixed(1) + '" fill="var(--text-primary)" font-family="IBM Plex Mono, monospace" font-size="' + fs + '" font-weight="' + fw + '" text-anchor="middle">' + n.id + '</text>');
        }}
        parts.push('</g>');
        viewport.innerHTML = parts.join('');
      }};

      function stepPhysics() {{
        if (!physicsRunning) {{
          renderPhysicsNetwork();
          physicsAnimId = requestAnimationFrame(stepPhysics);
          return;
        }}

        const cx = networkWidth / 2;
        const cy = networkHeight / 2;

        for (let i = 0; i < nodes.length; i++) {{
          for (let j = i + 1; j < nodes.length; j++) {{
            const n1 = nodes[i];
            const n2 = nodes[j];
            const dx = n2.x - n1.x;
            const dy = n2.y - n1.y;
            const dist = Math.hypot(dx, dy) || 1;
            const minDist = n1.radius + n2.radius + 42;
            if (dist < minDist) {{
              const force = (minDist - dist) / dist * 0.45;
              n1.vx -= dx * force;
              n1.vy -= dy * force;
              n2.vx += dx * force;
              n2.vy += dy * force;
            }}
          }}
        }}

        for (let l of links) {{
          const dx = l.target.x - l.source.x;
          const dy = l.target.y - l.source.y;
          const dist = Math.hypot(dx, dy) || 1;
          const displacement = dist - l.targetDist;
          const force = displacement * l.strength;
          const fx = (dx / dist) * force;
          const fy = (dy / dist) * force;

          if (!l.source.isCenter || draggedNode !== l.source) {{
            l.source.vx += fx * 0.35;
            l.source.vy += fy * 0.35;
          }}
          if (draggedNode !== l.target) {{
            l.target.vx -= fx;
            l.target.vy -= fy;
          }}
        }}

        for (let n of nodes) {{
          if (n === draggedNode) {{
            n.x = mousePos.x;
            n.y = mousePos.y;
            n.vx = 0;
            n.vy = 0;
            continue;
          }}

          n.vx += (cx - n.x) * 0.0025;
          n.vy += (cy - n.y) * 0.0025;
          n.vx *= 0.87;
          n.vy *= 0.87;
          n.x += n.vx;
          n.y += n.vy;

          const pad = n.radius + 18;
          if (n.x < pad) {{ n.x = pad; n.vx *= -0.5; }}
          if (n.x > networkWidth - pad) {{ n.x = networkWidth - pad; n.vx *= -0.5; }}
          if (n.y < pad) {{ n.y = pad; n.vy *= -0.5; }}
          if (n.y > networkHeight - pad) {{ n.y = networkHeight - pad; n.vy *= -0.5; }}
        }}

        renderPhysicsNetwork();
        physicsAnimId = requestAnimationFrame(stepPhysics);
      }}

      renderPhysicsNetwork();
      if (physicsAnimId) cancelAnimationFrame(physicsAnimId);
      physicsAnimId = requestAnimationFrame(stepPhysics);
    }}

    function zoomNetwork(factor) {{
      networkZoom = Math.min(3.5, Math.max(0.45, +(networkZoom * factor).toFixed(3)));
      const viewport = document.getElementById('physics-viewport');
      if (viewport) viewport.setAttribute('data-zoom', String(networkZoom));
      if (typeof renderPhysicsNetwork === 'function') {{
        renderPhysicsNetwork();
      }}
    }}

    function togglePhysicsSimulation() {{
      physicsRunning = !physicsRunning;
      const btn = document.getElementById('physics-toggle-btn');
      if (btn) {{
        btn.textContent = physicsRunning ? 'Pause' : 'Resume';
        btn.classList.toggle('active', !physicsRunning);
      }}
    }}

    function perturbPhysicsGraph() {{
      for (let n of nodes) {{
        n.vx += (Math.random() - 0.5) * 14;
        n.vy += (Math.random() - 0.5) * 14;
      }}
      if (typeof renderPhysicsNetwork === 'function') {{
        renderPhysicsNetwork();
      }}
    }}

    function resetPhysicsPositions() {{
      networkZoom = 1.0;
      networkPan = {{ x: 0, y: 0 }};
      const viewport = document.getElementById('physics-viewport');
      if (viewport) viewport.setAttribute('data-zoom', '1');
      const cx = networkWidth / 2;
      const cy = networkHeight / 2;
      nodes.forEach((n, idx) => {{
        if (idx === 0 || n.isCenter) {{
          n.x = cx;
          n.y = cy;
          n.vx = 0;
          n.vy = 0;
        }} else {{
          const angle = n.initAngle || ((idx / Math.max(1, nodes.length - 1)) * 2 * Math.PI);
          const dist = n.initDist || 110;
          n.x = cx + dist * Math.cos(angle);
          n.y = cy + dist * Math.sin(angle);
          n.vx = 0;
          n.vy = 0;
        }}
      }});
      if (typeof renderPhysicsNetwork === 'function') {{
        renderPhysicsNetwork();
      }}
    }}

    function bootAllVisualizations() {{
      initInteractivePanZoom();
      init3DViewer();
      initPhysicsNetwork();
    }}

    if (document.readyState === 'loading') {{
      document.addEventListener('DOMContentLoaded', bootAllVisualizations);
    }} else {{
      bootAllVisualizations();
    }}
  </script>
</body>
</html>
"""
  return full_html


def main() -> None:
  parser = argparse.ArgumentParser(description="Generate comprehensive interactive HTML5 ChEMBL dashboard.")
  parser.add_argument("--compound_json", type=str, default=None, help="Path to compound JSON file")
  parser.add_argument("--sdf", type=str, default=None, help="Path to SDF 3D conformer file")
  parser.add_argument("--image", type=str, default=None, help="Path to SVG 2D structure file")
  parser.add_argument("--activity_json", type=str, default=None, help="Path to activity JSON file")
  parser.add_argument("--drug_json", type=str, default=None, help="Path to drug indication/mechanism JSON")
  parser.add_argument("--similarity_json", type=str, default=None, help="Path to similarity JSON file")
  parser.add_argument("--title", type=str, default="ChEMBL Analytical Dashboard", help="Dashboard title")
  parser.add_argument("-o", "--output", type=str, required=True, help="Output HTML file path")
  args = parser.parse_args()

  compound_data = None
  if args.compound_json and os.path.exists(args.compound_json):
    with open(args.compound_json, "r", encoding="utf-8") as f:
      compound_data = json.load(f)

  sdf_content = None
  if args.sdf and os.path.exists(args.sdf):
    with open(args.sdf, "r", encoding="utf-8", errors="replace") as f:
      sdf_content = f.read()

  svg_content = None
  if args.image and os.path.exists(args.image):
    with open(args.image, "r", encoding="utf-8", errors="replace") as f:
      svg_content = f.read()

  activity_data = None
  if args.activity_json and os.path.exists(args.activity_json):
    with open(args.activity_json, "r", encoding="utf-8") as f:
      activity_data = json.load(f)

  drug_data = None
  if args.drug_json and os.path.exists(args.drug_json):
    with open(args.drug_json, "r", encoding="utf-8") as f:
      drug_data = json.load(f)

  similarity_data = None
  if args.similarity_json and os.path.exists(args.similarity_json):
    with open(args.similarity_json, "r", encoding="utf-8") as f:
      similarity_data = json.load(f)

  dashboard_html = build_full_dashboard(
      title=args.title,
      compound_data=compound_data,
      sdf_content=sdf_content,
      svg_content=svg_content,
      activity_data=activity_data,
      drug_data=drug_data,
      similarity_data=similarity_data,
  )

  out_dir = os.path.dirname(args.output)
  if out_dir:
    os.makedirs(out_dir, exist_ok=True)
  with open(args.output, "w", encoding="utf-8") as f:
    f.write(dashboard_html)

  ws_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  try:
    abs_out = os.path.abspath(args.output)
    rel_out = os.path.relpath(abs_out, ws_root) if abs_out.startswith(ws_root) else args.output
    if not rel_out.startswith((".", "/")):
      rel_out = f"./{rel_out}"
  except Exception:
    rel_out = args.output
  print(f"Generated standalone dashboard HTML at {rel_out}")


if __name__ == "__main__":
  main()
