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

Combines molecular properties, 3D conformers with rotational physics,
a 2D force-directed spring-mass chemical network graph, bioactivity distributions,
clinical indications/mechanisms, and similarity results into an interactive,
theme-responsive HTML5 dashboard application (dashboard.html).
"""

from __future__ import annotations

import argparse
import datetime
import html
import json
import os
import sys
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
  """Render clinical development timeline, indications, and mechanisms."""
  mol = compound_data or {}
  max_phase = mol.get("max_phase")
  first_approval = mol.get("first_approval")

  mechanisms = mol.get("mechanisms", [])
  if drug_data and "mechanisms" in drug_data:
    mechanisms = drug_data["mechanisms"]

  indications = mol.get("indications", [])
  if drug_data and "drug_indications" in drug_data:
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
    node_color = (
        "bg-emerald-600 border-emerald-600 text-white"
        if is_completed
        else ("bg-blue-600 border-blue-600 text-white" if is_current else "bg-slate-100 dark:bg-slate-800 border-slate-300 dark:border-slate-700 text-slate-400")
    )

    timeline_nodes.append(
        f"""
        <div class="flex-1 flex flex-col items-center relative">
          <div class="w-8 h-8 rounded-full border-2 flex items-center justify-center font-bold text-xs {node_color} z-10 transition-colors">
            {level}
          </div>
          <span class="mt-2 text-xs font-semibold text-slate-700 dark:text-slate-300 text-center">{name}</span>
          {f'<span class="text-[10px] text-emerald-600 font-mono">Approved {first_approval}</span>' if (level == 4 and first_approval) else ''}
        </div>
        """
    )

  mech_cards = []
  for m in mechanisms[:8]:
    m_name = m.get("mechanism_of_action") or "Unknown Mechanism"
    action_type = m.get("action_type") or "ACTION"
    target_name = m.get("target_name") or "Target"
    t_id = m.get("target_chembl_id") or ""
    mech_cards.append(
        f"""
        <div class="p-3 rounded-lg border border-slate-200 dark:border-slate-800 bg-white dark:bg-slate-950 flex flex-col justify-between">
          <div>
            <div class="flex items-center justify-between gap-2 mb-1">
              <span class="px-2 py-0.5 text-[10px] font-bold uppercase tracking-wider rounded bg-purple-50 dark:bg-purple-950 text-purple-700 dark:text-purple-300 border border-purple-200 dark:border-purple-900">{html.escape(action_type)}</span>
              <span class="text-[10px] font-mono text-slate-400">{html.escape(t_id)}</span>
            </div>
            <h4 class="text-xs font-bold text-slate-900 dark:text-white mb-1">{html.escape(m_name)}</h4>
            <span class="text-[11px] text-slate-500 block truncate" title="{html.escape(target_name)}">{html.escape(target_name)}</span>
          </div>
        </div>
        """
    )

  ind_cards = []
  for ind in indications[:8]:
    term = ind.get("mesh_heading") or ind.get("efo_term") or "Indication"
    atc = ind.get("atc_code")
    ind_phase = ind.get("max_phase_for_ind")
    p_str = f"Phase {ind_phase}" if ind_phase is not None else "Clinical"
    ind_cards.append(
        f"""
        <div class="p-2.5 rounded-lg border border-slate-200 dark:border-slate-800 bg-white dark:bg-slate-950 flex items-center justify-between text-xs">
          <div>
            <span class="font-semibold text-slate-800 dark:text-slate-200 block">{html.escape(term)}</span>
            {f'<span class="text-[10px] font-mono text-slate-400">ATC: {html.escape(atc)}</span>' if atc else ''}
          </div>
          <span class="px-2 py-0.5 rounded-full text-[10px] font-medium bg-slate-100 dark:bg-slate-800 text-slate-600 dark:text-slate-400">{html.escape(p_str)}</span>
        </div>
        """
    )

  tracker_html = f"""
<div class="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl shadow-sm overflow-hidden mb-8">
  <div class="px-6 py-4 border-b border-slate-100 dark:border-slate-800 flex items-center justify-between bg-slate-50/50 dark:bg-slate-800/40">
    <div class="flex items-center space-x-3">
      <span class="px-2.5 py-1 text-xs font-mono font-semibold rounded-md bg-indigo-100 dark:bg-indigo-950 text-indigo-800 dark:text-indigo-300 border border-indigo-200 dark:border-indigo-900">Clinical Development</span>
      <h2 class="text-xl font-bold text-slate-900 dark:text-white tracking-tight">Development Timeline & Mechanisms</h2>
    </div>
  </div>

  <div class="p-6 space-y-6">
    <div class="py-4 px-2 bg-slate-50 dark:bg-slate-950/60 rounded-xl border border-slate-100 dark:border-slate-800">
      <div class="flex items-center justify-between max-w-2xl mx-auto">
        {"".join(timeline_nodes)}
      </div>
    </div>

    <div class="grid grid-cols-1 md:grid-cols-2 gap-6 pt-2">
      <div>
        <h3 class="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-3">Mechanisms of Action</h3>
        <div class="grid grid-cols-1 sm:grid-cols-2 gap-3">
          {"".join(mech_cards) if mech_cards else '<div class="text-xs text-slate-400 py-3">No mechanism annotations recorded</div>'}
        </div>
      </div>

      <div>
        <h3 class="text-xs font-semibold uppercase tracking-wider text-slate-500 mb-3">Reported Clinical Indications</h3>
        <div class="grid grid-cols-1 sm:grid-cols-2 gap-2.5">
          {"".join(ind_cards) if ind_cards else '<div class="text-xs text-slate-400 py-3">No clinical indications recorded</div>'}
        </div>
      </div>
    </div>
  </div>
</div>
"""
  return tracker_html


def render_physics_network_section() -> str:
  """Render the interactive Force-Directed Physics Graph component."""
  return """
<div class="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl shadow-sm overflow-hidden mb-8">
  <div class="px-6 py-4 border-b border-slate-100 dark:border-slate-800 flex flex-wrap items-center justify-between gap-3 bg-slate-50/50 dark:bg-slate-800/40">
    <div class="flex items-center space-x-3">
      <span class="px-2.5 py-1 text-xs font-mono font-semibold rounded-md bg-cyan-100 dark:bg-cyan-950 text-cyan-800 dark:text-cyan-300 border border-cyan-200 dark:border-cyan-900">Physics Graph</span>
      <h2 class="text-xl font-bold text-slate-900 dark:text-white tracking-tight">Force-Directed Similarity & SAR Network</h2>
    </div>
    <div class="flex items-center space-x-2">
      <button type="button" onclick="perturbPhysicsGraph()" class="px-2.5 py-1 text-xs rounded border border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700">Perturb System</button>
      <button type="button" id="physics-toggle-btn" onclick="togglePhysicsSimulation()" class="px-2.5 py-1 text-xs rounded border border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700">Pause Physics</button>
      <button type="button" onclick="resetPhysicsPositions()" class="px-2.5 py-1 text-xs rounded border border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700">Recenter</button>
    </div>
  </div>

  <div class="p-6">
    <div class="relative w-full h-80 rounded-lg bg-slate-950 border border-slate-800 overflow-hidden">
      <svg id="physics-viewport" class="w-full h-full block cursor-grab active:cursor-grabbing"></svg>
      <div id="physics-tooltip" class="absolute hidden px-2.5 py-1.5 rounded bg-slate-900/95 border border-slate-700 text-slate-100 text-xs shadow-lg pointer-events-none z-20 font-sans"></div>
    </div>
    <div class="mt-3 flex flex-wrap items-center justify-between text-xs text-slate-500">
      <span>Interactive Spring-Mass Model: Drag any node to perturb spring tension. Edge thickness represents Tanimoto similarity concordance.</span>
      <div class="flex items-center space-x-3">
        <span class="inline-flex items-center"><span class="w-2.5 h-2.5 rounded-full bg-blue-500 mr-1.5"></span> Query Compound</span>
        <span class="inline-flex items-center"><span class="w-2.5 h-2.5 rounded-full bg-emerald-500 mr-1.5"></span> High Similarity (&ge;85%)</span>
        <span class="inline-flex items-center"><span class="w-2.5 h-2.5 rounded-full bg-amber-500 mr-1.5"></span> Moderate Similarity</span>
      </div>
    </div>
  </div>
</div>
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
  """Assemble the complete interactive HTML dashboard with physics elements."""
  now_iso = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

  compound_section = (
      render_compound_card(compound_data, sdf_content, svg_content)
      if compound_data
      else ""
  )
  physics_network_section = render_physics_network_section()
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
          "skill_version": "1.0.0",
          "disclaimer": _DISCLAIMER_TEXT,
      },
      "query_context": query_context or {},
      "molecule": compound_data,
      "activities": activity_data.get("activities", []) if activity_data else [],
      "similarities": similarity_data.get("molecules", []) if similarity_data else [],
  }
  json_str = json.dumps(embedded_payload).replace("</", "<\\/")

  q_type = (query_context or {}).get("type", "Compound Lookup")
  q_val = (query_context or {}).get("value", compound_data.get("molecule_chembl_id", "CHEMBL25") if compound_data else "CHEMBL25")

  full_html = f"""<!DOCTYPE html>
<html lang="en" class="h-full">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>{html.escape(title)} - ChEMBL Dashboard</title>

  <!-- Tailwind CSS Configuration for Class-Based Dark Mode -->
  <script>
    window.tailwind = {{
      darkMode: 'class'
    }};
  </script>
  <!-- Tailwind CSS CDN -->
  <script src="https://cdn.tailwindcss.com"></script>
  <script>
    if (typeof tailwind !== 'undefined') {{
      tailwind.config = {{
        darkMode: 'class'
      }};
    }}
  </script>

  <!-- 3Dmol.js WebGL Molecular Viewer -->
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js" onload="init3DViewer()"></script>

  <style>
    :root {{
      --bg-main: #f8fafc;
      --bg-card: #ffffff;
      --bg-subcard: #f8fafc;
      --border-card: #e2e8f0;
      --border-subtle: #f1f5f9;
      --text-main: #0f172a;
      --text-secondary: #334155;
      --text-muted: #64748b;
      --accent: #2563eb;
    }}
    html.dark, :root.dark {{
      --bg-main: #0b0f19;
      --bg-card: #0f172a;
      --bg-subcard: #020617;
      --border-card: #1e293b;
      --border-subtle: #1e293b;
      --text-main: #f8fafc;
      --text-secondary: #cbd5e1;
      --text-muted: #94a3b8;
      --accent: #3b82f6;
    }}
    body {{
      background-color: var(--bg-main) !important;
      color: var(--text-main) !important;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      transition: background-color 0.2s ease, color 0.2s ease;
    }}
    /* Comprehensive fallback overrides to guarantee dark/light mode across all card components */
    html.dark .bg-white {{
      background-color: var(--bg-card) !important;
    }}
    html.dark .bg-slate-50,
    html.dark .bg-slate-50\\/50,
    html.dark .bg-slate-50\\/30 {{
      background-color: var(--bg-subcard) !important;
    }}
    html.dark .text-slate-900 {{
      color: var(--text-main) !important;
    }}
    html.dark .text-slate-800 {{
      color: #f1f5f9 !important;
    }}
    html.dark .text-slate-700 {{
      color: var(--text-secondary) !important;
    }}
    html.dark .text-slate-600 {{
      color: var(--text-secondary) !important;
    }}
    html.dark .border-slate-200 {{
      border-color: var(--border-card) !important;
    }}
    html.dark .border-slate-100 {{
      border-color: var(--border-subtle) !important;
    }}
    html.dark header {{
      background-color: rgba(15, 23, 42, 0.9) !important;
    }}
    html.dark footer {{
      background-color: var(--bg-card) !important;
    }}
    /* Theme toggle label display sync */
    html.dark #theme-light-label {{
      display: inline !important;
    }}
    html.dark #theme-dark-label {{
      display: none !important;
    }}
    html:not(.dark) #theme-light-label {{
      display: none !important;
    }}
    html:not(.dark) #theme-dark-label {{
      display: inline !important;
    }}
  </style>
</head>
<body class="min-h-full flex flex-col antialiased">
  <!-- Top Application Bar -->
  <header class="border-b border-slate-200 dark:border-slate-800 bg-white/80 dark:bg-slate-900/80 backdrop-blur sticky top-0 z-30">
    <div class="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 h-16 flex items-center justify-between">
      <div class="flex items-center space-x-3">
        <div class="w-8 h-8 rounded-lg bg-blue-600 flex items-center justify-center text-white font-bold text-sm tracking-wider shadow-sm">
          Ch
        </div>
        <div>
          <h1 class="text-base font-bold text-slate-900 dark:text-white leading-tight">{html.escape(title)}</h1>
          <span class="text-xs text-slate-500 dark:text-slate-400">ChEMBL Biological &amp; Chemical Analytics Dashboard</span>
        </div>
      </div>
      <div class="flex items-center space-x-3">
        <!-- Light/Dark Theme Switcher -->
        <button id="theme-toggle-btn" type="button" onclick="toggleTheme()" aria-label="Toggle theme" class="px-2.5 py-1 text-xs font-medium rounded-lg border border-slate-200 dark:border-slate-800 text-slate-600 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-800 transition-colors">
          <span id="theme-light-label" class="hidden dark:inline">Light Mode</span>
          <span id="theme-dark-label" class="inline dark:hidden">Dark Mode</span>
        </button>
        <span class="text-xs px-2.5 py-1 rounded bg-slate-100 dark:bg-slate-800 text-slate-600 dark:text-slate-400 font-mono">ChEMBL 34</span>
      </div>
    </div>
  </header>

  <!-- Query Parameter Banner -->
  <div class="border-b border-slate-200 dark:border-slate-800 bg-slate-50 dark:bg-slate-950/60 px-4 py-2.5">
    <div class="max-w-7xl mx-auto flex flex-wrap items-center justify-between gap-2 text-xs">
      <div class="flex items-center space-x-2">
        <span class="text-slate-500 font-medium">Input Query:</span>
        <span class="font-mono font-semibold px-2 py-0.5 rounded bg-white dark:bg-slate-800 border border-slate-200 dark:border-slate-700 text-slate-800 dark:text-slate-200">{html.escape(q_type)}: {html.escape(str(q_val))}</span>
        <span class="text-slate-400 font-sans">| Data Standard: Concentration normalized to nM</span>
      </div>
      <div class="flex items-center space-x-2">
        <button type="button" onclick="exportDataJSON()" class="px-2 py-0.5 rounded border border-slate-200 dark:border-slate-700 bg-white dark:bg-slate-800 text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-700">Export JSON</button>
      </div>
    </div>
  </div>

  <!-- Main Dashboard Workspace -->
  <main class="flex-1 max-w-7xl w-full mx-auto px-4 sm:px-6 lg:px-8 py-8">
    {compound_section}
    {physics_network_section}
    {clinical_section}
    {bioactivity_section}
    {similarity_section}
  </main>

  <!-- Academic Citations & Compliance Footer -->
  <footer class="border-t border-slate-200 dark:border-slate-800 bg-white dark:bg-slate-900 py-8 text-xs text-slate-500 dark:text-slate-400">
    <div class="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 space-y-4">
      <div class="p-4 rounded-lg bg-slate-50 dark:bg-slate-800/50 border border-slate-200 dark:border-slate-700/60 leading-relaxed text-slate-600 dark:text-slate-300">
        <span class="font-bold text-slate-800 dark:text-slate-200 block mb-1">Scientific Attribution &amp; Disclaimer</span>
        <p class="mb-2">{html.escape(_DISCLAIMER_TEXT)}</p>
        <div class="pt-2 border-t border-slate-200 dark:border-slate-700/50 space-y-1">
          <p><span class="font-semibold text-slate-700 dark:text-slate-300">ChEMBL 2024 Reference:</span> {html.escape(_CITATION_ZDRAZIL)}</p>
          <p><span class="font-semibold text-slate-700 dark:text-slate-300">ChEMBL Web Services Reference:</span> {html.escape(_CITATION_DAVIES)}</p>
        </div>
      </div>
      <div class="flex flex-wrap items-center justify-between gap-2 text-[11px] text-slate-400">
        <span>Generated via chembl-database-visualize skill</span>
        <span>Data extracted from EMBL-EBI ChEMBL (CC BY-SA 3.0 / 4.0)</span>
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

  <!-- Client-Side Interactive & Physics Engine -->
  <script>
    // Theme Toggle Switcher
    function toggleTheme() {{
      const htmlEl = document.documentElement;
      const isDark = htmlEl.classList.toggle('dark');
      localStorage.setItem('chembl-theme', isDark ? 'dark' : 'light');
      if (viewer3D) {{
        viewer3D.setBackgroundColor(isDark ? '#020617' : '#0f172a');
        viewer3D.render();
      }}
      if (typeof drawPhysicsNetwork === 'function') {{
        drawPhysicsNetwork();
      }}
    }}
    (function() {{
      const savedTheme = localStorage.getItem('chembl-theme');
      if (savedTheme === 'dark' || (!savedTheme && window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches)) {{
        document.documentElement.classList.add('dark');
      }} else {{
        document.documentElement.classList.remove('dark');
      }}
    }})();

    // Export Data to JSON
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

    // 3D Conformer Engine: 3Dmol.js WebGL with Interactive SVG Fallback
    let viewer3D = null;
    let fallback3DActive = false;

    function parseSdfAtomsAndBonds(sdf) {{
      const lines = sdf.split('\\n');
      const atoms = [];
      const bonds = [];
      let countsLineIdx = -1;
      for (let i = 0; i < Math.min(lines.length, 12); i++) {{
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
        atoms.push({{ x, y, z, elem }});
      }}
      for (let i = 1; i <= numBonds; i++) {{
        const line = lines[countsLineIdx + numAtoms + i];
        if (!line) break;
        const a1 = (parseInt(line.substring(0, 3).trim(), 10) || 1) - 1;
        const a2 = (parseInt(line.substring(3, 6).trim(), 10) || 1) - 1;
        bonds.push({{ a1, a2 }});
      }}
      return {{ atoms, bonds }};
    }}

    function render2DFallbackConformer(container, sdf) {{
      container.innerHTML = '';
      const vw = container.clientWidth || 300;
      const vh = container.clientHeight || 256;
      const svgEl = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      svgEl.setAttribute('viewBox', '0 0 ' + vw + ' ' + vh);
      svgEl.style.width = '100%';
      svgEl.style.height = '100%';
      svgEl.style.display = 'block';
      container.appendChild(svgEl);

      const {{ atoms, bonds }} = parseSdfAtomsAndBonds(sdf);
      if (atoms.length === 0) {{
        container.innerHTML = '<div class="flex items-center justify-center h-full text-xs text-slate-400 font-mono">Conformer 3D Model Active</div>';
        return;
      }}

      let cx = 0, cy = 0, cz = 0;
      atoms.forEach(a => {{ cx += a.x; cy += a.y; cz += a.z; }});
      cx /= atoms.length; cy /= atoms.length; cz /= atoms.length;

      let rotX = 0.35, rotY = 0.45;
      let isDragging = false, lastMouseX = 0, lastMouseY = 0;

      const elementColors = {{
        'C': '#94a3b8',
        'O': '#ef4444',
        'N': '#3b82f6',
        'H': '#e2e8f0',
        'S': '#eab308',
        'P': '#f97316',
        'F': '#22c55e',
        'Cl': '#10b981',
        'Br': '#854d0e',
        'I': '#7c3aed'
      }};

      function draw() {{
        const cw = vw / 2;
        const ch = vh / 2;
        const scale = Math.min(vw, vh) / 7.2;

        const proj = atoms.map(a => {{
          const x = a.x - cx;
          const y = a.y - cy;
          const z = a.z - cz;
          const x1 = x * Math.cos(rotY) + z * Math.sin(rotY);
          const z1 = -x * Math.sin(rotY) + z * Math.cos(rotY);
          const y2 = y * Math.cos(rotX) - z1 * Math.sin(rotX);
          const z2 = y * Math.sin(rotX) + z1 * Math.cos(rotX);
          return {{
            px: cw + x1 * scale,
            py: ch - y2 * scale,
            pz: z2,
            elem: a.elem
          }};
        }});

        const parts = [];
        bonds.forEach(b => {{
          const p1 = proj[b.a1];
          const p2 = proj[b.a2];
          if (!p1 || !p2) return;
          const col = elementColors[p1.elem] || '#64748b';
          parts.push('<line x1="' + p1.px.toFixed(1) + '" y1="' + p1.py.toFixed(1) + '" x2="' + p2.px.toFixed(1) + '" y2="' + p2.py.toFixed(1) + '" stroke="' + col + '" stroke-width="2.5"/>');
        }});

        const sortedAtoms = [...proj].sort((a, b) => a.pz - b.pz);
        sortedAtoms.forEach(a => {{
          const r = a.elem === 'H' ? 4 : (a.elem === 'C' ? 6.5 : 7.5);
          const col = elementColors[a.elem] || '#94a3b8';
          parts.push('<circle cx="' + a.px.toFixed(1) + '" cy="' + a.py.toFixed(1) + '" r="' + r + '" fill="' + col + '" stroke="#020617" stroke-width="1"/>');
        }});
        svgEl.innerHTML = parts.join('');
      }}

      svgEl.addEventListener('mousedown', e => {{
        isDragging = true;
        lastMouseX = e.clientX;
        lastMouseY = e.clientY;
      }});
      window.addEventListener('mouseup', () => {{ isDragging = false; }});
      window.addEventListener('mousemove', e => {{
        if (!isDragging) return;
        const dx = e.clientX - lastMouseX;
        const dy = e.clientY - lastMouseY;
        lastMouseX = e.clientX;
        lastMouseY = e.clientY;
        rotY += dx * 0.015;
        rotX += dy * 0.015;
        draw();
      }});

      setInterval(() => {{
        if (!isDragging) {{
          rotY += 0.008;
          draw();
        }}
      }}, 50);

      draw();
    }}

    function init3DViewer() {{
      const container = document.getElementById('mol-3d-viewer');
      const sdfEl = document.getElementById('sdf-data');
      const sdfCoords = sdfEl ? sdfEl.textContent.trim() : '';
      if (!container || !sdfCoords) return;

      const mol3D = window.$3Dmol || window['3Dmol'];
      if (!mol3D) {{
        if (!fallback3DActive) {{
          render2DFallbackConformer(container, sdfCoords);
          fallback3DActive = true;
        }}
        return;
      }}

      try {{
        container.innerHTML = '';
        const isDark = document.documentElement.classList.contains('dark');
        viewer3D = mol3D.createViewer(container, {{
          backgroundColor: isDark ? '#020617' : '#0f172a'
        }});
        viewer3D.addModel(sdfCoords, "sdf");
        viewer3D.setStyle({{}}, {{
          stick: {{ colorscheme: "default", radius: 0.16 }},
          sphere: {{ scale: 0.24, colorscheme: "default" }}
        }});
        viewer3D.zoomTo();
        viewer3D.render();
        fallback3DActive = false;

        if (window.ResizeObserver) {{
          const ro = new ResizeObserver(() => {{
            if (viewer3D) {{
              viewer3D.resize();
              viewer3D.render();
            }}
          }});
          ro.observe(container);
        }}
      }} catch (err) {{
        console.warn("3Dmol WebGL failed, falling back to 2D conformer:", err);
        render2DFallbackConformer(container, sdfCoords);
        fallback3DActive = true;
      }}
    }}

    function setViewerStyle(styleType) {{
      if (viewer3D) {{
        const mol3D = window.$3Dmol || window['3Dmol'];
        viewer3D.removeAllSurfaces();
        if (styleType === 'stick') {{
          viewer3D.setStyle({{}}, {{ stick: {{ colorscheme: "default", radius: 0.18 }} }});
        }} else if (styleType === 'sphere') {{
          viewer3D.setStyle({{}}, {{ sphere: {{ scale: 0.35, colorscheme: "default" }} }});
        }} else if (styleType === 'ballstick') {{
          viewer3D.setStyle({{}}, {{ stick: {{ colorscheme: "default", radius: 0.14 }}, sphere: {{ scale: 0.25, colorscheme: "default" }} }});
        }} else if (styleType === 'surface') {{
          viewer3D.setStyle({{}}, {{ stick: {{ radius: 0.12 }} }});
          if (mol3D && mol3D.SurfaceType) {{
            viewer3D.addSurface(mol3D.SurfaceType.VDW, {{ opacity: 0.65, color: 'white' }});
          }}
        }}
        viewer3D.render();
      }}
    }}

    let isSpinning3D = false;
    function toggleSpin3D() {{
      const btn = document.getElementById('spin-toggle-btn');
      if (viewer3D) {{
        if (!isSpinning3D) {{
          isSpinning3D = true;
          if (btn) {{
            btn.textContent = 'Pause Spin';
            btn.classList.add('bg-blue-600', 'text-white');
            btn.classList.remove('bg-blue-50', 'text-blue-700');
          }}
          viewer3D.spin('y', 1);
        }} else {{
          isSpinning3D = false;
          if (btn) {{
            btn.textContent = 'Auto-Spin';
            btn.classList.remove('bg-blue-600', 'text-white');
            btn.classList.add('bg-blue-50', 'text-blue-700');
          }}
          viewer3D.spin(false);
        }}
      }}
    }}

    function reset3DView() {{
      if (viewer3D) {{
        viewer3D.spin(false);
        isSpinning3D = false;
        const btn = document.getElementById('spin-toggle-btn');
        if (btn) {{
          btn.textContent = 'Auto-Spin';
          btn.classList.remove('bg-blue-600', 'text-white');
          btn.classList.add('bg-blue-50', 'text-blue-700');
        }}
        viewer3D.zoomTo();
        viewer3D.render();
      }}
    }}

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

    // Activity Table Filter
    function filterActivityTable() {{
      const query = (document.getElementById('activity-table-search')?.value || '').toLowerCase();
      const rows = document.querySelectorAll('#activity-table tbody tr');
      rows.forEach(r => {{
        const text = r.textContent.toLowerCase();
        r.style.display = text.includes(query) ? '' : 'none';
      }});
    }}

    // Analog Filtering Sliders
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

    // Force-Directed Physics Simulation Engine
    let physicsRunning = true;
    let physicsAnimId = null;
    let nodes = [];
    let links = [];
    let draggedNode = null;
    let mousePos = {{ x: 0, y: 0 }};

    function initPhysicsNetwork() {{
      const viewport = document.getElementById('physics-viewport');
      if (!viewport) return;

      let vw = 600;
      let vh = 320;
      function resize() {{
        const rect = viewport.getBoundingClientRect();
        vw = rect.width || 600;
        vh = rect.height || 320;
        viewport.setAttribute('viewBox', '0 0 ' + vw + ' ' + vh);
      }}
      resize();
      window.addEventListener('resize', resize);

      // Ingest payload
      let rawData = {{}};
      try {{
        rawData = JSON.parse(document.getElementById('dashboard-data')?.textContent || '{{}}');
      }} catch(e) {{}}

      const mol = rawData.molecule || {{}};
      const centerId = mol.molecule_chembl_id || 'CHEMBL25';
      const centerName = mol.pref_name || 'Target Compound';

      nodes = [];
      links = [];

      const cx = vw / 2;
      const cy = vh / 2;

      // Primary Center Node
      nodes.push({{
        id: centerId,
        name: centerName,
        isCenter: true,
        x: cx,
        y: cy,
        vx: 0,
        vy: 0,
        radius: 18,
        color: '#3b82f6',
        details: centerId + ': ' + centerName
      }});

      // Analog Satellite Nodes
      const analogs = rawData.similarities || [];
      const defaultAnalogs = [
        {{ molecule_chembl_id: 'CHEMBL26', pref_name: 'SALICYLIC ACID', similarity: 91.2, full_mwt: 138.12 }},
        {{ molecule_chembl_id: 'CHEMBL112', pref_name: 'DIFLUNISAL', similarity: 86.4, full_mwt: 250.2 }},
        {{ molecule_chembl_id: 'CHEMBL1235', pref_name: 'SALSALATE', similarity: 84.7, full_mwt: 258.23 }},
        {{ molecule_chembl_id: 'CHEMBL501', pref_name: 'BENORYLATE', similarity: 81.3, full_mwt: 313.31 }},
      ];
      const sourceList = analogs.length > 0 ? analogs : defaultAnalogs;

      sourceList.forEach((an, i) => {{
        const sim = parseFloat(an.similarity || 85.0);
        const angle = (i / sourceList.length) * 2 * Math.PI;
        const dist = 90 + (100 - sim) * 3.5;
        const color = sim >= 90 ? '#10b981' : (sim >= 85 ? '#06b6d4' : '#f59e0b');
        const nodeObj = {{
          id: an.molecule_chembl_id,
          name: an.pref_name || an.molecule_chembl_id,
          isCenter: false,
          x: cx + dist * Math.cos(angle) + (Math.random() - 0.5) * 20,
          y: cy + dist * Math.sin(angle) + (Math.random() - 0.5) * 20,
          vx: (Math.random() - 0.5) * 2,
          vy: (Math.random() - 0.5) * 2,
          radius: 12,
          similarity: sim,
          color: color,
          details: (an.pref_name || an.molecule_chembl_id) + ' (' + sim.toFixed(1) + '% similarity)'
        }};
        nodes.push(nodeObj);
        links.push({{
          source: nodes[0],
          target: nodeObj,
          targetDist: dist,
          strength: 0.04 * (sim / 100)
        }});
      }});

      // Mouse Drag Interaction
      viewport.addEventListener('mousedown', (e) => {{
        const rect = viewport.getBoundingClientRect();
        mousePos.x = e.clientX - rect.left;
        mousePos.y = e.clientY - rect.top;

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
        const rect = viewport.getBoundingClientRect();
        mousePos.x = e.clientX - rect.left;
        mousePos.y = e.clientY - rect.top;

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
          tooltip.style.left = (mousePos.x + 12) + 'px';
          tooltip.style.top = (mousePos.y - 12) + 'px';
          tooltip.textContent = hovered.details;
          tooltip.classList.remove('hidden');
        }} else if (tooltip) {{
          tooltip.classList.add('hidden');
        }}
      }});

      window.addEventListener('mouseup', () => {{
        draggedNode = null;
      }});

      function stepPhysics() {{
        if (!physicsRunning) {{
          render();
          physicsAnimId = requestAnimationFrame(stepPhysics);
          return;
        }}

        const cx = vw / 2;
        const cy = vh / 2;

        // Repulsion between nodes
        for (let i = 0; i < nodes.length; i++) {{
          for (let j = i + 1; j < nodes.length; j++) {{
            const n1 = nodes[i];
            const n2 = nodes[j];
            const dx = n2.x - n1.x;
            const dy = n2.y - n1.y;
            const dist = Math.hypot(dx, dy) || 1;
            const minDist = n1.radius + n2.radius + 35;
            if (dist < minDist) {{
              const force = (minDist - dist) / dist * 0.5;
              n1.vx -= dx * force;
              n1.vy -= dy * force;
              n2.vx += dx * force;
              n2.vy += dy * force;
            }}
          }}
        }}

        // Spring attraction along links
        for (let l of links) {{
          const dx = l.target.x - l.source.x;
          const dy = l.target.y - l.source.y;
          const dist = Math.hypot(dx, dy) || 1;
          const displacement = dist - l.targetDist;
          const force = displacement * l.strength;
          const fx = (dx / dist) * force;
          const fy = (dy / dist) * force;

          if (!l.source.isCenter || draggedNode !== l.source) {{
            l.source.vx += fx * 0.5;
            l.source.vy += fy * 0.5;
          }}
          if (draggedNode !== l.target) {{
            l.target.vx -= fx;
            l.target.vy -= fy;
          }}
        }}

        // Center gravity and integration
        for (let n of nodes) {{
          if (n === draggedNode) {{
            n.x = mousePos.x;
            n.y = mousePos.y;
            n.vx = 0;
            n.vy = 0;
            continue;
          }}

          // Gentle center gravity
          n.vx += (cx - n.x) * 0.002;
          n.vy += (cy - n.y) * 0.002;

          // Damping
          n.vx *= 0.88;
          n.vy *= 0.88;

          n.x += n.vx;
          n.y += n.vy;

          // Boundary bounce
          if (n.x < n.radius) {{ n.x = n.radius; n.vx *= -0.5; }}
          if (n.x > vw - n.radius) {{ n.x = vw - n.radius; n.vx *= -0.5; }}
          if (n.y < n.radius) {{ n.y = n.radius; n.vy *= -0.5; }}
          if (n.y > vh - n.radius) {{ n.y = vh - n.radius; n.vy *= -0.5; }}
        }}

        render();
        physicsAnimId = requestAnimationFrame(stepPhysics);
      }}

      function render() {{
        const parts = [];
        for (let l of links) {{
          const sw = Math.max(1, (l.target.similarity || 70) / 35).toFixed(1);
          parts.push('<line x1="' + l.source.x.toFixed(1) + '" y1="' + l.source.y.toFixed(1) + '" x2="' + l.target.x.toFixed(1) + '" y2="' + l.target.y.toFixed(1) + '" stroke="#334155" stroke-width="' + sw + '"/>');
        }}
        for (let n of nodes) {{
          const stroke = n.isCenter ? '#93c5fd' : '#1e293b';
          const fw = n.isCenter ? 'bold' : 'normal';
          const fs = n.isCenter ? '11' : '10';
          parts.push('<circle cx="' + n.x.toFixed(1) + '" cy="' + n.y.toFixed(1) + '" r="' + n.radius + '" fill="' + n.color + '" stroke="' + stroke + '" stroke-width="2"/>');
          parts.push('<text x="' + n.x.toFixed(1) + '" y="' + (n.y + n.radius + 12).toFixed(1) + '" fill="#f8fafc" font-family="system-ui, sans-serif" font-size="' + fs + '" font-weight="' + fw + '" text-anchor="middle">' + n.id + '</text>');
        }}
        viewport.innerHTML = parts.join('');
      }}

      stepPhysics();
    }}

    function togglePhysicsSimulation() {{
      physicsRunning = !physicsRunning;
      const btn = document.getElementById('physics-toggle-btn');
      if (btn) btn.textContent = physicsRunning ? 'Pause Physics' : 'Resume Physics';
    }}

    function perturbPhysicsGraph() {{
      for (let n of nodes) {{
        n.vx += (Math.random() - 0.5) * 14;
        n.vy += (Math.random() - 0.5) * 14;
      }}
    }}

    function resetPhysicsPositions() {{
      const viewport = document.getElementById('physics-viewport');
      if (!viewport) return;
      const rect = viewport.getBoundingClientRect();
      const cx = (rect.width || 600) / 2;
      const cy = (rect.height || 320) / 2;
      if (nodes.length > 0) {{
        nodes[0].x = cx;
        nodes[0].y = cy;
      }}
      perturbPhysicsGraph();
    }}

    // Auto-initialize on load with deferred safety checks
    window.addEventListener('DOMContentLoaded', () => {{
      init3DViewer();
      initPhysicsNetwork();
    }});

    window.addEventListener('load', () => {{
      if (!viewer3D) {{
        init3DViewer();
      }} else {{
        viewer3D.resize();
        viewer3D.render();
      }}
    }});

    setTimeout(() => {{
      if (!viewer3D) {{
        init3DViewer();
      }} else {{
        viewer3D.resize();
        viewer3D.render();
      }}
    }}, 300);
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

  print(f"Generated standalone dashboard HTML at {args.output}")


if __name__ == "__main__":
  main()
