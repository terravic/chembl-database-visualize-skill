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

"""Chemical similarity and analog gallery component.

Transforms ChEMBL similarity and substructure search results into responsive
analog cards with Tanimoto similarity bars, real-time client-side filtering
sliders, and comparative chemical matrices.
"""

from __future__ import annotations

import argparse
import html
import json
import os
import sys
from typing import Any


def render_similarity_gallery(similarity_data: dict[str, Any]) -> str:
  """Render responsive analog gallery and filter controls."""
  molecules = similarity_data.get("molecules", [])
  total_count = len(molecules)

  cards = []
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
    max_phase = mol.get("max_phase")
    phase_str = f"Phase {max_phase}" if max_phase is not None else "Preclinical"

    # Color code similarity bar
    if sim >= 90.0:
      bar_color = "bg-emerald-500"
      badge_color = "bg-emerald-50 text-emerald-700 border-emerald-200 dark:bg-emerald-950 dark:text-emerald-300 dark:border-emerald-900"
    elif sim >= 80.0:
      bar_color = "bg-blue-500"
      badge_color = "bg-blue-50 text-blue-700 border-blue-200 dark:bg-blue-950 dark:text-blue-300 dark:border-blue-900"
    else:
      bar_color = "bg-amber-500"
      badge_color = "bg-amber-50 text-amber-700 border-amber-200 dark:bg-amber-950 dark:text-amber-300 dark:border-amber-900"

    cards.append(
        f"""
        <div class="analog-card border border-slate-200 dark:border-slate-800 rounded-lg p-4 bg-white dark:bg-slate-950 flex flex-col justify-between hover:shadow-md transition-shadow" data-similarity="{sim:.1f}" data-mwt="{mwt:.1f}">
          <div>
            <div class="flex items-center justify-between gap-2 mb-2">
              <span class="px-2 py-0.5 text-xs font-mono font-bold rounded bg-slate-100 dark:bg-slate-800 text-slate-800 dark:text-slate-200">{html.escape(cid)}</span>
              <span class="px-2 py-0.5 text-[11px] font-semibold rounded-full border {badge_color}">{sim_pct} Match</span>
            </div>
            <h4 class="text-sm font-semibold text-slate-900 dark:text-white truncate mb-1" title="{html.escape(name)}">{html.escape(name)}</h4>
            <div class="w-full bg-slate-100 dark:bg-slate-800 rounded-full h-1.5 mb-3 overflow-hidden">
              <div class="{bar_color} h-1.5 rounded-full" style="width: {min(100.0, sim):.1f}%"></div>
            </div>

            <!-- Structure sketch placeholder / vector -->
            <div class="h-28 bg-slate-50 dark:bg-slate-900 rounded border border-slate-100 dark:border-slate-800 flex items-center justify-center p-2 mb-3 text-center">
              <span class="text-xs font-mono text-slate-400 break-all line-clamp-3">{html.escape(smiles) if smiles else "Structure coordinates unavailable"}</span>
            </div>

            <!-- Chemical stats grid -->
            <div class="grid grid-cols-2 gap-2 text-xs mb-3">
              <div class="p-1.5 rounded bg-slate-50 dark:bg-slate-800/50">
                <span class="text-[10px] text-slate-500 block">MW</span>
                <span class="font-mono font-medium text-slate-800 dark:text-slate-200">{mwt:.1f} Da</span>
              </div>
              <div class="p-1.5 rounded bg-slate-50 dark:bg-slate-800/50">
                <span class="text-[10px] text-slate-500 block">AlogP</span>
                <span class="font-mono font-medium text-slate-800 dark:text-slate-200">{alogp:.2f}</span>
              </div>
            </div>
          </div>

          <div class="pt-2 border-t border-slate-100 dark:border-slate-800 flex items-center justify-between text-[11px]">
            <span class="text-slate-500 font-medium">{html.escape(phase_str)}</span>
            <span class="text-slate-400 font-mono text-[10px] truncate max-w-[120px]" title="{html.escape(inchi_key)}">{html.escape(inchi_key)}</span>
          </div>
        </div>
        """
    )

  cards_html = (
      "".join(cards)
      if cards
      else '<div class="col-span-full py-10 text-center text-slate-400 text-sm">No chemical analogs found for the specified query.</div>'
  )

  gallery_html = f"""
<div class="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl shadow-sm overflow-hidden mb-8">
  <div class="px-6 py-4 border-b border-slate-100 dark:border-slate-800 flex flex-wrap items-center justify-between gap-3 bg-slate-50/50 dark:bg-slate-800/40">
    <div class="flex items-center space-x-3">
      <span class="px-2.5 py-1 text-xs font-mono font-semibold rounded-md bg-emerald-100 dark:bg-emerald-950 text-emerald-800 dark:text-emerald-300 border border-emerald-200 dark:border-emerald-900">Analog Search</span>
      <h2 class="text-xl font-bold text-slate-900 dark:text-white tracking-tight">Chemical Similarity & Concordance Grid</h2>
    </div>
    <div class="flex items-center space-x-2 text-xs">
      <span class="text-slate-500">Showing <span id="visible-analog-count" class="font-bold text-slate-800 dark:text-slate-200">{total_count}</span> of {total_count} analogs</span>
    </div>
  </div>

  <!-- Interactive Controls Bar -->
  <div class="p-4 border-b border-slate-100 dark:border-slate-800 bg-slate-50/30 dark:bg-slate-800/20 flex flex-wrap items-center gap-6 text-xs">
    <div class="flex items-center space-x-3">
      <label for="similarity-filter" class="font-medium text-slate-700 dark:text-slate-300">Minimum Similarity:</label>
      <input type="range" id="similarity-filter" min="50" max="100" value="70" step="1" oninput="updateAnalogFilters()" class="w-32 accent-blue-600"/>
      <span id="similarity-filter-val" class="font-mono font-bold text-blue-600 dark:text-blue-400 w-10">70%</span>
    </div>

    <div class="flex items-center space-x-3">
      <label for="mwt-filter" class="font-medium text-slate-700 dark:text-slate-300">Max MW (Da):</label>
      <input type="range" id="mwt-filter" min="100" max="800" value="800" step="25" oninput="updateAnalogFilters()" class="w-32 accent-blue-600"/>
      <span id="mwt-filter-val" class="font-mono font-bold text-blue-600 dark:text-blue-400 w-14">800 Da</span>
    </div>

    <button type="button" onclick="resetAnalogFilters()" class="px-2.5 py-1 text-xs rounded border border-slate-200 dark:border-slate-700 text-slate-600 dark:text-slate-400 hover:bg-slate-100 dark:hover:bg-slate-800">Reset Filters</button>
  </div>

  <!-- Analog Card Grid -->
  <div id="analog-grid" class="p-6 grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-4">
    {cards_html}
  </div>
</div>
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
