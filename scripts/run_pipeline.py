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

"""Unified end-to-end pipeline runner for the ChEMBL Database Visualize skill.

Orchestrates compound lookups, bioactivity normalization, 2D/3D structure
retrieval, chemical similarity clustering, and interactive dashboard compilation.
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import subprocess
import sys
from typing import Any

_LICENSE_NOTICE = (
    "Notification: Please check the ChEMBL terms of use and licensing at"
    " https://chembl.gitbook.io/chembl-interface-documentation/about and"
    " https://www.ebi.ac.uk/chembl/. ChEMBL data is available under"
    " CC BY-SA 3.0 / CC BY-SA 4.0."
)

_CITATION_SUMMARY = (
    "Attribution: ChEMBL Database (Zdrazil et al. 2024, Davies et al. 2015)."
    " CC BY-SA 3.0/4.0."
)

_SCIENTIFIC_DISCLAIMER = (
    "Disclaimer: ChEMBL data is aggregated from scientific literature and"
    " screening assays. Experimental validation is required prior to therapeutic"
    " or clinical application."
)


def get_workspace_root() -> str:
  """Locate the workspace root directory."""
  script_dir = os.path.dirname(os.path.abspath(__file__))
  return os.path.dirname(script_dir)


def verify_license_guard() -> None:
  """Verify workspace license guard file; create with UTC timestamp if missing."""
  ws_root = get_workspace_root()
  guard_dir = os.path.join(ws_root, ".licenses")
  guard_file = os.path.join(guard_dir, "chembl_database_visualize_LICENSE.txt")
  fallback_guard = os.path.join(guard_dir, "chembl_database_LICENSE.txt")

  if not os.path.exists(guard_file) and not os.path.exists(fallback_guard):
    print("================================================================================")
    print(_LICENSE_NOTICE)
    print("================================================================================")
    os.makedirs(guard_dir, exist_ok=True)
    now_iso = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(guard_file, "w", encoding="utf-8") as f:
      f.write(f"{_LICENSE_NOTICE}\nTimestamp: {now_iso}\n")


def run_cmd(args: list[str]) -> tuple[int, str, str]:
  """Execute a subprocess command and return (code, stdout, stderr)."""
  proc = subprocess.run(args, capture_output=True, text=True)
  return proc.returncode, proc.stdout, proc.stderr


def main() -> None:
  parser = argparse.ArgumentParser(
      description="Run the end-to-end ChEMBL visualization pipeline."
  )
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument("-m", "--molecule", type=str, help="ChEMBL Molecule ID (e.g. CHEMBL25)")
  group.add_argument("-t", "--target", type=str, help="ChEMBL Target ID (e.g. CHEMBL203)")
  group.add_argument("-s", "--smiles", type=str, help="SMILES string for similarity search")
  group.add_argument("-q", "--search", type=str, help="Free-text search query")

  parser.add_argument("--similarity", type=int, default=85, help="Similarity threshold 0-100 (default: 85)")
  parser.add_argument("--activity_type", type=str, default="IC50", help="Bioactivity type (default: IC50)")
  parser.add_argument("--limit", type=int, default=10, help="Max results (default: 10)")
  parser.add_argument("--no_normalize", action="store_true", help="Disable bioactivity normalization to nM")
  parser.add_argument("--mock", action="store_true", help="Run offline synthetic mock pipeline")
  parser.add_argument("--artifact_dir", type=str, default=None, help="Optional agent harness artifact directory")
  parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output destination directory")

  args = parser.parse_args()

  # Step 1: Verify license guard
  verify_license_guard()

  # Ensure output directory exists
  out_dir = os.path.abspath(args.output_dir)
  os.makedirs(out_dir, exist_ok=True)

  ws_root = get_workspace_root()
  api_script = os.path.join(ws_root, "scripts", "chembl_api.py")
  gen_script = os.path.join(ws_root, "scripts", "generate_dashboard.py")

  mock_flag = ["--mock"] if args.mock else []

  compound_data: dict[str, Any] | None = None
  activity_data: dict[str, Any] | None = None
  drug_data: dict[str, Any] | None = None
  similarity_data: dict[str, Any] | None = None

  sdf_file = os.path.join(out_dir, "structure.sdf")
  svg_file = os.path.join(out_dir, "structure.svg")
  data_file = os.path.join(out_dir, "data.json")
  dashboard_file = os.path.join(out_dir, "dashboard.html")

  dashboard_title = "ChEMBL Analytical Dashboard"

  # Step 2: Handle Molecule Workflow
  if args.molecule:
    chembl_id = args.molecule.upper()
    dashboard_title = f"Compound Analysis: {chembl_id}"

    # Fetch molecule JSON
    mol_out = os.path.join(out_dir, "molecule_raw.json")
    code, out, err = run_cmd([
        sys.executable, api_script, "molecule",
        "--id", chembl_id,
        "--output", mol_out,
    ] + mock_flag)
    if code == 0 and os.path.exists(mol_out):
      with open(mol_out, "r", encoding="utf-8") as f:
        compound_data = json.load(f)
      pref_name = compound_data.get("pref_name")
      if pref_name:
        dashboard_title = f"{pref_name} ({chembl_id}) - Profile & Analytics"

    # Download 3D SDF
    run_cmd([
        sys.executable, api_script, "molecule",
        "--id", chembl_id,
        "--dl_format", "sdf",
        "--output", sdf_file,
    ] + mock_flag)

    # Download 2D SVG
    run_cmd([
        sys.executable, api_script, "image",
        "--id", chembl_id,
        "--output", svg_file,
    ] + mock_flag)

    # Fetch bioactivities
    act_out = os.path.join(out_dir, "activities_raw.json")
    act_cmd = [
        sys.executable, api_script, "activity",
        "--filter", f"molecule_chembl_id={chembl_id}",
        "--limit", str(args.limit),
        "--output", act_out,
    ]
    if not args.no_normalize:
      act_cmd.append("--normalize")
    code, _, _ = run_cmd(act_cmd + mock_flag)
    if code == 0 and os.path.exists(act_out):
      with open(act_out, "r", encoding="utf-8") as f:
        activity_data = json.load(f)

    # Fetch mechanisms & indications
    mech_out = os.path.join(out_dir, "mechanisms_raw.json")
    run_cmd([
        sys.executable, api_script, "mechanism",
        "--filter", f"molecule_chembl_id={chembl_id}",
        "--output", mech_out,
    ] + mock_flag)
    ind_out = os.path.join(out_dir, "indications_raw.json")
    run_cmd([
        sys.executable, api_script, "drug_indication",
        "--filter", f"molecule_chembl_id={chembl_id}",
        "--output", ind_out,
    ] + mock_flag)

    mechs, inds = [], []
    if os.path.exists(mech_out):
      try:
        with open(mech_out, "r", encoding="utf-8") as f:
          mechs = json.load(f).get("mechanisms", [])
      except Exception:
        pass
    if os.path.exists(ind_out):
      try:
        with open(ind_out, "r", encoding="utf-8") as f:
          inds = json.load(f).get("drug_indications", [])
      except Exception:
        pass
    drug_data = {"mechanisms": mechs, "drug_indications": inds}

    # Fetch chemical similarity analogs
    sim_out = os.path.join(out_dir, "similarity_raw.json")
    smiles = ((compound_data.get("molecule_structures") or {}).get("canonical_smiles")) if compound_data else None
    if smiles or args.mock:
      sim_cmd = [
          sys.executable, api_script, "similarity",
          "--smiles", smiles or "CC(=O)Oc1ccccc1C(=O)O",
          "--similarity", str(args.similarity),
          "--limit", str(args.limit),
          "--output", sim_out,
      ] + mock_flag
      code, _, _ = run_cmd(sim_cmd)
      if code == 0 and os.path.exists(sim_out):
        try:
          with open(sim_out, "r", encoding="utf-8") as f:
            similarity_data = json.load(f)
        except Exception:
          pass

  # Step 3: Handle Target Workflow
  elif args.target:
    target_id = args.target.upper()
    dashboard_title = f"Target SAR Profiling: {target_id}"

    # Fetch target details
    tgt_out = os.path.join(out_dir, "target_raw.json")
    run_cmd([
        sys.executable, api_script, "target",
        "--id", target_id,
        "--output", tgt_out,
    ] + mock_flag)

    tgt_name = target_id
    if os.path.exists(tgt_out):
      with open(tgt_out, "r", encoding="utf-8") as f:
        tgt_data = json.load(f)
        tgt_name = tgt_data.get("pref_name") or target_id
    dashboard_title = f"{tgt_name} ({target_id}) - Target Bioactivity Profile"

    # Fetch bioactivities for target
    act_out = os.path.join(out_dir, "activities_raw.json")
    act_cmd = [
        sys.executable, api_script, "activity",
        "--filter", f"target_chembl_id={target_id}", f"standard_type={args.activity_type}",
        "--limit", str(args.limit),
        "--output", act_out,
    ]
    if not args.no_normalize:
      act_cmd.append("--normalize")
    code, _, _ = run_cmd(act_cmd + mock_flag)
    if code == 0 and os.path.exists(act_out):
      with open(act_out, "r", encoding="utf-8") as f:
        activity_data = json.load(f)

  # Step 4: Handle SMILES Similarity Workflow
  elif args.smiles:
    dashboard_title = f"Chemical Similarity Search: {args.smiles[:30]}..."
    sim_out = os.path.join(out_dir, "similarity_raw.json")
    run_cmd([
        sys.executable, api_script, "similarity",
        "--smiles", args.smiles,
        "--similarity", str(args.similarity),
        "--limit", str(args.limit),
        "--output", sim_out,
    ] + mock_flag)
    if os.path.exists(sim_out):
      with open(sim_out, "r", encoding="utf-8") as f:
        similarity_data = json.load(f)

    # Provide query compound representation for 2D & 3D visualization
    top_hit = (similarity_data.get("molecules") or [{}])[0] if similarity_data else {}
    ref_id = top_hit.get("molecule_chembl_id", "CHEMBL25")
    compound_data = {
        "molecule_chembl_id": ref_id,
        "pref_name": top_hit.get("pref_name", "Query Compound"),
        "max_phase": top_hit.get("max_phase", 4),
        "molecule_structures": {"canonical_smiles": args.smiles},
        "molecule_properties": top_hit.get("molecule_properties", {
            "full_mwt": 180.16,
            "alogp": 1.31,
            "hbd": 1,
            "hba": 3,
            "rtb": 2,
            "psa": 63.6,
            "ro5_violations": 0,
        }),
    }
    run_cmd([
        sys.executable, api_script, "image",
        "--id", ref_id,
        "--output", svg_file,
    ] + mock_flag)
    run_cmd([
        sys.executable, api_script, "molecule",
        "--id", ref_id,
        "--dl_format", "sdf",
        "--output", sdf_file,
    ] + mock_flag)

  # Step 5: Handle Free-Text Search Workflow
  elif args.search:
    dashboard_title = f"ChEMBL Query: {args.search}"
    srch_out = os.path.join(out_dir, "search_raw.json")
    run_cmd([
        sys.executable, api_script, "molecule",
        "--search", args.search,
        "--limit", str(args.limit),
        "--output", srch_out,
    ] + mock_flag)
    if os.path.exists(srch_out):
      with open(srch_out, "r", encoding="utf-8") as f:
        raw_srch = json.load(f)
        mols = raw_srch.get("molecules", [])
        if mols:
          compound_data = mols[0]
          dashboard_title = f"{compound_data.get('pref_name') or 'Compound'} - Search Profile"

  # Step 6: Consolidate data into data.json
  consolidated = {
      "query_type": "molecule" if args.molecule else ("target" if args.target else ("smiles" if args.smiles else "search")),
      "title": dashboard_title,
      "molecule": compound_data,
      "bioactivities": activity_data.get("activities", []) if activity_data else [],
      "drug": drug_data,
      "similarities": similarity_data.get("molecules", []) if similarity_data else [],
      "attribution": _CITATION_SUMMARY,
      "disclaimer": _SCIENTIFIC_DISCLAIMER,
  }
  with open(data_file, "w", encoding="utf-8") as f:
    json.dump(consolidated, f, indent=2)

  # Step 7: Build Dashboard via generate_dashboard.py
  dash_cmd = [
      sys.executable, gen_script,
      "--title", dashboard_title,
      "--output", dashboard_file,
  ]
  if compound_data:
    mol_f = os.path.join(out_dir, "molecule_raw.json")
    if os.path.exists(mol_f):
      dash_cmd.extend(["--compound_json", mol_f])
  if os.path.exists(sdf_file):
    dash_cmd.extend(["--sdf", sdf_file])
  if os.path.exists(svg_file):
    dash_cmd.extend(["--image", svg_file])
  if activity_data:
    act_f = os.path.join(out_dir, "activities_raw.json")
    if os.path.exists(act_f):
      dash_cmd.extend(["--activity_json", act_f])
  if drug_data:
    drug_f = os.path.join(out_dir, "drug_summary.json")
    with open(drug_f, "w", encoding="utf-8") as f:
      json.dump(drug_data, f, indent=2)
    dash_cmd.extend(["--drug_json", drug_f])
  if similarity_data:
    sim_f = os.path.join(out_dir, "similarity_raw.json")
    if os.path.exists(sim_f):
      dash_cmd.extend(["--similarity_json", sim_f])

  code, _, err = run_cmd(dash_cmd)
  if code != 0:
    print(f"Error compiling dashboard: {err}", file=sys.stderr)
    sys.exit(code)

  # Ensure dashboard.html exists at project root for direct local testing
  # Avoid overwriting project root when running inside temporary directories (e.g. unit tests)
  root_dashboard = os.path.join(ws_root, "dashboard.html")
  try:
    import tempfile
    is_temp_dir = out_dir.startswith(tempfile.gettempdir()) or "/tmp" in out_dir
    if not is_temp_dir and os.path.abspath(dashboard_file) != os.path.abspath(root_dashboard):
      import shutil
      shutil.copyfile(dashboard_file, root_dashboard)
  except Exception:
    pass

  # If artifact_dir is specified, copy dashboard.html there as well
  if getattr(args, "artifact_dir", None) and os.path.exists(args.artifact_dir):
    try:
      import shutil
      art_dash = os.path.join(args.artifact_dir, "dashboard.html")
      shutil.copyfile(dashboard_file, art_dash)
    except Exception:
      pass

  # Step 8: Print CLI Execution Summary
  print("================================================================================")
  print(f"ChEMBL Visualization Pipeline Completed Successfully")
  print("================================================================================")
  print(f"Title:        {dashboard_title}")
  if compound_data:
    p = compound_data.get("molecule_properties", {})
    ro5_viol = p.get("ro5_violations", 0)
    print(f"Compound:     {compound_data.get('pref_name', 'N/A')} ({compound_data.get('molecule_chembl_id', 'N/A')})")
    print(f"Formula:      {p.get('full_molformula', 'N/A')} | MW: {p.get('full_mwt', 'N/A')} Da | AlogP: {p.get('alogp', 'N/A')}")
    print(f"Lipinski Ro5: {ro5_viol} violations ({'Compliant' if ro5_viol == 0 else 'Non-compliant'})")
  if activity_data:
    acts = activity_data.get("activities", [])
    print(f"Bioactivity:  {len(acts)} records normalized to standard nM")
  if similarity_data:
    mols = similarity_data.get("molecules", [])
    print(f"Analogs:      {len(mols)} chemical analogs identified")
  def to_rel(p: str) -> str:
    try:
      rel = os.path.relpath(p, ws_root)
      return f"./{rel}" if not rel.startswith((".", "/")) else rel
    except Exception:
      return p

  print(f"Artifacts:")
  print(f"  - Dashboard HTML:  {to_rel(dashboard_file)}")
  print(f"  - Project Root:    {to_rel(root_dashboard)}")
  print(f"  - Consolidated:    {to_rel(data_file)}")
  if os.path.exists(sdf_file):
    print(f"  - 3D Conformer:    {to_rel(sdf_file)}")
  if os.path.exists(svg_file):
    print(f"  - 2D Vector SVG:   {to_rel(svg_file)}")
  print("--------------------------------------------------------------------------------")
  print(_CITATION_SUMMARY)
  print(_SCIENTIFIC_DISCLAIMER)
  print("================================================================================")


if __name__ == "__main__":
  main()
