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

"""ChEMBL REST API client CLI with offline mock engine.

CLI tool covering all ChEMBL web services API endpoints. Writes JSON or binary
output to a file specified by --output. Enforces 5.0 QPS rate limiting and
retries transient errors (HTTP 429, 503) with exponential backoff. Provides an
offline mock mode (--mock) using synthetic fixtures.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from typing import Any
import urllib.error
import urllib.parse
import urllib.request

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
_LICENSE_NOTICE = (
    "Data from the ChEMBL Database. You MUST notify the user"
    " that this data comes from ChEMBL and advise them to"
    " review the ChEMBL licensing terms."
)

SEARCHABLE_ENDPOINTS = frozenset([
    "activity",
    "assay",
    "chembl_id_lookup",
    "document",
    "molecule",
    "protein_classification",
    "target",
])

ENDPOINT_MAP = {
    "activity": "activity",
    "activity_supp": "activity_supplementary_data_by_activity",
    "assay": "assay",
    "assay_class": "assay_class",
    "atc_class": "atc_class",
    "binding_site": "binding_site",
    "biotherapeutic": "biotherapeutic",
    "cell_line": "cell_line",
    "chembl_id_lookup": "chembl_id_lookup",
    "chembl_release": "chembl_release",
    "compound_record": "compound_record",
    "compound_structural_alert": "compound_structural_alert",
    "document": "document",
    "document_similarity": "document_similarity",
    "drug": "drug",
    "drug_indication": "drug_indication",
    "drug_warning": "drug_warning",
    "go_slim": "go_slim",
    "mechanism": "mechanism",
    "metabolism": "metabolism",
    "molecule": "molecule",
    "molecule_form": "molecule_form",
    "organism": "organism",
    "protein_classification": "protein_classification",
    "source": "source",
    "target": "target",
    "target_component": "target_component",
    "target_relation": "target_relation",
    "tissue": "tissue",
    "xref_source": "xref_source",
}

UNIT_CONVERSION_TO_NM = {
    "nm": 1.0,
    "um": 1e3,
    "µm": 1e3,
    "mm": 1e6,
    "m": 1e9,
    "pm": 1e-3,
    "fm": 1e-6,
}


class StandardHttpClient:
  """Built-in polite HTTP client using standard library urllib."""

  def __init__(self, base_url: str, qps: float = 5.0):
    self.base_url = base_url.rstrip("/")
    self.min_interval = 1.0 / qps if qps > 0 else 0.0
    self.last_request_time = 0.0

  def _rate_limit(self) -> None:
    now = time.time()
    elapsed = now - self.last_request_time
    if elapsed < self.min_interval:
      time.sleep(self.min_interval - elapsed)
    self.last_request_time = time.time()

  def fetch_json(self, url: str) -> dict[str, Any]:
    content = self.fetch_bytes(url)
    return json.loads(content.decode("utf-8"))

  def fetch_bytes(self, url: str) -> bytes:
    max_retries = 3
    base_backoff = 1.0
    for attempt in range(max_retries + 1):
      self._rate_limit()
      req = urllib.request.Request(
          url,
          headers={
              "User-Agent": (
                  "chembl-database-visualize/1.0.0 (Automated Scientific Agent)"
              ),
              "Accept": "*/*",
          },
      )
      try:
        with urllib.request.urlopen(req, timeout=30) as resp:
          return resp.read()
      except urllib.error.HTTPError as e:
        if e.code in (429, 503) and attempt < max_retries:
          sleep_time = base_backoff * (2**attempt) + random.uniform(0.1, 0.5)
          time.sleep(sleep_time)
          continue
        body = e.read() if hasattr(e, "read") else b""
        raise RuntimeError(
            f"HTTP error {e.code}: {e.reason} ({body.decode('utf-8', errors='replace')[:200]})"
        ) from e
      except Exception as e:
        if attempt < max_retries:
          time.sleep(base_backoff * (2**attempt))
          continue
        raise RuntimeError(f"Network error: {e}") from e
    raise RuntimeError(f"Max retries exceeded for {url}")


# Initialize client: use polite_http if available, otherwise StandardHttpClient
try:
  from polite_http import http_client
  _CLIENT = http_client.HttpClient(BASE_URL, qps=5.0)
except Exception:
  _CLIENT = StandardHttpClient(BASE_URL, qps=5.0)


def _get_project_root() -> str:
  """Resolve project root directory."""
  script_dir = os.path.dirname(os.path.abspath(__file__))
  return os.path.dirname(script_dir)


def _load_sample_fixture(filename: str) -> Any:
  """Load a synthetic fixture file from sample_data/ or assets/."""
  root = _get_project_root()
  candidates = [
      os.path.join(root, "sample_data", filename),
      os.path.join(root, "assets", filename),
      os.path.join(root, "skills", "chembl_database_visualize", "sample_data", filename),
  ]
  for path in candidates:
    if os.path.exists(path):
      if filename.endswith(".json"):
        with open(path, "r", encoding="utf-8") as f:
          return json.load(f)
      else:
        with open(path, "rb") as f:
          return f.read()
  return None


def _write_json(data: Any, output_path: str) -> None:
  """Write a Python object as indented JSON to output_path."""
  out_dir = os.path.dirname(output_path)
  if out_dir:
    os.makedirs(out_dir, exist_ok=True)

  if isinstance(data, dict):
    data["_license_notice"] = (
        "Data from the ChEMBL Database. Please review the licensing terms at"
        " https://www.ebi.ac.uk/chembl/"
    )

  with open(output_path, "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2)
    f.write("\n")

  print(
      json.dumps(
          {
              "status": "success",
              "output_file": output_path,
              "size_bytes": os.path.getsize(output_path),
              "license_notice": _LICENSE_NOTICE,
          },
          indent=2,
      )
  )


def _normalize_activity_records(data: dict[str, Any]) -> dict[str, Any]:
  """Normalize activity values to nM."""
  activities = data.get("activities", [])
  for act in activities:
    std_val = act.get("standard_value")
    std_units = (act.get("standard_units") or "").strip().lower()
    if std_val is not None:
      try:
        val_float = float(std_val)
        if std_units in UNIT_CONVERSION_TO_NM:
          factor = UNIT_CONVERSION_TO_NM[std_units]
          norm_nm = val_float * factor
          act["normalized_value_nM"] = round(norm_nm, 4)
        else:
          act["normalized_value_nM"] = None
      except (ValueError, TypeError):
        act["normalized_value_nM"] = None
    else:
      act["normalized_value_nM"] = None
  return data


def _make_request(url: str) -> dict[str, Any]:
  """Send HTTP GET request to url and return parsed JSON."""
  try:
    if hasattr(_CLIENT, "fetch_json"):
      return _CLIENT.fetch_json(url)
    content = _CLIENT.fetch_bytes(url)
    return json.loads(content.decode("utf-8"))
  except Exception as e:
    return {
        "status": "error",
        "message": f"Request failed: {e}",
    }


def _download_binary(url: str, output_path: str) -> dict[str, Any]:
  """Download binary content from url and save to output_path."""
  try:
    content = _CLIENT.fetch_bytes(url)
    out_dir = os.path.dirname(output_path)
    if out_dir:
      os.makedirs(out_dir, exist_ok=True)
    with open(output_path, "wb") as f:
      f.write(content)
    return {
        "status": "success",
        "message": f"Saved to {output_path}",
        "size_bytes": len(content),
        "license_notice": _LICENSE_NOTICE,
    }
  except Exception as e:
    return {"status": "error", "message": f"Download failed: {e}"}


def _build_url(
    endpoint: str,
    item_id: str | None = None,
    item_ids: list[str] | None = None,
    search_query: str | None = None,
    limit: int | None = None,
    offset: int | None = None,
    filters: list[str] | None = None,
) -> str:
  """Construct a ChEMBL API URL based on arguments."""
  api_name = ENDPOINT_MAP.get(endpoint, endpoint)
  params: dict[str, Any] = {}

  if item_id:
    url = f"{BASE_URL}/{api_name}/{item_id}.json"
  elif item_ids:
    ids_str = ";".join(item_ids)
    url = f"{BASE_URL}/{api_name}/set/{ids_str}.json"
  elif search_query:
    url = f"{BASE_URL}/{api_name}/search.json"
    params["q"] = search_query
  else:
    url = f"{BASE_URL}/{api_name}.json"

  if limit is not None:
    params["limit"] = limit
  if offset is not None:
    params["offset"] = offset
  if filters:
    for f in filters:
      if "=" in f:
        k, v = f.split("=", 1)
        params[k] = v

  if params:
    url += "?" + urllib.parse.urlencode(params)
  return url


def cmd_status(args: argparse.Namespace) -> None:
  """Check ChEMBL API status."""
  if getattr(args, "mock", False):
    result = {
        "status": "UP",
        "chembl_release": "ChEMBL_34",
        "mock": True,
        "_license_notice": _LICENSE_NOTICE,
    }
  else:
    url = f"{BASE_URL}/status.json"
    result = _make_request(url)
  _write_json(result, args.output)


def cmd_generic(args: argparse.Namespace) -> None:
  """Handle generic queries for standard ChEMBL endpoints."""
  cmd = args.command
  is_mock = getattr(args, "mock", False)

  if is_mock:
    result = _handle_generic_mock(args)
  else:
    ids_list = args.ids.split(";") if args.ids else None
    url = _build_url(
        endpoint=cmd,
        item_id=args.id,
        item_ids=ids_list,
        search_query=args.search,
        limit=args.limit,
        offset=args.offset,
        filters=args.filter,
    )
    result = _make_request(url)

  if cmd == "activity" and getattr(args, "normalize", False):
    result = _normalize_activity_records(result)

  _write_json(result, args.output)


def _handle_generic_mock(args: argparse.Namespace) -> dict[str, Any]:
  """Synthesize or load mock response for generic endpoints."""
  cmd = args.command
  cid = (args.id or "").upper()

  if cmd == "molecule":
    if cid == "CHEMBL25" or not cid:
      fixture = _load_sample_fixture("chembl25_aspirin.json")
      if fixture:
        return fixture
    return {
        "molecule_chembl_id": cid or "CHEMBL_MOCK_1",
        "pref_name": f"SYNTHETIC COMPOUND {cid}",
        "max_phase": 2,
        "molecule_structures": {
            "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "standard_inchi_key": "MOCK-INCHIKEY-SYNTHETIC",
        },
        "molecule_properties": {
            "full_mwt": 320.4,
            "alogp": 2.5,
            "hbd": 2,
            "hba": 4,
            "psa": 65.0,
            "rtb": 3,
            "ro5_violations": 0,
        },
    }

  if cmd == "activity":
    fixture = _load_sample_fixture("egfr_activities.json")
    if fixture:
      return fixture
    return {
        "activities": [
            {
                "activity_id": 99901,
                "assay_chembl_id": "CHEMBL_ASSAY_MOCK",
                "target_chembl_id": "CHEMBL203",
                "target_pref_name": "Epidermal growth factor receptor erbB1",
                "target_organism": "Homo sapiens",
                "standard_type": "IC50",
                "standard_value": 42.0,
                "standard_units": "nM",
                "standard_relation": "=",
            }
        ],
        "page_meta": {"limit": 10, "offset": 0, "total_count": 1},
    }

  if cmd == "target":
    return {
        "target_chembl_id": cid or "CHEMBL203",
        "pref_name": "Epidermal growth factor receptor erbB1",
        "target_type": "SINGLE PROTEIN",
        "organism": "Homo sapiens",
        "target_components": [
            {
                "accession": "P00533",
                "component_id": 1,
                "component_name": "Epidermal growth factor receptor",
            }
        ],
    }

  if cmd == "mechanism":
    return {
        "mechanisms": [
            {
                "mechanism_of_action": "Cyclooxygenase-1 inhibitor",
                "action_type": "INHIBITOR",
                "target_chembl_id": "CHEMBL221",
                "target_name": "Cyclooxygenase-1",
                "direct_interaction": True,
                "disease_efficacy": True,
            }
        ]
    }

  if cmd == "drug_indication":
    return {
        "drug_indications": [
            {
                "efo_term": "pain",
                "mesh_heading": "Pain",
                "max_phase_for_ind": 4,
            },
            {
                "efo_term": "fever",
                "mesh_heading": "Fever",
                "max_phase_for_ind": 4,
            },
        ]
    }

  # Generic fallback for other endpoints
  plural = f"{cmd}s"
  return {
      plural: [],
      "page_meta": {"limit": args.limit or 5, "offset": args.offset or 0, "total_count": 0},
  }


def cmd_similarity(args: argparse.Namespace) -> None:
  """Run server-side similarity search by SMILES."""
  if getattr(args, "mock", False):
    fixture = _load_sample_fixture("similarity_aspirin.json")
    result = fixture if fixture else {"molecules": [], "page_meta": {"total_count": 0}}
  else:
    smiles_encoded = urllib.parse.quote(args.smiles, safe="")
    url = f"{BASE_URL}/similarity/{smiles_encoded}/{args.similarity}.json"
    params: dict[str, Any] = {}
    if args.limit:
      params["limit"] = args.limit
    if args.offset:
      params["offset"] = args.offset
    if params:
      url += "?" + urllib.parse.urlencode(params)
    result = _make_request(url)
  _write_json(result, args.output)


def cmd_substructure(args: argparse.Namespace) -> None:
  """Run server-side substructure search by SMILES."""
  if getattr(args, "mock", False):
    fixture = _load_sample_fixture("similarity_aspirin.json")
    result = fixture if fixture else {"molecules": [], "page_meta": {"total_count": 0}}
  else:
    smiles_encoded = urllib.parse.quote(args.smiles, safe="")
    url = f"{BASE_URL}/substructure/{smiles_encoded}.json"
    params: dict[str, Any] = {}
    if args.limit:
      params["limit"] = args.limit
    if args.offset:
      params["offset"] = args.offset
    if params:
      url += "?" + urllib.parse.urlencode(params)
    result = _make_request(url)
  _write_json(result, args.output)


def cmd_image(args: argparse.Namespace) -> None:
  """Download compound 2D structure image."""
  out_dir = os.path.dirname(args.output)
  if out_dir:
    os.makedirs(out_dir, exist_ok=True)

  if getattr(args, "mock", False):
    fixture_bytes = _load_sample_fixture("chembl25_aspirin.svg")
    if not fixture_bytes:
      fixture_bytes = _load_sample_fixture("structure_fallback.svg")
    if not fixture_bytes:
      fixture_bytes = b"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><text y='50'>Structure</text></svg>"
    with open(args.output, "wb") as f:
      f.write(fixture_bytes)
    result = {
        "status": "success",
        "message": f"Saved mock image to {args.output}",
        "size_bytes": len(fixture_bytes),
        "license_notice": _LICENSE_NOTICE,
    }
  else:
    url = f"{BASE_URL}/image/{args.id}"
    params: dict[str, Any] = {}
    if args.dimensions:
      params["dimensions"] = args.dimensions
    if args.engine:
      params["engine"] = args.engine
    if args.img_format:
      params["format"] = args.img_format
    if params:
      url += "?" + urllib.parse.urlencode(params)
    result = _download_binary(url, args.output)

  print(json.dumps(result, indent=2))


def cmd_molecule_download(args: argparse.Namespace) -> None:
  """Download molecule structure file in SDF or MOL format."""
  fmt = args.dl_format
  chembl_id = args.id
  if not chembl_id:
    error = {
        "status": "error",
        "message": "--id is required for --dl_format",
    }
    print(json.dumps(error, indent=2))
    sys.exit(1)

  output = args.output if args.output else f"{chembl_id}.{fmt}"
  out_dir = os.path.dirname(output)
  if out_dir:
    os.makedirs(out_dir, exist_ok=True)

  if getattr(args, "mock", False):
    fixture_bytes = _load_sample_fixture("chembl25_aspirin.sdf")
    if not fixture_bytes:
      fixture_bytes = f"{chembl_id}\nMock 3D\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n$$$$\n".encode("utf-8")
    with open(output, "wb") as f:
      f.write(fixture_bytes)
    result = {
        "status": "success",
        "message": f"Saved mock {fmt.upper()} to {output}",
        "size_bytes": len(fixture_bytes),
        "license_notice": _LICENSE_NOTICE,
    }
  else:
    url = f"{BASE_URL}/molecule/{chembl_id}.{fmt}"
    result = _download_binary(url, output)

  print(json.dumps(result, indent=2))


def _add_common_args(parser: argparse.ArgumentParser) -> None:
  """Add shared arguments to a subparser."""
  parser.add_argument("--id", type=str, help="Single ChEMBL ID or numeric ID")
  parser.add_argument(
      "--ids",
      type=str,
      help="Semicolon-separated list of IDs for batch fetch",
  )
  parser.add_argument(
      "--search",
      type=str,
      help="Free-text search query (only for searchable endpoints)",
  )
  parser.add_argument(
      "--limit",
      type=int,
      default=5,
      help="Max results to return (default: 5)",
  )
  parser.add_argument(
      "--offset",
      type=int,
      default=None,
      help="Pagination offset",
  )
  parser.add_argument(
      "--filter",
      type=str,
      nargs="*",
      help="Filter as KEY=VALUE pairs",
  )
  parser.add_argument(
      "--output",
      type=str,
      required=True,
      help="Output JSON file path (required)",
  )
  parser.add_argument(
      "--mock",
      action="store_true",
      help="Enable offline synthetic mock execution",
  )


def build_parser() -> argparse.ArgumentParser:
  """Build the top-level argparse parser with all subcommands."""
  parser = argparse.ArgumentParser(
      description=(
          "ChEMBL REST API client with offline mock mode. Query bioactive"
          " molecules, targets, activities, and more. All output is written to"
          " --output file."
      )
  )
  subparsers = parser.add_subparsers(
      dest="command", help="API endpoint to query"
  )

  for cmd_name in sorted(ENDPOINT_MAP.keys()):
    api_name = ENDPOINT_MAP[cmd_name]
    searchable = " (searchable)" if api_name in SEARCHABLE_ENDPOINTS else ""
    sp = subparsers.add_parser(
        cmd_name, help=f"Query {api_name} endpoint{searchable}"
    )
    _add_common_args(sp)
    if cmd_name == "activity":
      sp.add_argument(
          "--normalize",
          action="store_true",
          help="Normalize bioactivity values to nM",
      )
    if cmd_name == "molecule":
      sp.add_argument(
          "--dl_format",
          type=str,
          choices=["sdf", "mol"],
          help="Download molecule structure file (SDF or MOL)",
      )
    sp.set_defaults(func=cmd_generic)

  sp_status = subparsers.add_parser("status", help="Check ChEMBL API status")
  sp_status.add_argument(
      "--output",
      type=str,
      required=True,
      help="Output JSON file path (required)",
  )
  sp_status.add_argument(
      "--mock",
      action="store_true",
      help="Enable offline synthetic mock execution",
  )
  sp_status.set_defaults(func=cmd_status)

  sp_sim = subparsers.add_parser(
      "similarity", help="Server-side similarity search by SMILES"
  )
  sp_sim.add_argument("--smiles", type=str, required=True, help="SMILES string")
  sp_sim.add_argument(
      "--similarity",
      type=int,
      required=True,
      help="Similarity threshold (0-100)",
  )
  sp_sim.add_argument(
      "--limit",
      type=int,
      default=5,
      help="Max results (default: 5)",
  )
  sp_sim.add_argument(
      "--offset",
      type=int,
      default=None,
      help="Pagination offset",
  )
  sp_sim.add_argument(
      "--output",
      type=str,
      required=True,
      help="Output JSON file path (required)",
  )
  sp_sim.add_argument(
      "--mock",
      action="store_true",
      help="Enable offline synthetic mock execution",
  )
  sp_sim.set_defaults(func=cmd_similarity)

  sp_sub = subparsers.add_parser(
      "substructure", help="Server-side substructure search by SMILES"
  )
  sp_sub.add_argument("--smiles", type=str, required=True, help="SMILES string")
  sp_sub.add_argument(
      "--limit",
      type=int,
      default=5,
      help="Max results (default: 5)",
  )
  sp_sub.add_argument(
      "--offset",
      type=int,
      default=None,
      help="Pagination offset",
  )
  sp_sub.add_argument(
      "--output",
      type=str,
      required=True,
      help="Output JSON file path (required)",
  )
  sp_sub.add_argument(
      "--mock",
      action="store_true",
      help="Enable offline synthetic mock execution",
  )
  sp_sub.set_defaults(func=cmd_substructure)

  sp_img = subparsers.add_parser(
      "image", help="Download compound image (SVG by default)"
  )
  sp_img.add_argument(
      "--id",
      type=str,
      required=True,
      help="ChEMBL ID or InChI Key",
  )
  sp_img.add_argument(
      "--output",
      type=str,
      required=True,
      help="Output file path",
  )
  sp_img.add_argument(
      "--dimensions",
      type=int,
      help="Image size in pixels (max 500, default 500)",
  )
  sp_img.add_argument(
      "--engine",
      type=str,
      default=None,
      help="Rendering engine (default: rdkit)",
  )
  sp_img.add_argument(
      "--img_format",
      type=str,
      choices=["svg", "png"],
      default=None,
      help="Image format: svg (default) or png",
  )
  sp_img.add_argument(
      "--mock",
      action="store_true",
      help="Enable offline synthetic mock execution",
  )
  sp_img.set_defaults(func=cmd_image)

  return parser


def main() -> None:
  main_parser = build_parser()
  main_args = main_parser.parse_args()

  if not main_args.command:
    main_parser.print_help()
    sys.exit(1)

  if getattr(main_args, "dl_format", None):
    cmd_molecule_download(main_args)
  else:
    main_args.func(main_args)


if __name__ == "__main__":
  main()
