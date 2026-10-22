# Product Requirement Document (PRD)
## Project: ChEMBL Database Visualize Skill (`chembl-database-visualize`)

### Document Metadata
- Document Version: 2.0.0
- Publication Date: 2026-10-22
- License: Apache License 2.0
- Target Platform: Standard AI Agent Platforms and CLI Environments

---

### 1. User Stories and Project Scope

#### 1.1 Target Personas
1. **Computational Chemist and Pharmacologist (Computational Domain Specialist)**:
   - Needs programmatic access to precise bioactivity metrics, standardized concentrations, and chemical structures.
   - Requires automated pIC50 and pKi calculations, Lipinski drug-likeness scoring, and interactive 3D conformer exploration to prioritize candidate compounds.
2. **Clinical Researcher and Translational Analyst (Translational Domain Analyst)**:
   - Needs rapid cross-referencing between target proteins, known approved drugs, clinical trial development stages, and mechanism of action annotations.
   - Requires consolidated clinical phase timelines and disease indication mappings linked to standard vocabularies (EFO, MeSH, ATC).
3. **Agentic System and Autonomous LLM Workflow (Autonomous Consumer)**:
   - Executes multi-step discovery workflows without human intervention.
   - Requires strict CLI contract adherence, deterministic JSON outputs, offline synthetic mocking capability for sandboxed testing, and self-contained interactive HTML dashboard artifacts.

#### 1.2 In-Scope Capabilities
- Backward compatibility with the baseline ChEMBL API client covering 30 REST endpoints (molecules, targets, bioactivities, assays, mechanisms, indications, documents, and image downloads).
- Rate-limited HTTP query engine enforcing 5.0 queries per second (QPS) with exponential backoff and jitter on HTTP 429 and 503 response codes.
- Concentration unit normalization converting diverse units (M, mM, uM, pM, fM) to standardized nanomolar (nM) units.
- Calculation of negative logarithmic potency: `pIC50 = 9 - log10(normalized_value_nM)`.
- Computation of Lipinski Rule of 5 parameters: Molecular Weight (`MW <= 500 Da`), AlogP (`<= 5.0`), Hydrogen Bond Donors (`HBD <= 5`), Hydrogen Bond Acceptors (`HBA <= 10`), Rotatable Bonds (`<= 10`), and Topological Polar Surface Area (`TPSA <= 140 A^2`), with status alerts displayed only when rule violations occur (`violations > 0`).
- Pure SVG Lipinski Rule of 5 radar chart (`#radar-svg-container`) with `Zoom In`, `Zoom Out`, and `Reset View` controls alongside mouse-wheel zooming and drag panning, paired with a flat 6-cell hairline-divided descriptor matrix.
- Interactive 2D chemical topology visualization (`#chem-svg-container`) with element-specific color coding (carbonyl/hydroxyl oxygen `#E11D48`, ester oxygen `#D97706`, methyl CH3 `#4F46E5`, aromatic carbon skeleton `#0D9488`), `Zoom In`, `Zoom Out`, `Reset View`, and `Toggle Highlights` controls.
- Interactive 3D molecular conformer viewer (`Structural Data Diagram`, `#mol-3d-viewer`) parsing embedded SDF atomic coordinates and bond orders with interactive mouse drag rotation, `Zoom In`, `Zoom Out`, `Reset View`, `Auto-Spin` toggle, and synchronized active button highlighting across four rendering styles (`Sticks`, `Ball & Stick`, `Spheres`, and translucent VDW `Surface`).
- Interactive force-directed relationship graph (`Interactive Relationship Node Map`, `#physics-viewport`) representing query compound (`#4F46E5`) and structural analog nodes (`#059669`, `#0D9488`, `#D97706`) with spring-repulsion simulation, node dragging, hover tooltips, `Zoom In`, `Zoom Out`, `Reset View`, `Perturb`, and `Pause` controls.
- High-contrast light and dark theme architecture (`#theme-toggle-btn`) with defensive storage access (`safeGetStorage` / `safeSetStorage`) for sandboxed iframe compatibility and synchronized CSS custom property tokens across all surfaces, tables, and viewports.
- Chemical similarity explorer featuring a high-contrast `Normalized Descriptor Heatmap Matrix` (interpolating Vivid Indigo `#4F46E5` and Vibrant Teal `#0D9488` across Tanimoto similarity, molecular weight, AlogP, TPSA, and HBD/HBA), hairline-divided analog cards, and real-time client-side sliders (`Min Similarity` and `Max MW`) with a `Reset` action.
- End-to-end unified workflow pipeline (`scripts/run_pipeline.py`) executing compound, target, or SMILES searches, automatically integrating chemical similarity analogs into compound profiles, and compiling standalone `dashboard.html` files using relative workspace paths.
- Deterministic offline mock mode (`--mock`) utilizing synthetic data without external network egress.
- License verification guard ensuring compliance with EMBL-EBI terms and CC BY-SA licensing terms.

#### 1.3 Out-of-Scope Capabilities
- De novo molecular generation or generative chemistry optimization.
- Molecular dynamics simulations, quantum mechanics property calculations, or binding free energy simulations (FEP/TI).
- Direct modification or write-back mutations to the upstream EMBL-EBI ChEMBL databases.
- Integration of proprietary, non-public chemical libraries without explicit adapter configuration.

---

### 2. Technical Architecture

#### 2.1 System Component Interaction
The following diagram details the interaction between the host agent harness, `SKILL.md`, the Python command-line utilities, the mock engine, and generated visual artifacts:

```
+-------------------------------------------------------------------------+
|                              Agent Harness                              |
|                   (Standard Agent Harness, Headless CLI)                |
+-------------------------------------------------------------------------+
                                     |
                                     v
+-------------------------------------------------------------------------+
|                                SKILL.md                                 |
|            (System Prompts, Tool Registry, Execution Rules)             |
+-------------------------------------------------------------------------+
                                     |
                    +----------------+----------------+
                    |                                 |
                    v                                 v
+---------------------------------------+ +-------------------------------+
|          scripts/run_pipeline.py      | |     scripts/chembl_api.py     |
|   (Unified Orchestrator & Dispatch)   | | (Rate-Limited REST Client)    |
+---------------------------------------+ +-------------------------------+
       |                    |                            |
       | License Check      | Pipeline Calls             | Network / Mock
       v                    v                            v
+--------------+   +-----------------------+   +--------------------------+
| .licenses/   |   | scripts/              |   | EMBL-EBI ChEMBL REST API |
| LICENSE.txt  |   | visualize_compound.py |   | (or sample_data/ mock)   |
+--------------+   | visualize_bioactivity |   +--------------------------+
                   | visualize_similarity  |
                   +-----------------------+
                               |
                               v
                   +-----------------------+
                   | scripts/              |
                   | generate_dashboard.py |
                   +-----------------------+
                               |
                               v
                   +-----------------------+
                   |    dashboard.html     |
                   | (Self-Contained HTML5 |
                   |  Bento Analytics App) |
                   +-----------------------+
```

#### 2.2 Sequence Diagram: Compound Analysis Workflow

```
Agent Harness           run_pipeline.py        chembl_api.py       generate_dashboard.py
     |                        |                      |                        |
     |--- run_pipeline.py --->|                      |                        |
     |    (--molecule CHEMBL25|                      |                        |
     |     -o ./output)       |                      |                        |
     |                        |-- verify license --->|                        |
     |                        |   guard file         |                        |
     |                        |                      |                        |
     |                        |-- fetch molecule --->|                        |
     |                        |   metadata           |                        |
     |                        |<-- return JSON ------|                        |
     |                        |                      |                        |
     |                        |-- download SDF ----->|                        |
     |                        |<-- return .sdf ------|                        |
     |                        |                      |                        |
     |                        |-- download SVG ----->|                        |
     |                        |<-- return .svg ------|                        |
     |                        |                      |                        |
     |                        |-- fetch bioactivity->|                        |
     |                        |<-- return norm nM ---|                        |
     |                        |                      |                        |
     |                        |-- compile data ------------------------------>|
     |                        |   (properties, SAR, structures, analogs)      |
     |                        |                                               |
     |                        |<-- write self-contained dashboard.html -------|
     |                        |                                               |
     |<-- stdout summary -----|                                               |
     |    (exit code 0)       |                                               |
```

---

### 3. Functional Requirements

#### 3.1 Input Specifications
- Primary Query Parameters (Mutually Exclusive):
  - `--molecule <CHEMBL_ID>`: Alphanumeric ChEMBL compound identifier (e.g., `CHEMBL25`).
  - `--target <CHEMBL_ID>`: Alphanumeric ChEMBL target identifier (e.g., `CHEMBL203`).
  - `--smiles <STRING>`: Valid SMILES string for structural comparison (e.g., `CC(=O)Oc1ccccc1C(=O)O`).
  - `--search <QUERY>`: Free-text search phrase for target or molecule resolution.
- Modifiers:
  - `--similarity <INT>`: Integer threshold between 0 and 100 (default: 85).
  - `--activity_type <STRING>`: Bioactivity filter (e.g., `IC50`, `Ki`, `EC50`, `Kd`).
  - `--limit <INT>`: Upper bound for fetched records (default: 10).
  - `--normalize`: Boolean flag enforcing conversion to nanomolar units.
  - `--mock`: Boolean flag enabling local synthetic execution.
  - `-o / --output_dir <PATH>`: Relative or local destination directory for all generated artifacts.

#### 3.2 Analytical Processing Stages
1. **License Verification**:
   - The pipeline checks for `.licenses/chembl_database_visualize_LICENSE.txt` in the workspace root.
   - If missing, it prints the EMBL-EBI terms notification and initializes the file with the current ISO 8601 UTC timestamp.
2. **Data Acquisition**:
   - Dispatches requests to `scripts/chembl_api.py`.
   - In live mode: applies 5.0 QPS rate limiting with jitter backoff.
   - In mock mode: bypasses sockets and extracts matching payloads from `sample_data/`.
3. **Physicochemical Evaluation**:
   - Parses `full_mwt`, `alogp`, `hbd`, `hba`, `psa`, and `rtb` from compound properties.
   - Calculates Lipinski Rule of 5 violations and renders alert badges strictly when `violations > 0`.
4. **Bioactivity Standardization & SAR**:
   - Normalizes concentration units into standard nM.
   - Computes `pIC50 = 9 - log10(nM)`.
   - Groups records into high (`< 100 nM`), moderate (`100` to `10,000 nM`), and weak (`> 10,000 nM`) potency bins.
5. **Dashboard Compilation & Component Integration**:
   - Compiles five synchronized flat Bento visual analytics sections into a single self-contained HTML5 application:
     - **Compound Overview Panel (`#compound-overview-section`)**:
       - `2D Chemical Structure` (`#chem-svg-container`) with non-overlapping element-colored topology, `Zoom In`, `Zoom Out`, `Reset View`, and `Toggle Highlights`.
       - `Structural Data Diagram` (`#mol-3d-viewer`) rendering 3D atomic conformers from embedded SDF coordinates with `Zoom In`, `Zoom Out`, `Reset View`, `Auto-Spin`, and synchronized active button highlighting across `Sticks`, `Ball & Stick`, `Spheres`, and `Surface` modes.
       - `Physicochemical Profile (Lipinski Ro5)` (`#radar-svg-container`) 6-axis SVG radar chart with `Zoom In`, `Zoom Out`, and `Reset View` controls, paired with a 6-cell hairline-divided descriptor matrix (`MW`, `AlogP`, `HBD`, `HBA`, `Rotatable Bonds`, `TPSA`) and SMILES/InChIKey identifier strip.
     - **Interactive Relationship Node Map (`#relationship-map-section`)**:
       - Force-directed SVG graph (`#physics-viewport`) modeling compound-analog structural relationships with node dragging, hover tooltips, `Zoom In`, `Zoom Out`, `Reset View`, `Perturb`, and `Pause` controls.
     - **Clinical Development & Target Mechanisms (`#clinical-section`)**:
       - 5-stage milestone timeline tracking regulatory progression from Preclinical (Phase 0) through Approved (Phase 4) alongside horizontal-divided monospace tables for mechanisms of action and clinical indications.
     - **Bioactivity & Potency Distribution (`#bioactivity-section`)**:
       - 4-metric summary strip (`Assays Screened`, `Median pIC50`, `Peak Potency`, `Sub-100 nM Hits`), standardized nanomolar pIC50 histogram (`#hist-svg-container`) with `Zoom In`, `Zoom Out`, and `Reset View` controls, and a real-time search-filterable monospace assay table (`#activity-table`) with a `Reset` button.
     - **Chemical Similarity & Descriptor Matrix (`#similarity-section`)**:
       - High-contrast `Normalized Descriptor Heatmap Matrix` and hairline-divided analog cards controlled simultaneously by real-time `Min Similarity` and `Max MW` range sliders and a `Reset` button.
   - Embeds a single minimalist light/dark theme switcher (`#theme-toggle-btn`) with synchronized CSS token updates across all DOM elements, charts, and 3D viewports.

#### 3.3 Dashboard Design System and Visual Specifications
- **Typography**: `DM Sans` for UI headers, labels, and navigation; `IBM Plex Mono` for ChEMBL IDs, SMILES, InChIKeys, coordinates, and tabular data.
- **Light Mode Tokens**:
  - Page Background (`--bg-page`): `#FAFAF9`
  - Cards and Surfaces (`--bg-surface`): `#FFFFFF`
  - Borders and Hairlines (`--border`): `#E7E5E4` (`1px solid`, `6px` to `8px` border-radius)
  - Primary Text (`--text-primary`): `#1C1917`
  - Secondary Text (`--text-secondary`): `#787774`
- **Dark Mode Tokens**:
  - Page Background (`--bg-page`): `#0D0E11`
  - Cards and Surfaces (`--bg-surface`): `#16171D`
  - Borders and Hairlines (`--border`): `#292A30` (`1px solid`, `6px` to `8px` border-radius)
  - Primary Text (`--text-primary`): `#ECE8E2`
  - Secondary Text (`--text-secondary`): `#9A9895`
- **Data Visualization Palette**:
  - Vibrant Teal (`#0D9488`): Primary data series, carbon skeleton bonds, compliant radar polygons, and primary analog nodes.
  - Energetic Amber (`#D97706`): Secondary data series, ester functional groups, and moderate potency/similarity tiers.
  - Vivid Indigo (`#4F46E5`): Query compound focal node, methyl group highlights, and molecular weight / lipophilicity heatmap scales.
  - Functional Emerald (`#059669`): Lipinski Ro5 boundary polygon, high-potency assay bins (`pIC50 >= 8.0`), completed clinical phases, and high-similarity analogs (`>= 90%`).
  - Deep Rose (`#E11D48`): Lipinski rule violation alerts, carbonyl/hydroxyl oxygen atoms, and low-potency threshold alerts (`> 10 uM`).

#### 3.4 Output Specifications
The output directory must contain:
- `data.json`: Complete serialized ChEMBL records and analytical results.
- `structure.sdf`: 3D coordinates file (when available).
- `structure.svg`: 2D vector chemical graphic (when available).
- `dashboard.html`: Self-contained, responsive single-page web application.

---

### 4. Non-Functional Requirements

#### 4.1 Latency and Performance
- Under `--mock` mode, the complete pipeline execution must complete in under 5.0 seconds.
- In live mode, requests must strictly respect the 5.0 QPS ceiling (200 ms minimum inter-request interval).
- Generated `dashboard.html` files must render in standard modern web browsers in under 1.5 seconds without layout shifts.

#### 4.2 Security, Privacy, and Healthcare Compliance
- **Zero PHI/PII Mandate**: All code, test fixtures, and mock files must contain exclusively synthetic or public reference chemistry data. No patient-level or clinical research participant data is permitted.
- **No Credential Exposure**: The skill must not require, store, or log API keys or internal credentials.
- **Relative Path Portability**: All CLI summaries, documentation references, and script resolution paths use relative workspace paths rooted at the repository top level.
- **Bounded Network Egress**: All network traffic is bounded to `https://www.ebi.ac.uk/chembl/api/data` and standard public CDNs (`cdn.jsdelivr.net`, `cdnjs.cloudflare.com`).

#### 4.3 Reliability and Extensibility
- Dual-engine HTTP support: Automatic fallback to standard library `urllib.request` when `polite-http` is unavailable.
- Defensive structure and storage handling: When 3D SDF coordinates are missing or when `window.localStorage` is restricted inside a sandboxed iframe, the dashboard degrades gracefully without throwing uncaught JavaScript exceptions.
- Modular code architecture: Separate scripts for compound profiling, bioactivity analysis, similarity filtering, and dashboard generation permit independent invocation or integration into larger computational pipelines.

---

### 5. Error Handling and Edge Cases Matrix

| Error Condition | Trigger Event | System Response | Output State |
| :--- | :--- | :--- | :--- |
| Missing License Guard | Workspace root lacks `.licenses/chembl_database_visualize_LICENSE.txt` | Displays licensing notice; writes guard file with ISO 8601 UTC timestamp | Execution proceeds normally |
| Invalid Molecule ID | ChEMBL ID format invalid or non-existent in database | Returns HTTP 404 or structured error payload | Exit code 1; JSON error message written to stderr/stdout |
| Missing 3D Conformer | Molecule has 2D topology but no external 3D SDF file | Pipeline embeds fallback 3D coordinates or flags `has_3d=false` | `dashboard.html` renders interactive 3D projection without exceptions |
| Sandboxed Iframe Storage Restriction | Host iframe blocks `window.localStorage` (`SecurityError`) | `safeGetStorage` and `safeSetStorage` catch storage exceptions | Theme toggle and UI controls operate normally in memory |
| Non-standard Bioactivity Units | Assay reports qualitative or atypical units (e.g., `% inhibition`, `AU`) | Normalization engine sets `normalized_value_nM=null` and `pIC50=null` | Records preserved in table; omitted from pIC50 histogram |
| Upstream Rate Limit (HTTP 429) | Exceeded EMBL-EBI concurrency or quota thresholds | Exponential backoff retry with random jitter (up to 3 retries) | Retries transparently; fails with error JSON if exhausted |
| Upstream Outage (HTTP 503) | EMBL-EBI server maintenance or transient failure | Exponential backoff retry with random jitter (up to 3 retries) | Retries transparently; fails with error JSON if exhausted |
| Network Isolation / Sandboxed CI | Pipeline invoked in environment without external internet | Pipeline fails unless `--mock` is specified | `--mock` uses local synthetic fixtures; exit code 0 |
| Empty Similarity Search | SMILES search yields no analogs above threshold | Pipeline returns empty similarity array | Dashboard displays notice: "No chemical analogs found for the specified query." |
