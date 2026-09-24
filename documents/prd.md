# Product Requirement Document (PRD)
## Project: ChEMBL Database Visualize Skill (`chembl-database-visualize`)

### Document Metadata
- Document Version: 1.2.0
- Publication Date: 2026-09-18
- License: Apache License 2.0
- Target Platform: Standard AI Agent Platforms and CLI Environments

---

### 1. User Stories and Project Scope

#### 1.1 Target Personas
1. **Computational Chemist and Pharmacologist (Technical User)**:
   - Needs programmatic access to precise bioactivity metrics, standardized concentrations, and chemical structures.
   - Requires automated pIC50 and pKi calculations, Lipinski drug-likeness scoring, and interactive 3D conformer exploration to prioritize candidate compounds.
2. **Clinical Researcher and Translational Analyst (Semi-Technical User)**:
   - Needs rapid cross-referencing between target proteins, known approved drugs, clinical trial development stages, and mechanism of action annotations.
   - Requires consolidated clinical phase timelines and disease indication mappings linked to standard vocabularies (EFO, MeSH, ATC).
3. **Agentic System and Autonomous LLM Workflow (Autonomous Consumer)**:
   - Executes multi-step discovery workflows without human intervention.
   - Requires strict CLI contract adherence, deterministic JSON outputs, offline synthetic mocking capability for sandboxed testing, and self-contained interactive HTML dashboard artifacts.

#### 1.2 In-Scope Capabilities
- Backward compatibility with the baseline ChEMBL API client covering 30 REST endpoints (molecules, targets, bioactivities, assays, mechanisms, indications, documents, and image downloads).
- Rate-limited HTTP query engine enforcing 5.0 queries per second (QPS) with exponential backoff and jitter on HTTP 429 and 503 response codes.
- Concentration unit normalization converting diverse units (M, mM, uM, uM, pM, fM) to standardized nanomolar (nM) units.
- Calculation of negative logarithmic potency: pIC50 = 9 - log10(normalized_value_nM).
- Computation of Lipinski Rule of 5 parameters: Molecular Weight (MW <= 500 Da), AlogP (<= 5.0), Hydrogen Bond Donors (HBD <= 5), Hydrogen Bond Acceptors (HBA <= 10), Rotatable Bonds (<= 10), and Topological Polar Surface Area (TPSA <= 140 A^2).
- Pure SVG Spider/Radar chart rendering depicting chemical space boundaries and compound property compliance.
- Interactive 2D chemical topology visualization with element-specific color coding (Oxygen red, Ester Oxygen amber, Methyl CH3 indigo, Carbon skeleton dark slate/white), non-overlapping skeletal geometry, interactive "Toggle Highlights" button for functional group halos and ring shading, and SVG tooltip elements.
- Interactive 3D molecular conformer viewer using 3Dmol.js rendering embedded SDF coordinates with automated interactive HTML5 3D rotational canvas fallback, Auto-Spin rotational physics, Reset View, and rendering style selectors (Sticks, Ball & Stick, Spheres, Surface).
- Interactive physics-driven force-directed graph modeling representing compound, target, and analog nodes with spring forces, charge repulsion, drag-and-drop interactions, and Perturb, Pause, and Recenter controls.
- Minimalist icon-only light and dark theme toggle featuring Sun and Moon SVG icons with class-based Tailwind styling and CSS variable overrides across all components.
- Tanimoto similarity filtering grid with responsive analog cards, match percentage progress bars, and client-side interactive sliders for real-time similarity and molecular weight thresholding.
- End-to-end unified workflow pipeline (`scripts/run_pipeline.py`) executing compound, target, or SMILES searches, automatically integrating chemical similarity analogs into compound profiles, and compiling standalone `dashboard.html` files.
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
                   |  Dashboard / 3Dmol)   |
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
     |                        |   (properties, SAR, structures)               |
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
  - `-o / --output_dir <PATH>`: Destination directory for all generated artifacts.

#### 3.2 Analytical Processing Stages
1. **License Verification**:
   - The pipeline checks for `.licenses/chembl_database_visualize_LICENSE.txt` in the workspace root.
   - If missing, it prints the EMBL-EBI terms notification and initializes the file with the current ISO 8601 UTC timestamp.
2. **Data Acquisition**:
   - Dispatches requests to `chembl_api.py`.
   - In live mode: applies 5.0 QPS rate limiting with jitter backoff.
   - In mock mode: bypasses sockets and extracts matching payloads from `sample_data/`.
3. **Physicochemical Evaluation**:
   - Parses `full_mwt`, `alogp`, `hbd`, `hba`, `psa`, and `rtb` from compound properties.
   - Calculates Lipinski Rule of 5 violations and assigns compliance badges.
4. **Bioactivity Standardization & SAR**:
   - Normalizes concentration units into standard nM.
   - Computes `pIC50 = 9 - log10(nM)`.
   - Groups records into high (<100 nM), moderate (100 to 10,000 nM), and weak (>10,000 nM) potency bins.
5. **Dashboard Compilation & Component Integration**:
   - Compiles five synchronized visual analytics sections into a single self-contained HTML5 application:
     - **Compound Card**: 2D chemical topology with non-overlapping geometry and highlight toggle, 3D interactive conformer with WebGL/HTML5 canvas fallback and spin/reset/style controls, SVG Lipinski radar chart, and 6-property metric grid.
     - **Force-Directed Physics Graph**: Interactive HTML5 canvas spring-mass simulation modeling compound-analog structural relationships with drag-and-drop mechanics and simulation controls.
     - **Clinical Tracker**: Milestone timeline tracking regulatory progression from Phase 0 through Phase 4 with mechanism of action and clinical indication annotations.
     - **Bioactivity Analytics**: Standardized nanomolar potency histogram and real-time search-filterable assay records table with pIC50 scores.
     - **Chemical Similarity Grid**: Interactive analog cards with Tanimoto similarity progress bars and client-side dual sliders for real-time similarity and molecular weight thresholding.
   - Embeds a minimalist icon-only light/dark theme switch with synchronized palette styling across all DOM elements, charts, and 3D viewports.

#### 3.3 Output Specifications
The output directory must contain:
- `data.json`: Complete serialized ChEMBL records and analytical results.
- `structure.sdf`: 3D coordinates file (when available).
- `structure.svg`: 2D vector chemical graphic (when available).
- `dashboard.html`: Fully self-contained, responsive single-page web application.

---

### 4. Non-Functional Requirements

#### 4.1 Latency and Performance
- Under `--mock` mode, the complete pipeline execution must complete in under 5.0 seconds.
- In live mode, requests must strictly respect the 5.0 QPS ceiling (200 ms minimum inter-request interval).
- Generated `dashboard.html` files must render in standard modern web browsers (Chrome, Firefox, Safari) in under 1.5 seconds without layout shifts.

#### 4.2 Security, Privacy, and Healthcare Compliance
- **Zero PHI/PII Mandate**: All code, test fixtures, and mock files must contain exclusively synthetic or public reference chemistry data. No patient-level or clinical research participant data is permitted.
- **No Credential Exposure**: The skill must not require, store, or log API keys or internal credentials.
- **No Unsanctioned Egress**: All network traffic is strictly bounded to `https://www.ebi.ac.uk/chembl/api/data` and approved public CDNs (Tailwind CDN, cdnjs 3Dmol.js).

#### 4.3 Reliability and Extensibility
- Dual-engine HTTP support: Automatic fallback to standard library `urllib.request` when `polite-http` is unavailable.
- Defensive structure parsing: When 3D SDF coordinates are missing, the dashboard gracefully disables the 3D WebGL viewport and displays the 2D vector structure without throwing JavaScript exceptions.
- Modular code architecture: Separate scripts for compound profiling, bioactivity analysis, similarity filtering, and dashboard generation permit independent invocation or integration into larger computational pipelines.

---

### 5. Error Handling and Edge Cases Matrix

| Error Condition | Trigger Event | System Response | Output State |
| :--- | :--- | :--- | :--- |
| Missing License Guard | Workspace root lacks `.licenses/chembl_database_visualize_LICENSE.txt` | Displays licensing notice; writes guard file with ISO 8601 UTC timestamp | Execution proceeds normally |
| Invalid Molecule ID | ChEMBL ID format invalid or non-existent in database | Returns HTTP 404 or structured error payload | Exit code 1; JSON error message written to stderr/stdout |
| Missing 3D Conformer | Molecule has 2D topology but no 3D SDF coordinates | Pipeline flags `has_3d=false`; embeds fallback placeholder | `dashboard.html` displays 2D structure and disables 3D viewport |
| Non-standard Bioactivity Units | Assay reports qualitative or atypical units (e.g., `% inhibition`, `AU`) | Normalization engine sets `normalized_value_nM=null` and `pIC50=null` | Records preserved in table; omitted from pIC50 histogram |
| Upstream Rate Limit (HTTP 429) | Exceeded EMBL-EBI concurrency or quota thresholds | Exponential backoff retry with random jitter (up to 3 retries) | Retries transparently; fails with error JSON if exhausted |
| Upstream Outage (HTTP 503) | EMBL-EBI server maintenance or transient failure | Exponential backoff retry with random jitter (up to 3 retries) | Retries transparently; fails with error JSON if exhausted |
| Network Isolation / Sandboxed CI | Pipeline invoked in environment without external internet | Pipeline fails unless `--mock` is specified | `--mock` uses local synthetic fixtures; exit code 0 |
| Empty Similarity Search | SMILES search yields no analogs above threshold | Pipeline returns empty similarity array | Dashboard displays notice: "No analogs found matching criteria" |
