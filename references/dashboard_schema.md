# Dashboard Data Contract Specification

This document details the JSON schema and data structures expected by `scripts/generate_dashboard.py` to render the interactive web dashboard (`dashboard.html`).

## Schema Overview

The dashboard payload is a structured JSON object embedded within the compiled HTML file. It aggregates molecule metadata, structural coordinate files, normalized bioactivity assays, clinical indications, mechanisms of action, and analog similarity records.

```json
{
  "metadata": {
    "title": "Compound Analysis: CHEMBL25 (Aspirin)",
    "generated_at": "2026-09-18T15:00:00Z",
    "chembl_release": "ChEMBL_34",
    "skill_version": "1.0.0",
    "disclaimer": "ChEMBL bioactivity values, molecular properties, and predicted data are aggregated from literature and high-throughput screens and must be experimentally validated before medicinal or therapeutic use."
  },
  "molecule": {
    "molecule_chembl_id": "CHEMBL25",
    "pref_name": "ASPIRIN",
    "synonyms": ["Acetylsalicylic acid", "Ecotrin", "Empirin"],
    "molecule_type": "Small molecule",
    "max_phase": 4,
    "first_approval": 1899,
    "structure": {
      "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    },
    "properties": {
      "full_mwt": 180.16,
      "alogp": 1.31,
      "hbd": 1,
      "hba": 3,
      "psa": 63.6,
      "rtb": 2,
      "ro5_violations": 0,
      "aromatic_rings": 1,
      "heavy_atoms": 13
    }
  },
  "structures": {
    "has_3d": true,
    "sdf": "...",
    "svg": "..."
  },
  "bioactivities": [
    {
      "activity_id": 10001,
      "assay_chembl_id": "CHEMBL829598",
      "target_chembl_id": "CHEMBL203",
      "target_pref_name": "Epidermal growth factor receptor erbB1",
      "target_organism": "Homo sapiens",
      "standard_type": "IC50",
      "standard_value": 45.0,
      "standard_units": "nM",
      "normalized_value_nM": 45.0,
      "pIC50": 7.35,
      "relation": "=",
      "document_chembl_id": "CHEMBL1137932"
    }
  ],
  "mechanisms": [
    {
      "mechanism_of_action": "Cyclooxygenase-1 inhibitor",
      "action_type": "INHIBITOR",
      "target_chembl_id": "CHEMBL221",
      "target_name": "Cyclooxygenase-1",
      "direct_interaction": true,
      "disease_efficacy": true
    }
  ],
  "indications": [
    {
      "efo_term": "fever",
      "mesh_heading": "Fever",
      "max_phase_for_ind": 4,
      "atc_code": "N02BA01"
    }
  ],
  "similarities": [
    {
      "molecule_chembl_id": "CHEMBL26",
      "pref_name": "SALICYLIC ACID",
      "similarity": 89.5,
      "canonical_smiles": "O=C(O)c1ccccc1O",
      "full_mwt": 138.12,
      "alogp": 1.48,
      "max_phase": 4
    }
  ]
}
```

## Field Definitions

### 1. Metadata Object
- `title` (string, required): Title displayed in the primary dashboard navigation header.
- `generated_at` (string, ISO 8601 UTC): Timestamp of generation.
- `chembl_release` (string, optional): ChEMBL data release version (for example, `ChEMBL_34`).
- `skill_version` (string, required): Semantic version of the `chembl-database-visualize` skill.
- `disclaimer` (string, required): Mandatory scientific disclaimer statement.

### 2. Molecule Object
- `molecule_chembl_id` (string, required): Primary identifier for the compound.
- `pref_name` (string, nullable): Preferred chemical or generic drug name.
- `synonyms` (array of strings): Common brand names or chemical synonyms.
- `molecule_type` (string): Chemical classification (for example, `Small molecule`, `Oligosaccharide`).
- `max_phase` (integer, 0 to 4): Maximum clinical development phase achieved.
- `properties` (object):
  - `full_mwt` (float): Molecular weight in Daltons (Ro5 threshold <= 500.0).
  - `alogp` (float): Octanol-water partition coefficient (Ro5 threshold <= 5.0).
  - `hbd` (integer): Hydrogen bond donors (Ro5 threshold <= 5).
  - `hba` (integer): Hydrogen bond acceptors (Ro5 threshold <= 10).
  - `psa` (float): Topological polar surface area in square Angstroms (threshold <= 140.0).
  - `rtb` (integer): Rotatable bond count (threshold <= 10).
  - `ro5_violations` (integer, 0 to 4): Count of Lipinski rule violations.

### 3. Structures Object
- `has_3d` (boolean): Indicates whether 3D conformer coordinates are present.
- `sdf` (string): Raw MDL/SDF string for 3D conformer WebGL rendering.
- `svg` (string): Raw SVG markup for 2D chemical structure rendering.

### 4. Bioactivities Array
- `activity_id` (integer): Unique ChEMBL activity record identifier.
- `target_chembl_id` (string): ChEMBL ID of the biological target.
- `target_pref_name` (string): Name of the target protein or complex.
- `standard_type` (string): Measurement type (`IC50`, `Ki`, `EC50`, `Kd`).
- `normalized_value_nM` (float): Standardized bioactivity concentration in nanomolar.
- `pIC50` (float): Negative logarithm of molar affinity, computed as `9 - log10(normalized_value_nM)`.

### 5. Similarities Array
- `molecule_chembl_id` (string): ChEMBL identifier of the structural analog.
- `similarity` (float, 0.0 to 100.0): Tanimoto similarity percentage score.
- `canonical_smiles` (string): SMILES chemical representation of the analog.
- `full_mwt` (float): Molecular weight of the analog.
- `alogp` (float): Partition coefficient of the analog.
