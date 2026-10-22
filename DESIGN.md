# Dashboard Design Specification

## 1. Layout and Typography Architecture

The interactive dashboard (`dashboard.html`) implements a structured, high-density visual layout designed for medicinal chemistry and bioactivity analysis.

| Design Dimension | Specification | Implementation Details |
| :--- | :--- | :--- |
| **Typography** | `DM Sans` and `IBM Plex Mono` | `DM Sans` is used for section headers, navigation, and labels. `IBM Plex Mono` is used for ChEMBL IDs, SMILES strings, InChIKeys, coordinates, and tabular assay values. |
| **Surface Architecture** | Flat Solid Surfaces | Solid backgrounds with `1px solid` borders (`6px` to `8px` border radius). No decorative background gradients or drop shadows. |
| **Container Hierarchy** | Flat Bento Grid | Sections are organized as flat Bento panels divided by `1px solid` hairlines without nested card-in-card containers. |
| **Status Indicators** | Exception-Only Alerts | Status badges are rendered only when Lipinski Rule of 5 violations (`violations > 0`) or low-potency threshold alerts (`> 10 uM`) occur. |
| **Chart Gridlines** | Data-Scoped Gridlines | Gridlines are restricted strictly to coordinate plotting areas (Lipinski radar axes and pIC50 histogram reference lines). |
| **Header Hierarchy** | Text-Based Section Headers | Panel and subcell headers use plain typographic hierarchy without decorative icons next to titles. |

---

## 2. Runtime and Execution Specifications

* **Self-Contained Delivery:** The compiled `dashboard.html` file embeds all CSS rules, SVG structures, 3D SDF coordinate blocks, serialized JSON payloads, and interactive JavaScript routines in a single standalone HTML5 document.
* **Sandboxed Execution Compatibility:** Client-side state management wraps `window.localStorage` access in defensive `try/catch` blocks (`safeGetStorage` and `safeSetStorage`) so theme toggling and interactive controls operate reliably inside restricted iframe sandboxes.
* **Interactive Viewport Controls:** Every visual diagram and chart viewport provides `Zoom In`, `Zoom Out`, and `Reset View` controls alongside mouse-wheel zooming and click-and-drag panning or rotation.
* **Data Integrity:** All visual components bind directly to the processed ChEMBL compound, bioactivity, mechanism, indication, and similarity datasets embedded in `#dashboard-data` and `#sdf-data`.

---

## 3. Light and Dark Theme Specifications

The dashboard natively supports both light mode and high-contrast dark mode via CSS custom properties and the `html.dark` / `data-theme` root selector.

* **Theme Switcher (`#theme-toggle-btn`):** A single minimalist button in the top application header that switches between light and dark palettes and dynamically updates all 2D/3D viewports and network graphs.
* **Color Tokens:**
  * **Light Mode:**
    * Page Background (`--bg-page`): `#FAFAF9`
    * Cards and Surfaces (`--bg-surface`): `#FFFFFF`
    * Borders and Hairlines (`--border`): `#E7E5E4` (`1px solid`, `6px` to `8px` border-radius)
    * Text Primary (`--text-primary`): `#1C1917`
    * Text Secondary (`--text-secondary`): `#787774`
  * **Dark Mode:**
    * Page Background (`--bg-page`): `#0D0E11`
    * Cards and Surfaces (`--bg-surface`): `#16171D`
    * Borders and Hairlines (`--border`): `#292A30` (`1px solid`, `6px` to `8px` border-radius)
    * Text Primary (`--text-primary`): `#ECE8E2`
    * Text Secondary (`--text-secondary`): `#9A9895`

---

## 4. Component Layout and Table Rules

* **Direct Component Titles:** Component headers use concise, descriptive titles:
  * `2D Chemical Structure`
  * `Structural Data Diagram`
  * `Physicochemical Profile (Lipinski Ro5)`
  * `Interactive Relationship Node Map`
  * `Clinical Development & Target Mechanisms`
  * `Bioactivity & Potency Distribution`
  * `Chemical Similarity & Descriptor Matrix`
  * `Normalized Descriptor Heatmap Matrix`
* **Monospace Data Tables:** All tabular views (mechanisms, clinical indications, bioactivity assays, and the descriptor heatmap matrix) use `IBM Plex Mono` for aligned numerical reading, with `1px solid` horizontal row dividers and no vertical cell borders.
* **3D Structural Data Diagram:** Renders 3D atomic coordinates and bond orders from the embedded SDF block with interactive mouse rotation, `Zoom In`, `Zoom Out`, `Reset View`, `Auto-Spin`, and synchronized active button highlighting across `Sticks`, `Ball & Stick`, `Spheres`, and `Surface` representations.

---

## 5. Data Visualization Color Palette

All charts, molecular topologies, 3D conformers, network graphs, and heatmap matrices use a semantic 5-color visualization palette designed for high contrast against both light and dark surfaces:

| Color Role | Hex Code | Visual Application |
| :--- | :--- | :--- |
| **Vibrant Teal** | `#0D9488` | Primary data series, carbon skeleton bonds, compliant radar polygons, moderate potency bins (`pIC50 7.0-7.9`), and primary analog nodes (`85-89%` similarity). |
| **Energetic Amber** | `#D97706` | Secondary data series, ester functional groups, moderate potency bins (`pIC50 5.0-6.9`), and moderate similarity analogs (`< 85%`). |
| **Vivid Indigo** | `#4F46E5` | Query compound focal node in the relationship map, methyl group highlights, and molecular weight / lipophilicity heatmap scales. |
| **Functional Emerald** | `#059669` | Lipinski Ro5 boundary polygon, high-potency assay bins (`pIC50 >= 8.0`), completed clinical phases, and high-similarity analogs (`>= 90%`). |
| **Deep Rose** | `#E11D48` | Lipinski rule violation alerts, carbonyl/hydroxyl oxygen atoms, and low-potency assay bins (`pIC50 < 5.0` or `> 10 uM`). |

### Visualization Guidelines
* **Matrices and Heatmaps:** The `Normalized Descriptor Heatmap Matrix` applies a continuous opacity scale transitioning into saturated **Vibrant Teal** (`#0D9488`) and **Vivid Indigo** (`#4F46E5`) to represent normalized Tanimoto similarity, molecular weight, lipophilicity (`AlogP`), polar surface area (`TPSA`), and hydrogen-bond counts (`HBD / HBA`).
* **Relationship Node Map and Molecular Graphs:** Nodes use distinct palette colors corresponding to similarity tiers, connected by semi-transparent edges (`rgba`) weighted by Tanimoto similarity.
* **Surface Inheritance:** Chart viewports use transparent backgrounds so they inherit `--bg-surface` dynamically across theme transitions, with faint internal gridlines (`#E7E5E4` in light mode, `#292A30` in dark mode).
