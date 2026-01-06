# Single-Cell Analysis Results

**Objective:** To infer tissue-level effects by projecting bulk signatures onto cell-type specific markers derived from single-cell literature.

## Data Sources
*   **Markers:** Derived from **Voigt et al. (2019)** "Single-cell transcriptomics of the human retinal pigment epithelium and choroid".
*   **Signatures:** Internal H2O2 and PRG4 Rescue DEGs.

## Methods
*   **Enrichment:** Fisher's Exact Test to determine if "Rescue" genes are enriched for specific cell-type markers (RPE, Endothelial, Fibroblast).

## Key Findings
*   **RPE Identity:** PRG4 rescue is associated with the upregulation of key RPE identity markers (**MITF, TYR, SERPINF1**), suggesting it prevents dedifferentiation (EMT).
*   **Fibrosis:** PRG4 rescue is associated with the downregulation of fibroblast markers (**COL1A1, FBLN1**), indicating an anti-fibrotic effect.
*   **Receptors:** Confirmed that **CD44** and **TLR2** (PRG4 receptors) are expressed and modulated in the bulk data.

## Files
*   **`cell_type_enrichment_heatmap.png`**: Visual summary of marker enrichment.
*   **`cell_type_enrichment.csv`**: Statistical results (p-values) for each cell type.