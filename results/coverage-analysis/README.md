# Coverage Analysis Results (Task 2)

**Objective:** To quantify the limitations of using "Essential Gene" screens and justify the use of the cross-cell-type K562 Genome-Wide dataset.

## Data Sources
All Perturb-seq datasets are from **Replogle et al. (2022)**:
1.  **RPE1 Essential:** `rpe1_normalized_bulk_01.h5ad` (~2,500 genes).
2.  **K562 Essential:** `K562_essential_normalized_bulk_01.h5ad` (~2,000 genes).
3.  **K562 Genome-Wide (GWPS):** `K562_gwps_normalized_bulk_01.h5ad` (~11,000 genes).

## Methods
*   **Query Set:** The **PRG4 Rescue Signature** (DEGs with FDR < 0.05).
*   **Metric:** Percent Coverage = (Number of DEGs present in Perturb-seq Library) / (Total DEGs).

## Key Findings
*   **Critical Gap:** The **RPE1 Essential** library covers only **11-16%** of the genes in our AMD signature. Crucial pathways (e.g., GPCR signaling, ECM remodeling) are completely absent.
*   **Solution:** The **K562 GWPS** library covers **~65%** of the signature.
*   **Decision:** We must use the K562 GWPS dataset for the primary Virtual Screen to avoid missing >80% of the relevant biology.

## Files
*   **`coverage_comparison.pdf`**: Bar chart contrasting the coverage of the three libraries.
*   **`enrichment_missing_h2o2.pdf`**: Enrichment analysis of the genes *missed* by the RPE1 Essential screen (proving we would miss important biology).
*   **`coverage_summary.csv`**: Raw counts of measured and perturbed genes.