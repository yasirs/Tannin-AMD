# Baseline Expression Results (Task 3)

**Objective:** To establish which genes are actively expressed in the cell lines used for screening, ensuring we don't predict targets that are irrelevant (not expressed) in the RPE context.

## Data Sources
*   **Raw Perturb-seq Data:** `rpe1_raw_bulk_01.h5ad` and `K562_gwps_raw_bulk_01.h5ad` (Replogle et al. 2022). Note: These "raw" files contain mean UMI counts per cell.

## Methods
1.  **Baseline Calculation:** We calculated the mean expression (UMI/cell) of all **Non-Targeting Control** populations in each dataset.
2.  **Expression Threshold:** A gene was considered "Expressed" if its mean UMI count was **> 0.05**.
3.  **Filtering:** This "Expressed in RPE1" list was used to filter the final Virtual Screen hits.

## Key Findings
*   **Common Core:** We identified **7,094 genes** that are robustly expressed in both RPE1 and K562.
*   **Marker Validation:** Confirmed that *ACTB* and *GAPDH* are high in both, while *HBG1/2* (Hemoglobin) are exclusively in K562.

## Files
*   **`expressed_genes_RPE1.csv`**: The master list of druggable RPE targets.
*   **`all_baseline_expression.csv`**: Full table of mean UMI counts for all genes.
*   **`marker_heatmap.pdf`**: Visual validation of cell-type markers.
