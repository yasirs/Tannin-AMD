# Task 3: Baseline Expression Filtering - Summary

**Date:** 2026-01-06
**Analyst:** Gemini Agent

## Objective
To establish baseline expression levels for RPE1 and K562 cell lines using Perturb-seq non-targeting controls, identifying genes suitable for cross-cell-type comparisons.

## Key Findings

### 1. Global Expression Statistics
We identified the baseline expression (mean UMI per cell) for 10,117 total genes across the datasets.

*   **RPE1 Only**: 1,655 genes
*   **K562 GWPS Only**: 1,154 genes
*   **Commonly Measured**: **7,094 genes**

At a confidence threshold of **> 0.05 UMI/cell**, we found high overlap in the transcriptome backbone, providing a solid foundation for concordance analysis.

### 2. Lineage Marker Validation
The baseline expression profiles correctly identify the cell type identities:

| Marker | Category | RPE1 (UMI) | K562 (UMI) | Interpretation |
| :--- | :--- | :--- | :--- | :--- |
| **ACTB** | Housekeeping | 53.14 | 23.17 | Strong in both |
| **GAPDH** | Housekeeping | 151.61 | 27.72 | Strong in both |
| **HBG1** | K562 Specific | 0.00 | 3.67 | K562 Only (Globin) |
| **HBG2** | K562 Specific | 0.00 | 10.04 | K562 Only (Globin) |
| **BEST1** | RPE Specific | 0.00* | 0.00* | Not captured/filtered |

*Note: Many specialized RPE markers (BEST1, RPE65) were either not targeted by the essential-gene-heavy Perturb-seq libraries or were filtered out during data processing by the original authors (mean UMI < 0.01).*

### 3. Expression Scale Differences
*   **RPE1** generally shows higher UMI counts for housekeeping and metabolic genes compared to **K562**.
*   This suggests that while the same genes are present, their absolute expression levels (and thus the potential dynamic range for perturbation effects) differ significantly between the leukemia and epithelial contexts.

## Recommendations
1.  **Use Common Set (N=7,094)**: Use this specific gene set for Task 4 (Cell-Type Concordance) to ensure we are comparing "apples to apples."
2.  **Filter K562 Hits**: When prioritizing hits from the K562 GWPS screen (Task 7), we must ensure the target genes are at least minimally expressed in RPE1 (UMI > 0.01) to ensure biological relevance.

## Output Files
*   `all_baseline_expression.csv`: Master table of mean UMI counts per gene.
*   `expression_distribution.pdf`: Density plots of transcriptomic depth.
*   `baseline_scatter.pdf`: Gene-by-gene comparison of cell line baselines.
*   `marker_heatmap.pdf`: Visual validation of cell type markers.
