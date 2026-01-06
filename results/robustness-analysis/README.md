# Robustness Analysis Results (Task 1)

**Objective:** To determine the optimal statistical threshold for defining the "PRG4 Rescue" gene signature from the internal bulk RNA-seq data, ensuring high-confidence matches with the external Perturb-seq library.

## Data Sources
*   **Internal Data:** `RPE_gene pvals.xlsx` (DESeq2 results for RPE cells treated with H2O2 +/- PRG4).
*   **External Data:** **RPE1 Essential Perturb-seq** (Replogle et al. 2022, `rpe1_normalized_bulk_01.h5ad`).

## Methods
1.  **Threshold Comparison:** We tested four statistical cutoffs for defining Differentially Expressed Genes (DEGs):
    *   Nominal P-value < 0.05 (Loose)
    *   Nominal P-value < 0.01
    *   FDR < 0.1
    *   **FDR < 0.05 (Strict)** - *Selected for downstream analysis.*
2.  **Validation Metric:** For each threshold, we calculated the Spearman correlation between the resulting signature and the perturbation profiles in the RPE1 dataset. Higher correlation implies reduced noise and stronger biological signal.

## Key Findings
*   **Selection:** The **FDR < 0.05** threshold yielded the highest consistency and correlation scores (Spearman rho ~0.32 vs ~0.29 for loose p-values).
*   **Impact:** Using this stricter threshold reduced the signature size to ~3,900 genes but significantly increased the reliability of the "Virtual Screen" hits.

## Files
*   **`threshold_comparison.csv`**: Metrics (N_DEGs, Correlation) for all tested thresholds.
*   **`threshold_sensitivity.pdf`**: Visualization of correlation strength vs. threshold.
*   **`jaccard_stability.pdf`**: Heatmap demonstrating that top hits (e.g., *KEAP1*) are robust across stringent thresholds.