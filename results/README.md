# Results Directory - Tannin-AMD Project

This directory contains the outputs of the comprehensive computational analysis performed on January 6, 2026.

## Main Report
*   **`Tannin_AMD_Analysis_Report.md`**: The master document summarizing all findings, figures, and recommendations.

## Sub-Analysis Folders

### 1. `robustness-analysis/` (Task 1)
*   Figures showing how signature thresholds (p < 0.05 vs FDR < 0.05) affect correlation strength.
*   Enrichment plots for the defined signatures.

### 2. `coverage-analysis/` (Task 2)
*   **`coverage_comparison.pdf`**: Bar chart proving the superiority of the K562 GWPS dataset over Essential screens.
*   Data on genes missing from the Essential screen (e.g., GPCRs).

### 3. `baseline-expression/` (Task 3)
*   **`expression_distribution.pdf`**: Comparison of RPE1 and K562 transcriptome depth.
*   **`marker_heatmap.pdf`**: Validation of cell-type markers (ACTB, Globins, etc.).

### 4. `concordance-analysis/` (Task 4)
*   **`concordance_distribution.pdf`**: Histogram of correlation coefficients between RPE1 and K562 knockdowns.
*   **`gene_concordance.csv`**: Full table of correlations for ~7000 genes.

### 5. `virtual-screen/` (Task 7)
*   **`prg4_virtual_screen_results.csv`**: The ranked list of ~11,000 potential targets.
*   **`enrichment_top_mimetics.csv`**: Pathway analysis of the top hits (shows NRF2 link).

### 6. `cohort-GSE135092/` (Task 8 - RNAseq)
*   **`gse135092_pca.pdf`**: Visualizing the 537-sample cohort structure.
*   **`gse135092_risk_violins.pdf`**: Risk gene expression in AMD vs Control.

### 7. `cohort-GSE29801/` (Task 8 - Microarray)
*   **`gse29801_pca.pdf`**: PCA of the validation cohort (293 samples).
*   **`gse29801_risk_violins.pdf`**: Validation of risk gene downregulation.

### 8. `external-validation/` (Other Validation)
*   **`amd_rescue_scatter.png`**: Scatter plot showing PRG4 reversal of known AMD risk genes.
*   **`serum_rescue_scatter.pdf`**: Plot showing PRG4 opposing the Serum Starvation phenotype.
*   **`amd_risk_gene_validation.csv`**: The numerical data for the scatter plot.

### 8. `sc-analysis/` (Single-Cell Investigation)
*   **`cell_type_enrichment_heatmap.png`**: Heatmap showing how PRG4 rescue affects RPE vs Fibroblast markers.
*   **`prg4_receptor_bulk.csv`**: Expression changes of CD44, TLR2, PRG4 in the bulk data.
