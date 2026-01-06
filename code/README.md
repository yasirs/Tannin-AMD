# Code Directory - Tannin-AMD Project

This directory contains all analysis scripts used to investigate the role of PRG4 in AMD.

## 1. Bridge Analysis (Core)
Scripts connecting internal RPE bulk RNA-seq data with external Perturb-seq and GEO datasets.

*   **`robustness_check.py`**: (Task 1) Tests various p-value thresholds to find the most stable signature for finding hits in Perturb-seq. Recommended: FDR < 0.05.
*   **`coverage_check.py`**: (Task 2) Quantifies how many of our DEGs are present in RPE1 vs K562 Perturb-seq datasets. Found K562 GWPS is essential.
*   **`extract_baseline.py`**: (Task 3) Calculates baseline expression from raw Perturb-seq files to define the "addressable transcriptome".
*   **`analyze_concordance.py`**: (Task 4) Correlates knockdown profiles between RPE1 and K562 to validate cell-type transferability.
*   **`transfer_model_analysis.py`**: (Task 6) Builds a linear regression model ($LFC_{RPE} \sim LFC_{K562}$) to quantify the scaling relationship between cell types.
*   **`virtual_screen.py`**: (Task 7) The main driver. Correlates PRG4 Rescue signature with ~11,000 K562 knockdowns to find mimetics (e.g., KEAP1) and antagonists (e.g., ATR).
*   **`validate_hits.py`**: (Task 5) Runs pathway enrichment on the top hits from the virtual screen.
*   **`run_gsea_custom.py`**: Implements a custom Pre-ranked GSEA algorithm to compare our signatures against external datasets (Serum Starvation).
*   **`analyze_leading_edge.py`**: Analyzes the core genes driving the GSEA enrichment and checks if mimetics target them.

## 2. External Data Processing
Scripts to fetch and process public GEO datasets.

*   **`fetch_external_data.py`**: Attempts to download and parse GEO Series Matrix files. (Note: Metadata extraction works, but expression data often requires raw processing).
*   **`process_gse135092_raw.py`**: Processes 500+ raw TSV files from the GSE135092 (AMD RNA-seq) dataset to build a count matrix and run differential expression.
*   **`analyze_gse129964.py`**: Analyzes the GSE129964 (Serum Starvation) count matrix to define a "Stress" signature.
*   **`validate_amd_risk.py`**: Checks if known AMD risk genes (CFH, RPE65) are rescued in our internal data.

## 3. Single-Cell Analysis
*   **`cell_type_projection.py`**: Projects bulk signatures onto cell-type specific markers (RPE, Endothelial, Fibroblast) to infer tissue-level effects.

## 4. Visualization
R scripts for high-quality, Tufte-style plotting.

*   **`visualize_robustness.R`**: Plots for Task 1 (Threshold vs Correlation).
*   **`visualize_coverage.R`**: Plots for Task 2 (Bar charts of gene coverage).
*   **`visualize_baseline.R`**: Plots for Task 3 (Expression distributions).
*   **`visualize_concordance.R`**: Plots for Task 4 (Histogram of correlations).
*   **`plot_gsea_standard.R`**: Generates standard enrichment plots (green curve/black rug) for the GSEA results.
*   **`make_final_plots.R`**: Generates the complex external validation plots (Serum Scatter, AMD Patient Densities).