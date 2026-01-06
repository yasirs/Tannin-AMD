# ATAC-seq Integration Results

This directory contains results from integrating the **GSE99287 RPE ATAC-seq Atlas** with internal PRG4/H2O2 bulk RNA-seq data.

## Key Findings
1.  **Global Accessibility Loss:** All differentially accessible regions (DARs) in the GSE99287 dataset represent a loss of accessibility in AMD patients.
2.  **Rescue Candidates:** We identified **104 genes** that both lose chromatin accessibility in AMD and are rescued (upregulated) by PRG4 treatment in our H2O2 stress model.
3.  **Top Hits:**
    *   **PURG**: Strongly rescued by PRG4.
    *   **TLR3**: Known involvement in RPE inflammation and AMD.
    *   **FAM161A**: Linked to retinitis pigmentosa and ciliary function.

## Files
*   **`atac_integration_summary.md`**: Statistical summary and top 20 candidates.
*   **`top_prg4_rescue_dar_candidates.csv`**: Full list of the 104 candidate genes.
*   **`atac_vs_h2o2_scatter.pdf`**: Correlation between AMD chromatin changes and H2O2 transcriptional response.
*   **`prg4_rescue_on_dars_boxplot.pdf`**: Comparison of PRG4 rescue effect on DAR-linked vs non-DAR genes.

## Analysis Script
`code/bridge-analysis/analyze_atac_integration.py`
