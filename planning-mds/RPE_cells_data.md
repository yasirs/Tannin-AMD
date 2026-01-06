# RPE Cells Transcriptomics Data

This directory contains transcriptomics data and analysis for RPE (Retinal Pigment Epithelium) cells under various conditions, including Hydrogen Peroxide (H2O2) stress and Lubricin (PRG4) treatment.

## Experimental Design

The study analyzed **12 samples** categorized into four conditions (3 biological replicates each):

| Condition | Samples | Description |
| :--- | :--- | :--- |
| **CTRL** | TS1, TS2, TS3 | Control (Untreated) |
| **H2O2** | TS5, TS6, TS7 | Hydrogen Peroxide Stress |
| **PRG4** | TS9, TS11, TS12 | Lubricin (PRG4) Treatment |
| **H2O2PRG4** | TS13, TS15, TS16 | Co-treatment (H2O2 + PRG4) |

## Data Structure

### Raw Data
*   `data/count/`: Contains `.counts` files (featureCounts output) for each sample.

### Processed Results (in `code/`)
The primary processed files are located in the `code/` subdirectory:

1.  **`RPE_gene pvals.xlsx`**: Differential expression results calculated using DESeq2.
    *   Contains Log2 Fold Change and p-values (nominal and adjusted) for three key contrasts:
        *   `H2O2_vs_CTRL`
        *   `PRG4_vs_CTRL`
        *   `H2O2PRG4_vs_H2O2`
    *   Includes Ensembl Gene IDs and HGNC Symbols.

2.  **`RPE_TPMS.xlsx`**: Transcripts Per Million (TPM) normalized expression values.
    *   Columns correspond to individual samples (`TS1`, `TS2`, etc.).
    *   Used for heatmaps and comparative expression analysis.

3.  **`RPE_RPMs.xlsx`**: Reads Per Million (RPM) normalized values.

4.  **`TanninReadCounts.xlsx`**: Raw read counts extracted from the featureCounts files.

5.  **`RPE_gsa.xlsx`**: Gene Set Analysis (GSA) results for KEGG and GO terms across the specified contrasts.

## Analysis Logic
The analysis was performed using the `code/analysis.R` script, which:
*   Imports raw counts from `../data/count/`.
*   Performs PCA and variance-stabilizing transformation (VST).
*   Conducts differential expression analysis using **DESeq2**.
*   Calculates TPM/RPM normalization.
*   Generates heatmaps (e.g., `RPE_ARE_heatmap.tiff`) for specific gene sets related to oxidative stress response (e.g., MT2A, HMOX1, G6PD).
*   Runs GSA for functional enrichment.
