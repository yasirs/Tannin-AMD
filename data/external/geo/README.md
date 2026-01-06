# GEO Datasets

Organized raw and processed data from NCBI Gene Expression Omnibus.

## Datasets

### GSE135092 (AMD RNA-seq)
*   **Source:** [Lydia Parker et al.]
*   **Contents:**
    *   `raw_data/`: Extracted raw count files for each sample.
    *   `*_expression_matrix.csv.gz`: Combined count matrix (Samples x Genes).
    *   `*_amd_signature.csv`: Differential expression results (AMD vs Control).

### GSE29801 (AMD Microarray)
*   **Source:** [Unpublished?]
*   **Contents:**
    *   `raw_data/`: Extracted raw probe data.
    *   `*_expression_matrix.csv.gz`: Probe-level expression matrix.
    *   `*_amd_signature.csv`: Differential expression results (AMD vs Control).

### GSE99287 (AMD RPE ATAC-seq)
*   **Source:** Wang et al., 2018 (PMID: 29634994)
*   **Description:** Chromatin accessibility landscape of human RPE in AMD vs Control.
*   **Contents:**
    *   `GSE99287_RPE_ATACSeq_peak_counts.txt.gz`: Peak count matrix.
    *   `GSE99287_peak_counts_annotation.txt.gz`: Peak annotations (nearest gene, genomic location).

### GSE129964 (Serum Starvation)
*   **Description:** Timecourse of RPE differentiation/atrophy.
