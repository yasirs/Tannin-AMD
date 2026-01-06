# External Data

This directory houses third-party datasets used for validation and expansion of the core analysis.

## Datasets

### 1. Transcriptomics (GEO)
Located in `geo/`.

*   **GSE135092**: The primary human validation cohort. RNA-seq of Macula/Periphery/Retina from AMD and Control donors.
*   **GSE29801**: Independent validation cohort. Microarray data.
*   **GSE99287**: Epigenomic landscape of AMD RPE (ATAC-seq). Used to map PRG4 effects to chromatin regulatory elements.

### 2. Functional Genomics (Perturb-seq)
Located in `perturbseq/`.

*   **K562 GWPS**: Genome-wide CRISPRi screen (Replogle et al., 2022). Used for the Virtual Screen.
*   **RPE1 Essential**: Essential gene screen in RPE1 cells. Used for cell-type concordance checks.

## Scripts
*   Download scripts are located in `code/external-validation/`.
