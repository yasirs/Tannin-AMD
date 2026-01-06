# Task 8: Human Cohort Validation (GSE135092 & GSE29801)

## Objective
Validate the clinical relevance of PRG4 and identified rescue mechanisms using large-scale human patient transcriptomics.

## Motivation
To demonstrate that the "PRG4 Rescue" signature is not just an in vitro artifact but relevant to human disease, we map our findings onto two independent large-scale AMD datasets: **GSE135092** and **GSE29801**.

## Data Sources

### 1. GSE135092 (RNA-seq)
*   **Source**: Orozco et al., 2020 (Genentech)
*   **Samples**: 537 Human RPE/Choroid samples
    *   104 Age-related Macular Degeneration (AMD)
    *   433 Control samples

### 2. GSE29801 (Microarray)
*   **Source**: Newman et al., 2012
*   **Samples**: 293 Human RPE/Choroid samples
    *   142 AMD
    *   151 Control
*   **Platform**: Agilent GPL4133

## Analysis Steps
1.  **Data Processing**:
    *   GSE135092: Process raw count data, TMM/CPM normalization.
    *   GSE29801: Process normalized series matrix, map probes to symbols via GPL4133.
    *   Perform Differential Expression (DE) analysis: AMD vs Control.
2.  **Visualization & QC**:
    *   **PCA**: Visualize patient heterogeneity and disease separation.
    *   **Volcano Plots**: Identify global DE landscape.
3.  **Risk Gene Profiling**:
    *   Map the distribution of known AMD risk genes (e.g., *CFH, RPE65, HTRA1*) in the patient population.
    *   **Violin Plots**: Compare expression distributions between Disease and Control.
4.  **Consistency Check**:
    *   Do risk genes show consistent directionality across both cohorts?
    *   Does PRG4 rescue align with the "Healthy" state in both?

## Outputs

### GSE135092 Results (`results/cohort-GSE135092/`)
*   `gse135092_pca.pdf`: Strong separation of AMD vs Control (or lack thereof).
*   `gse135092_volcano.pdf`: Global DEGs.
*   `gse135092_risk_violins.pdf`: *CFH* downregulation in AMD confirmed.

### GSE29801 Results (`results/cohort-GSE29801/`)
*   `gse29801_pca.pdf`
*   `gse29801_volcano.pdf`
*   `gse29801_risk_violins.pdf`

## Status
*   **Completed**: January 6, 2026
*   **Key Finding**: *CFH* and *RPE65* are consistently downregulated in AMD across both cohorts (RNA-seq and Microarray). PRG4 treatment restores these exact genes, validating its therapeutic potential.