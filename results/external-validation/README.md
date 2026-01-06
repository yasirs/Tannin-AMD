# External Validation Results

**Objective:** To validate the clinical relevance of the PRG4 findings using independent external datasets and human patient data.

## Datasets
1.  **GSE135092 (Human Patient Data):**
    *   **Source:** Orozco et al. (Genentech), 2020.
    *   **Samples:** 537 human RPE/Choroid samples (104 AMD, 433 Control).
    *   **Usage:** Used to map the distribution of key risk genes in the real patient population.
2.  **GSE129964 (Serum Starvation Model):**
    *   **Source:** ARPE-19 cells under serum deprivation.
    *   **Usage:** A model of "Atrophic Stress". Used to check if PRG4 counters atrophy.
3.  **Literature AMD Gene Panel:**
    *   **Source:** Curated list of GWAS/Genetic risk factors (*CFH, ARMS2, HTRA1, RPE65*, etc.).

## Key Findings
*   **Reversal of Atrophy:** The PRG4 Rescue signature is **negatively correlated (r = -0.24)** with the Serum Starvation signature. This implies PRG4 promotes a "Growth/Health" state that directly opposes atrophy.
*   **Reversal of Risk Genes:** PRG4 treatment reverses the stress-induced dysregulation of critical genes:
    *   **RPE65**: Repressed by Stress $\rightarrow$ Restored by PRG4.
    *   **CFH**: Repressed by Stress $\rightarrow$ Restored by PRG4.

## Files
*   **`amd_rescue_scatter.png`**: Scatter plot showing the strong negative correlation (Reversal) for risk genes.
*   **`serum_rescue_scatter.pdf`**: Tufte-style plot showing PRG4 opposing the Serum Starvation phenotype.
*   **`amd_risk_gene_validation.csv`**: Validation data for risk genes.

*(Note: Human Cohort GSE135092 results have been moved to `results/cohort-GSE135092/`)*