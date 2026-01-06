# Cell-Type Concordance Results (Task 4)

**Objective:** To determine if K562 (Leukemia) cells are a valid proxy for RPE1 (Retinal) cells by quantitatively comparing the transcriptional effects of knocking down the same genes in both lines.

## Data Sources
*   **RPE1 Essential** (Replogle et al. 2022).
*   **K562 Essential** (Replogle et al. 2022).

## Methods
1.  **Alignment:** We extracted the perturbation profiles (Z-normalized expression changes) for the **~2,000 shared essential genes** targeted in both libraries.
2.  **Correlation:** For each gene knockdown, we calculated the **Pearson correlation** between its RPE1 profile and its K562 profile.

## Key Findings
*   **Global Correlation:** The mean concordance is **r ~ 0.19**, indicating significant cell-type specificity.
*   **High-Concordance Core:** We identified a subset of **~500 genes** (e.g., Proteasome, Ribosome subunits) with r > 0.3. These "Housekeeping" functions are highly conserved and safe to model in K562.
*   **Implication:** K562 is a valid model for *core* cellular stress responses (like the NRF2 pathway) but less suitable for lineage-specific factors.

## Files
*   **`gene_concordance.csv`**: The lookup table of transferability scores.
*   **`concordance_distribution.pdf`**: Histogram of correlation coefficients.