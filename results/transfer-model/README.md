# Transfer Model Results (Task 6)

**Objective:** To quantify the *scaling relationship* between perturbation effects in K562 vs RPE1. Does a hit in K562 imply a stronger or weaker effect in RPE?

## Data
*   **Training Set:** The Top 500 Concordant Genes identified in Task 4.
*   **Features:** Transcriptome-wide Z-scores (~3.6 million data points).

## Methods
*   **Linear Regression:** We fit the model $LFC_{RPE} = \beta \cdot LFC_{K562} + \alpha$.

## Key Findings
*   **Scaling Factor ($\beta$):** **0.705**.
*   **Interpretation:** A gene knockdown in RPE1 typically produces **~70%** of the effect magnitude observed in K562.
*   **Application:** When predicting RPE efficacy from K562 screens, we effectively "discount" the signal strength by ~30% to avoid overestimation.

## Files
*   **`transfer_scatter_density.pdf`**: Plot visualizing the global relationship between the cell types.
*   **`model_stats.txt`**: Regression coefficients.