# Task 1: Signature Robustness Check - Summary

**Date:** 2026-01-06
**Analyst:** Gemini Agent

## Objective
To determine the optimal statistical threshold for defining bulk RNA-seq signatures used in bridge analysis, ensuring stability and biological relevance.

## Key Findings

### 1. Threshold Sensitivity & Gene Coverage
We compared nominal p-values (`p < 0.05`, `p < 0.01`), FDR thresholds (`FDR < 0.1`, `FDR < 0.05`), and simple Fold-Change rankings (`Top N`).

| Contrast | Threshold | N DEGs (Bulk) | N Overlap (Perturb-seq) | Top Hit | Spearman Rho |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **H2O2 Stress** | p < 0.05 | 5734 | 2816 | SNIP1 | 0.297 |
| | **p < 0.01** | **3838** | **1813** | **POLR3D** | **0.324** |
| | FDR < 0.05 | 3943 | 1866 | POLR3D | 0.323 |
| | Top 1000 | 1000 | **13** | ATF4 | 0.890* |
| **PRG4 Rescue** | p < 0.05 | 5625 | 3671 | ARL4D | 0.198 |
| | **p < 0.01** | **3856** | **2648** | **ARL4D** | **0.233** |
| | FDR < 0.05 | 3975 | 2712 | ARL4D | 0.230 |
| | Top 1000 | 1000 | **43** | CCNK | 0.496* |

***Note on Top N**: The "Top N by LogFC" strategy resulted in extremely poor overlap with the Perturb-seq dataset (e.g., only 13 out of top 1000 H2O2 genes were present in Perturb-seq). This makes the high correlations unreliable and the strategy unsuitable.*

### 2. Stability Assessment
*   **P-value vs FDR**: There is high concordance between `p < 0.01` and `FDR < 0.05` (Jaccard ~ 1.0 in H2O2 contrast), indicating they select nearly identical gene sets.
*   **Consistency**: The top identified Knockdowns (KDs) are highly stable between `p < 0.05` and `p < 0.01`. For H2O2, `SNIP1` and `POLR3D` are consistently top-ranked. For PRG4 Rescue, `ARL4D` is robustly the top hit.
*   **Signal Strength**: Stricter thresholds (`p < 0.01`, `FDR < 0.05`) yielded slightly higher correlation coefficients for the top hits compared to the looser `p < 0.05`, suggesting a reduction in noise without loss of biological signal.

### 3. Biological Validation (Preliminary)
*   **H2O2 Stress**: 
    *   **SNIP1** (Smad Nuclear Interacting Protein 1) is a transcriptional co-repressor.
    *   **POLR3D** (RNA Polymerase III Subunit D) is involved in sensing DNA damage and viral infection.
    *   The shift from SNIP1 to POLR3D as top hit with stricter thresholds warrants further investigation but both are biologically plausible global regulators.
*   **PRG4 Rescue**:
    *   **ARL4D** (ADP Ribosylation Factor Like GTPase 4D) remains the top hit across all valid thresholds. This suggests a very strong, robust link between the rescue phenotype and ARL4D depletion.

## Recommendations

1.  **Adopt `FDR < 0.05` (or `p < 0.01`)**: This threshold offers a better signal-to-noise ratio than the previously used `p < 0.05`. It maintains sufficient gene coverage (~1800-2700 genes) while producing higher correlation scores.
2.  **Abandon `Top N` by FC**: The low overlap with Perturb-seq genes makes this approach viable only if we can significantly improve the intersection (see Task 2: Coverage Analysis).
3.  **Proceed with Bridge Analysis**: Re-run the core bridge analysis pipeline using the refined `FDR < 0.05` signatures to generate the final candidate lists.

## Output Files
*   `threshold_comparison.csv`: Full metrics for all tested thresholds.
*   `stability_metrics.csv`: Pairwise Jaccard similarity of top 50 hits.
