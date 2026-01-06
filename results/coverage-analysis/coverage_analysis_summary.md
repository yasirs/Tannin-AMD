# Task 2: Coverage Analysis - Summary

**Date:** 2026-01-06
**Analyst:** Gemini Agent

## Objective
To quantify the coverage limitations of the "Essential-only" Perturb-seq datasets and determine if the Genome-Wide (GWPS) K562 dataset provides a necessary advantage for analyzing AMD-related signatures.

## Key Findings

### 1. Massive Coverage Gap in Essential Screens
The RPE1 Essential dataset, while cell-type matched, covers a very small fraction of the relevant biology found in our bulk RNA-seq signatures.
*   **H2O2 Stress:** RPE1 Essential covers only **11.5%** of the DEGs.
*   **PRG4 Rescue:** RPE1 Essential covers only **16.7%** of the DEGs.

### 2. K562 GWPS Provides Superior Coverage
Switching to the genome-wide K562 dataset results in a ~3x to 5x increase in addressable biology.
*   **H2O2 Stress:** Coverage jumps to **51.3%** (vs 11.5%).
*   **PRG4 Rescue:** Coverage jumps to **65.6%** (vs 16.7%).

| Contrast | RPE1 Essential | K562 GWPS | **Gain** |
| :--- | :--- | :--- | :--- |
| **H2O2 Stress** | 11.5% | 51.3% | **+39.8%** |
| **PRG4 Baseline** | 16.3% | 61.1% | **+44.8%** |
| **PRG4 Rescue** | 16.7% | 65.6% | **+48.9%** |

### 3. Missing Biology in RPE1 Essential
By analyzing the genes *missing* from the RPE1 Essential screen (but present in our signatures), we identified critical pathways that would be completely ignored if we stayed with RPE1 only:
*   **Neuroactive Ligand-Receptor Interaction**: Top missing pathway (p < 1e-18).
*   **Calcium Signaling**: Critical for RPE physiology (p < 1e-12).
*   **GPCR Signaling**: Major drug targets (p < 1e-10).
*   **ECM Interactions**: Key for AMD/Drusen pathology.

These are "non-essential" for cell survival in a dish but **essential for RPE function and disease**.

## Conclusion & Recommendation

**The RPE1 Essential dataset is insufficient for this project.** 
It misses >80% of the differentially expressed genes, including those in pathways central to AMD pathology (ECM, signaling).

**Recommendation:**
1.  **Adopt K562 GWPS**: We must use the K562 Genome-Wide dataset as our primary discovery engine to capture the full spectrum of AMD-related biology.
2.  **Validate Transferability**: Since K562 is a leukemia line, we must rigorously assess its concordance with RPE1 for the shared essential genes (Task 4) and potentially build a transfer model (Task 6) to correct for cell-type differences.

## Output Files
*   `coverage_summary.csv`: Detailed counts of measured/perturbed genes.
*   `coverage_comparison.pdf`: Visualization of the coverage gap.
*   `enrichment_missing_h2o2.pdf`: Pathways lost in the Essential-only screen.
*   `missing_genes_*.csv`: Lists of specific genes missing from each screen.
