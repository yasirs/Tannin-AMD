
# Analysis Preferences & Methodological Standards

This document outlines the standard methodologies to be used for statistical comparisons in this project. All future planning and analysis agents must adhere to these guidelines to ensure rigor and correctness.

## 1. Comparing Conditions with Continuous Scores

When comparing two conditions where continuous scores (e.g., expression levels, pathway scores, perturbation effects) are available, **do not use simple correlation plots** of gene sets. Instead, use the following approaches:

### A. Threshold-Based GSEA (One Condition defines a Set)
If one condition allows for defining a specific gene set (e.g., top K genes, significant DEGs), treat it as a **Gene Set** and compare it against the continuous **Score** of the second condition using **Gene Set Enrichment Analysis (GSEA)**.

*   **Method**: Run GSEA (e.g., via R package `fgsea`).
*   **Metrics to Record**:
    *   GSEA P-value
    *   Enrichment Score (ES)
    *   Normalized Enrichment Score (NES)
*   **Visualizations**:
    *   Standard GSEA enrichment plot.
    *   Heatmap of "Leading Edge" genes (if the number is appropriate, e.g., 4–50 genes).

### B. Rank-Based Comparison (Score vs. Score)
If comparing the distribution of scores directly without defining a set:

*   **Method**: Use a **Non-parametric Rank Comparison Test** (e.g., Wilcoxon Rank-Sum test / Mann-Whitney U test) to assess if the distributions or ranks differ significantly or are shifted.

## 2. Comparing Gene Sets (Discrete vs. Discrete)

When both conditions result in discrete gene sets (e.g., "Upregulated genes in Condition A" vs. "Upregulated genes in Condition B"):

### Over-Representation Analysis (ORA)
*   **Method**: Perform Over-Representation Analysis (e.g., Fisher’s Exact Test, Hypergeometric Test, or `fgsea::ora`).
*   **Metrics to Record**:
    *   P-value (Hypergeometric/Fisher)
    *   Odds Ratio
    *   Intersection Size
*   **Visualizations**:
    *   Venn Diagram.
    *   Heatmap of the overlapping genes.

## 3. Pathway Scoring & Sample-Level Analysis

When analyzing pathway activity, choice of method depends on the experimental context (Sample Size N, Homogeneity).

### A. Controlled Experiments (Small N, Clear Contrasts)
**Goal**: Determine if a pathway is induced or reversed by a treatment (e.g., Drug vs Control, H2O2 vs PBS).
*   **Avoid**: Computing mean Z-scores per sample and performing t-tests. This reduces complex data to a crude average and ignores gene rank/importance.
*   **Preferred Method**: **GSEA Normalized Enrichment Score (NES)**.
    *   Use the **NES** from standard GSEA (e.g., `fgsea`) as the primary "effect size" metric for the pathway change between conditions.
    *   Use the GSEA **FDR/P-value** for significance.

### B. Heterogeneous Population Cohorts (Large N)
**Goal**: Correlate pathway activity with clinical traits (Age, Disease Severity) across many individuals (e.g., TCGA, GEO clinical cohorts).
*   **Preferred Method**: **GSVA (Gene Set Variation Analysis)** or **ssGSEA**.
    *   **Why**: These methods calculate a **rank-based score per sample** that is robust to outliers and technical noise common in large datasets.
    *   **Usage**: Calculate GSVA scores for all samples $\rightarrow$ Use scores in standard regression models (`lm(Score ~ Age + Sex + Disease)`).

