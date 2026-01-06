# Task 1: Signature Robustness Check

## Objective
Validate the stability and biological relevance of the bulk RNA-seq signatures used for bridge analysis. Determine if the current threshold (p < 0.05 nominal) is appropriate or if more stringent criteria should be used.

## Motivation
The current bridge analysis uses **nominal p < 0.05**, resulting in ~5,000-6,000 DEGs per contrast. This may be inflated and could introduce noise. We need to assess:
1. How stable are the top correlations across different significance thresholds?
2. Are we capturing real biology or statistical artifacts?
3. What's the optimal threshold for maximizing signal while minimizing noise?

## Data Sources

### Input Data
- **Bulk RNA-seq DEG Results**: `data/RPE_cells/code/RPE_gene pvals.xlsx`
  - Contains: `log2FoldChange`, `pvalue`, `padj` (FDR-adjusted p-value)
  - Three contrasts:
    - `H2O2_vs_CTRL`
    - `PRG4_vs_CTRL`
    - `H2O2PRG4_vs_H2O2`

- **Existing Bridge Results**: `results/bridge-results/`
  - `H2O2_Stress_bridge_correlations.csv`
  - `PRG4_Baseline_bridge_correlations.csv`
  - `PRG4_Rescue_bridge_correlations.csv`

- **Perturb-seq Data**: `data/external/perturbseq/rpe1_normalized_bulk_01.h5ad`

## Analysis Tasks

### 1. Threshold Sensitivity Analysis
Compare signatures at multiple thresholds:
- **p < 0.05** (nominal, current)
- **p < 0.01** (nominal, more stringent)
- **FDR < 0.1** (standard for publications)
- **FDR < 0.05** (very stringent)
- **Top N genes by |log2FC|** (e.g., top 1000, 2000, 3000)

For each threshold, record:
- Number of DEGs
- Number of genes overlapping with Perturb-seq
- Summary statistics of log2FC distribution

### 2. Re-run Bridge Correlations
For each threshold variant:
- Calculate Spearman correlations with RPE1 Perturb-seq
- Rank knockdowns by correlation
- Extract top 20-50 hits

### 3. Stability Assessment
Quantify:
- **Correlation stability**: Do the same KDs appear in top hits across thresholds?
  - Jaccard similarity between top-N hit lists
  - Rank correlation of KD rankings across thresholds
  
- **Effect size relationship**: Plot correlation coefficient vs signature size
  - Does using more genes dilute the signal or enrich it?

### 4. Biological Validation
For each threshold:
- Perform GO/KEGG enrichment on the signature genes
- Check if top-ranked KDs from bridge analysis hit known pathways
  - H2O2 stress → oxidative stress, NRF2 targets, DNA damage
  - PRG4 → inflammation, ECM remodeling, cytoprotection

## Expected Outputs

### Visualizations
1. **Threshold comparison table**: 
   - Rows = thresholds, Columns = [N_DEGs, N_overlap, top_KD, top_correlation]

2. **Stability heatmap**:
   - Jaccard similarity between top-N hit lists across thresholds

3. **Enrichment comparison**:
   - Bar charts showing top pathways for each threshold's signature

4. **Rank stability plot**:
   - Scatter plots showing KD rank consistency across threshold pairs

### Files to Generate
- `results/robustness-analysis/threshold_comparison.csv`
- `results/robustness-analysis/stability_metrics.csv`
- `results/robustness-analysis/signature_enrichment_summary.csv`
- `results/robustness-analysis/robustness_summary.md` (interpretive summary)

## Success Criteria
- Identify the **optimal threshold** that balances:
  - Sufficient gene coverage
  - Strong, reproducible correlations
  - Biologically coherent enrichment
  
- Document which threshold to use for downstream analyses (Tasks 2-9)

## Dependencies
**None** - This is a foundational task.

## Next Steps
Results inform which signature definition to use in:
- Task 2: Coverage Analysis
- Task 4: Cell-Type Concordance
- Task 7: Virtual PRG4 Screen
