# Task 4: Cell-Type Concordance Analysis

## Objective
Quantify how similar or different perturbation effects are between RPE1 and K562 cell lines to determine if cross-cell-type transfer is biologically plausible.

## Motivation
Before investing in a transfer model (Task 6), we need to answer:
1. **Are perturbation effects conserved across cell types?**
   - Do the same gene knockdowns produce similar transcriptional changes?
2. **Which genes/pathways show conserved vs cell-type-specific responses?**
3. **Is the concordance strong enough to justify using K562 data to inform RPE biology?**

## Data Sources

### Perturb-seq Data
- **RPE1 Essential**: `data/external/perturbseq/rpe1_normalized_bulk_01.h5ad` (~2,679 KDs)
- **K562 Essential**: `data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad` (~2,285 KDs)

### Filtered Gene Set
- **Task 3 output**: `results/baseline-expression/commonly_expressed_genes.csv`
  - Use ONLY genes expressed in both cell types for valid comparisons

### Overlapping Perturbations
- Extract the intersection of KDs present in both RPE1 and K562 Essential datasets
- Expected overlap: ~2,000-2,500 genes

## Analysis Tasks

### 1. Identify Overlapping Knockdowns
- Match perturbations by **target gene** (not guide ID, since guides differ)
- For genes with multiple guides: average the expression profiles, or analyze per-guide
- Create a list of perturbations to compare

### 2. Per-Perturbation Concordance
For each overlapping knockdown:
- Extract the expression profile (Z-normalized pseudo-bulk) from RPE1
- Extract the expression profile from K562
- Filter to **commonly expressed genes** (from Task 3)
- Calculate **Spearman correlation** between the two profiles

Result: A distribution of correlations (one per KD) showing how conserved each perturbation is.

### 3. Overall Concordance Metrics
- **Median concordance**: Median of all per-KD correlations
- **Distribution analysis**: 
  - What % of KDs have ρ > 0.5 (strong conservation)?
  - What % have ρ < 0.2 (cell-type-specific)?
- **Gene-wise concordance**: For each measured gene, how correlated is its response across all KDs?

### 4. Functional Stratification
Partition KDs by concordance level and perform enrichment:
- **Highly conserved KDs** (e.g., ρ > 0.6): Which biological processes?
  - Hypothesis: Core cellular functions (ribosome, proteasome, DNA repair)
- **Cell-type-specific KDs** (e.g., ρ < 0.2): Which processes?
  - Hypothesis: Lineage-specific pathways, metabolism, signaling

Perform GO/KEGG enrichment on each concordance bin.

### 5. Gene-Level Response Patterns
For each gene (in commonly expressed set):
- Calculate correlation of its **response vector** across all KDs in RPE1 vs K562
- Identify genes whose responses are:
  - **Universally conserved** (high correlation)
  - **Context-dependent** (low correlation)

Example: TF targets might be conserved, while metabolic genes might differ.

### 6. Scatter Plot Analysis
Create scatter plots for select KDs:
- X-axis: K562 expression change for gene i
- Y-axis: RPE1 expression change for gene i
- Each point = one measured gene
- Color by pathway or gene category

Визуальный QC: Do points cluster along diagonal (conservation) or scatter (divergence)?

## Expected Outputs

### Visualizations
1. **Concordance distribution histogram**:
   - X-axis = Spearman ρ (RPE1 vs K562 for each KD)
   - Y-axis = Frequency
   
2. **Concordance heatmap**:
   - Rows = KDs, Columns = Concordance bins
   - Color = correlation strength
   
3. **Enrichment bar charts**:
   - Conserved vs divergent KD pathway enrichment
   
4. **Scatter plots**: Representative KDs (high, medium, low concordance)

5. **Gene-wise conservation heatmap**: 
   - Genes vs concordance metric

### Files to Generate
- `results/concordance-analysis/overlapping_knockdowns.csv`
- `results/concordance-analysis/per_KD_concordance.csv`
- `results/concordance-analysis/conserved_KDs.csv` (high concordance)
- `results/concordance-analysis/divergent_KDs.csv` (low concordance)
- `results/concordance-analysis/gene_response_conservation.csv`
- `results/concordance-analysis/concordance_enrichment.csv`
- `results/concordance-analysis/concordance_summary.md`

### Key Metrics
- **Median concordance**: Overall ρ across all KDs
- **% conserved**: Fraction with ρ > 0.5
- **Top conserved pathways**: GO terms for high-concordance KDs
- **Top divergent pathways**: GO terms for low-concordance KDs

## Success Criteria
- **If median concordance ρ > 0.4**: Transfer modeling is plausible → Proceed to Task 6
- **If concordance is pathway-specific**: Build pathway-stratified transfer models
- **If concordance is low (ρ < 0.3)**: K562 may not be suitable for RPE predictions → Focus on RPE1-only analyses

## Dependencies
- **Task 3 (Baseline Expression Filter)**: Must complete to get commonly expressed gene set

## Next Steps
- **If concordance is strong**: Proceed to Task 6 (Transfer Model Development)
- **If concordance is weak**: Re-evaluate strategy; may need to find other RPE Perturb-seq data or focus on literature-based inference
- Results inform **Task 5** (Known Biology Validation) to verify findings with ground truth
