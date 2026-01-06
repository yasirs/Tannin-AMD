# Task 3: Baseline Expression Filtering

## Objective
Establish which genes are expressed in both RPE1 and K562 cell lines to ensure valid cross-cell-type comparisons. Genes not expressed in a cell type won't show meaningful perturbation effects.

## Motivation
K562 (myeloid leukemia) and RPE1 (retinal epithelial) are very different cell types:
- Different baseline transcriptomes
- Different lineage-specific programs
- Different biological contexts (cancer vs immortalized normal)

**Key Insight**: A gene KD can only reveal functional effects if that gene is actually expressed. Comparing perturbation profiles for a gene expressed in K562 but not RPE1 (or vice versa) is meaningless.

## Data Sources

### Baseline Expression Data
We need to infer baseline expression from the Perturb-seq data itself:

- **RPE1**: `data/external/perturbseq/rpe1_raw_bulk_01.h5ad`
  - Use **raw pseudo-bulk** data (not Z-normalized)
  - Expression across all perturbations can proxy for baseline

- **K562 Essential**: `data/external/perturbseq/K562_essential_raw_bulk_01.h5ad`
- **K562 GWPS**: `data/external/perturbseq/K562_gwps_raw_bulk_01.h5ad`

### Approach
Since we don't have pure unperturbed RNA-seq for these cell lines, we can:
1. Use **non-targeting controls** from each dataset as baseline
2. Use the **median/mean expression** across all perturbations (most KDs don't affect most genes)
3. Define "expressed" as: mean UMI/TPM > threshold (e.g., >0.1 or >1.0)

## Analysis Tasks

### 1. Extract Baseline Expression
For each dataset (RPE1, K562_Ess, K562_GWPS):
- Identify **non-targeting control** perturbations (if available in `.obs`)
- If controls exist: Average expression across controls
- If not: Use median expression across all perturbations

### 2. Define Expression Thresholds
Test multiple thresholds to define "expressed":
- **Loose**: Mean expression > 0.1 UMI/cell
- **Moderate**: Mean expression > 1.0 UMI/cell
- **Stringent**: Mean expression in top 50% of genes, or detected in >50% of perturbations

Document the number of genes expressed in each cell type at each threshold.

### 3. Create Filtered Gene Sets
For cross-cell-type analyses, create gene sets that are:
- **Commonly expressed**: Expressed in BOTH RPE1 and K562 (at a chosen threshold)
- **RPE1-specific**: Expressed in RPE1 only
- **K562-specific**: Expressed in K562 only

### 4. Characterize Expression Differences
- Compare baseline expression distributions between cell types
- Identify lineage-specific genes (e.g., visual cycle genes in RPE1, heme biosynthesis in K562)
- Perform GO enrichment on cell-type-specific expression sets

### 5. Validation with Known Markers
Check if expected markers are correctly classified:
- **RPE-specific**: BEST1, RPE65, RLBP1, TTR
- **K562-specific**: HBE1, HBG1/2 (globins), CD34
- **Housekeeping**: ACTB, GAPDH, TBP (should be in both)

## Expected Outputs

### Gene Lists
- `results/baseline-expression/expressed_genes_RPE1.csv`
- `results/baseline-expression/expressed_genes_K562_essential.csv`
- `results/baseline-expression/expressed_genes_K562_gwps.csv`
- `results/baseline-expression/commonly_expressed_genes.csv`
- `results/baseline-expression/cell_type_specific_genes.csv`

### Visualizations
1. **Expression distribution plots**: Histogram or density plots of baseline expression in each cell type
2. **Venn diagram**: Overlap of expressed genes between RPE1 and K562
3. **Scatter plot**: RPE1 baseline vs K562 baseline expression for commonly measured genes
4. **Heatmap**: Expression of known markers across cell types

### Summary Statistics
- Number of genes expressed in each cell type
- Number commonly expressed
- Correlation of baseline expression for common genes

### Files to Generate
- `results/baseline-expression/baseline_expression_summary.csv`
- `results/baseline-expression/cell_type_enrichment.csv`
- `results/baseline-expression/baseline_expression_report.md`

## Success Criteria
- Clear definition of which genes can be validly compared across cell types
- Validation that expected lineage markers behave correctly
- Filtered gene sets for use in downstream concordance and transfer analyses

## Dependencies
**None** - Can run in parallel with Tasks 1 and 2.

## Next Steps
Results directly feed into:
- **Task 4 (Cell-Type Concordance)**: Use commonly expressed genes for valid comparisons
- **Task 6 (Transfer Model)**: Train only on genes expressed in both contexts
- **Task 7 (Virtual Screen)**: Filter K562 GWPS results to RPE-relevant genes
