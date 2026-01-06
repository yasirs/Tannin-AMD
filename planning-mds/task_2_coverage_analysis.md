# Task 2: Coverage Analysis

## Objective
Quantify the limitation of using **essential-gene-only** Perturb-seq datasets by measuring what fraction of our bulk RNA-seq biology is actually addressable in each Perturb-seq screen.

## Motivation
The RPE1 Perturb-seq dataset contains only **essential genes** (~2,679 KDs), which may miss critical AMD/PRG4 biology. The K562 GWPS has broader coverage (~11,258 KDs). 

**Key Question**: What percentage of our differentially expressed genes are even *perturbed* in the available datasets?

## Data Sources

### Input Data
- **Bulk RNA-seq DEG Results**: `data/RPE_cells/code/RPE_gene pvals.xlsx`
  - Use the **optimal threshold** determined in Task 1
  - Three contrasts: H2O2_vs_CTRL, PRG4_vs_CTRL, H2O2PRG4_vs_H2O2

- **Perturb-seq Datasets**:
  - RPE1 Essential: `data/external/perturbseq/rpe1_normalized_bulk_01.h5ad`
  - K562 Essential: `data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad`
  - K562 GWPS: `data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad`

### Gene ID Mapping
- Use **Ensembl Gene IDs** for matching
- Both bulk and Perturb-seq data include Ensembl IDs

## Analysis Tasks

### 1. Per-Contrast Coverage
For each of the 3 bulk RNA-seq contrasts:
- Extract DEG list (using Task 1's optimal threshold)
- Count how many DEGs are:
  - **Measured** in each Perturb-seq dataset (in `.var`)
  - **Perturbed** in each dataset (in `.obs`, i.e., actually knocked down)

Create a summary table:

| Contrast | N_DEGs | Measured (RPE1) | KD'd (RPE1) | Measured (K562_Ess) | KD'd (K562_Ess) | Measured (K562_GWPS) | KD'd (K562_GWPS) |
|:---------|:-------|:----------------|:------------|:--------------------|:----------------|:---------------------|:-----------------|
| H2O2 Stress | ... | ... | ... | ... | ... | ... | ... |
| PRG4 Baseline | ... | ... | ... | ... | ... | ... | ... |
| PRG4 Rescue | ... | ... | ... | ... | ... | ... | ... |

### 2. Venn Diagrams
For each contrast, create 3-way Venn diagrams showing overlap between:
- Bulk DEGs
- RPE1 KDs
- K562 GWPS KDs

### 3. Biological Characterization of "Missing" Genes
For genes that are:
- **DE in bulk** but **NOT perturbed in RPE1**
- **DE in bulk** but **NOT perturbed in K562 GWPS**

Characterize each "missing set":
- Perform GO/KEGG enrichment
- Identify key pathways not addressable by essential-only screens
- Check if missing genes are known AMD risk factors (literature search optional)

### 4. Gain from K562 GWPS
Quantify the additional coverage provided by K562 GWPS:
- How many more DEGs are perturbed in K562 GWPS vs RPE1?
- Are the additionally covered genes enriched for specific biology?

## Expected Outputs

### Visualizations
1. **Coverage bar charts**: 
   - X-axis = Dataset (RPE1, K562_Ess, K562_GWPS)
   - Y-axis = % of DEGs covered
   - Grouped by contrast

2. **Venn diagrams**: One per contrast

3. **Enrichment bar charts**: Top pathways in "missing gene" sets

### Files to Generate
- `results/coverage-analysis/coverage_summary.csv`
- `results/coverage-analysis/missing_genes_H2O2.csv`
- `results/coverage-analysis/missing_genes_PRG4_baseline.csv`
- `results/coverage-analysis/missing_genes_PRG4_rescue.csv`
- `results/coverage-analysis/missing_genes_enrichment.csv`
- `results/coverage-analysis/coverage_analysis_summary.md` (interpretive report)

### Key Metrics
- **Coverage rate**: % of DEGs with KDs in each dataset
- **Incremental gain**: Additional coverage from K562 GWPS over RPE1
- **Biological gap**: Pathways missing from essential-only screens

## Success Criteria
- Clear quantification of the essential-gene limitation
- Evidence-based justification for using K562 GWPS (if coverage is substantially higher)
- Identification of specific biological processes not addressable with current data

## Dependencies
- **Task 1 (Signature Robustness)**: Must complete first to determine which DEG threshold to use

## Next Steps
If K562 GWPS provides substantial additional coverage:
- Proceed to **Task 4** (Cell-Type Concordance) to assess if K562 data is transferable
- Justify investment in **Task 6** (Transfer Model) and **Task 7** (Virtual Screen)

If coverage is similar across datasets:
- May not need transfer modeling; can focus on RPE1-only analyses
