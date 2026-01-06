# Task 5: Known Biology Validation

## Objective
Validate the cell-type concordance findings (Task 4) using **known biological ground truths** from the literature. Verify that well-established pathways behave as expected in both cell types.

## Motivation
Statistical concordance (Task 4) is necessary but not sufficient. We need to check:
1. **Do known oxidative stress regulators behave similarly in RPE1 and K562?**
2. **Are core cellular processes (DNA repair, ribosome, proteasome) conserved?**
3. **Do lineage-specific pathways correctly show divergence?**

This provides **biological validation** that the concordance metrics are meaningful.

## Data Sources

### Perturb-seq Data
- **RPE1 Essential**: `data/external/perturbseq/rpe1_normalized_bulk_01.h5ad`
- **K562 Essential**: `data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad`

### Task 4 Output
- `results/concordance-analysis/per_KD_concordance.csv`

### Known Biology Gene Sets
Curate or obtain from MSigDB/literature:

#### Core/Conserved Processes (expect HIGH concordance)
- **Ribosome biogenesis**: RPL/RPS genes
- **Proteasome**: PSMA/PSMB/PSMC genes
- **DNA repair**: BRCA1, RAD51, XRCC family, RFC2/4
- **Cell cycle checkpoints**: TP53, CDC family, PLK1

#### Oxidative Stress Response (critical for AMD biology)
- **NRF2 pathway**: NFE2L2, KEAP1, NQO1, HMOX1, GCLM
- **Antioxidant enzymes**: SOD1/2, GPX1-4, CAT, PRDX family
- **ROS sensors**: TXNRD1, GSR

#### Lineage-Specific (expect LOW concordance)
- **RPE-specific**: RPE65, BEST1, RLBP1, MERTK
- **Hematopoietic**: HBB, HBG1/2, CD34, GATA1

## Analysis Tasks

### 1. Extract Concordance for Known Gene Sets
For each curated gene set:
- Filter to genes that are KD'd in both RPE1 and K562
- Extract their concordance scores from Task 4 results
- Calculate summary statistics: mean, median, 25th/75th percentiles

### 2. Compare Concordance Across Categories
Statistical test (e.g., Mann-Whitney U or Kruskal-Wallis):
- H0: Conserved processes do NOT have higher concordance than lineage-specific
- Expected: Core processes >> lineage-specific

Create a **boxplot** comparing concordance distributions across gene set categories.

### 3. Oxidative Stress Deep Dive
Since this is central to AMD biology:
- For each NRF2/antioxidant gene KD:
  - Plot RPE1 vs K562 expression profiles (scatter plot)
  - Highlight canonical targets (NQO1, HMOX1, etc.)
  - Check if KD effects are conserved

- Calculate pathway-level concordance:
  - Average concordance of all oxidative stress-related KDs
  - Compare to genome-wide baseline

### 4. Individual Case Studies
Select 3-5 well-studied genes from literature and manually inspect:
- **Example 1**: KEAP1 KD → Should upregulate NRF2 targets in both cell types
- **Example 2**: TP53 KD → Should affect apoptosis/cell cycle genes
- **Example 3**: RPL5 (ribosomal) KD → Should trigger p53-independent nucleolar stress

For each:
- Visualize expression changes in RPE1 vs K562
- Check top affected genes and pathways
- Verify alignment with published findings

### 5. Literature Comparison
If published KD/KO RNA-seq exists for specific genes:
- Compare Perturb-seq profiles to bulk RNA-seq from literature
- Validate that Perturb-seq pseudo-bulk recapitulates expected biology

## Expected Outputs

### Visualizations
1. **Concordance boxplot by category**:
   - X-axis = Gene set category (Core, Oxidative Stress, Lineage-Specific)
   - Y-axis = Concordance (ρ)
   
2. **Oxidative stress heatmap**:
   - Rows = Oxidative stress KDs
   - Columns = Key target genes (NQO1, HMOX1, etc.)
   - Color = Expression change (Z-score)
   - Side-by-side for RPE1 and K562

3. **Case study scatter plots**: 
   - One per selected gene (e.g., KEAP1, TP53, RPL5)

4. **Pathway concordance bar chart**:
   - Top biological pathways ranked by average KD concordance

### Files to Generate
- `results/validation-analysis/gene_set_concordance.csv`
- `results/validation-analysis/oxidative_stress_validation.csv`
- `results/validation-analysis/case_study_details.csv`
- `results/validation-analysis/validation_summary.md`

### Key Metrics
- **Mean concordance for core processes**: Should be high (ρ > 0.5)
- **Mean concordance for lineage-specific**: Should be low (ρ < 0.3)
- **Oxidative stress concordance**: Critical benchmark for AMD relevance

## Success Criteria
- Core processes show significantly higher concordance than lineage-specific (p < 0.05)
- Oxidative stress pathway concordance ≥ 0.4 (moderate conservation)
- At least 2/3 case studies align with published findings

**If validation fails**: Concordance may be spurious; reconsider transfer modeling approach.

## Dependencies
- **Task 4 (Cell-Type Concordance)**: Must complete to get per-KD concordance scores
- **Task 3 (Baseline Expression)**: Optional, for filtering genes

## Next Steps
- **If validation succeeds**: High confidence to proceed with Task 6 (Transfer Model)
- **If validation fails**: May need pathway-specific models or alternative approaches
- Findings inform which pathways are safe to transfer from K562 to RPE1
