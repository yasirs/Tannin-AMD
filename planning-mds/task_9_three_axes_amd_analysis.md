# Task 9: Three-Axes AMD Pathobiology Analysis

## Objective
Systematically demonstrate that PRG4 comprehensively reverses the three core pathobiological axes of AMD: **Oxidative Stress**, **Inflammation**, and **Senescence**. Produce publication-ready figures with detailed captions and quantitative validation.

## Scientific Rationale

AMD pathobiology comprises three interconnected axes:

| Axis | Disease State (↑ in AMD) | Healthy State (PRG4 Rescue) |
|:-----|:-------------------------|:----------------------------|
| **1. Oxidative/Redox** | ROS accumulation, oxidative damage, mitochondrial dysfunction | Antioxidant response, NRF2 activation, ROS detoxification |
| **2. Inflammatory** | Complement activation, cytokine release, NFκB signaling | Anti-inflammatory, resolution, complement regulation |
| **3. Senescence** | Cell cycle arrest, SASP, p16/p21 induction | Anti-senescence, proliferative capacity, cellular fitness |

**Central Hypothesis**: PRG4 treatment reverses all three axes, providing comprehensive protection against AMD-like phenotypes.

**H2O2 as Dual Model**: H2O2 treatment induces both oxidative stress AND stress-induced premature senescence (SIPS), making it a clinically relevant model for AMD pathobiology.

---

## Data Sources

### Internal Experimental Data
- **Bulk RNA-seq DEGs**: `data/RPE_cells/code/RPE_gene pvals.xlsx`
  - Contrasts: H2O2_vs_CTRL, PRG4_vs_CTRL, H2O2PRG4_vs_H2O2
  - Threshold: FDR < 0.05 (3,975 DEGs for PRG4 rescue)

- **TPM Expression**: `data/RPE_cells/code/RPE_TPMS.xlsx`
  - For heatmaps and z-score normalization

### External Validation Cohorts
- **GSE135092**: AMD vs Control (537 samples, RNA-seq)
  - DEG results: `results/cohort-GSE135092/GSE135092_DE_results.csv`
- **GSE29801**: AMD vs Control (293 samples, microarray)

### Perturb-seq Datasets
- **RPE1 Essential**: `data/external/perturbseq/rpe1_normalized_bulk_01.h5ad` (2,679 KDs)
- **K562 Essential**: `data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad` (2,285 KDs)
- **K562 GWPS**: `data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad` (11,258 KDs)

### Gene Set Sources

| Axis | Gene Set | Source | Size |
|:-----|:---------|:-------|:-----|
| **Senescence** | SenMayo | Saul et al. 2022 (PMID: 35974106), MSigDB | 125 genes |
| **Senescence** | Gold standards | CDKN1A, CDKN2A (excluded from SenMayo) | 2 genes |
| **Senescence** | SASP subset | Curated from SenMayo | ~40 genes |
| **Senescence** | GO:0090398 | Cellular Senescence | Variable |
| **Oxidative** | HALLMARK_OXIDATIVE_PHOSPHORYLATION | MSigDB Hallmarks | ~200 genes |
| **Oxidative** | NRF2 targets | KEAP1-NRF2 pathway genes | ~50 genes |
| **Oxidative** | GO:0006979 | Response to Oxidative Stress | Variable |
| **Oxidative** | GO:0034599 | Cellular Response to Oxidative Stress | Variable |
| **Inflammatory** | HALLMARK_INFLAMMATORY_RESPONSE | MSigDB Hallmarks | ~200 genes |
| **Inflammatory** | HALLMARK_COMPLEMENT | MSigDB Hallmarks | ~200 genes |
| **Inflammatory** | NFκB targets | Literature-curated | ~100 genes |
| **Inflammatory** | GO:0006954 | Inflammatory Response | Variable |
| **Known AMD** | GWAS genes | IAMDGC 2016 | 50 genes |
| **Known AMD** | Literature AMD markers | CFH, RPE65, BEST1, etc. | ~20 genes |

---

## Analysis Tasks

### Part 1: Gene Set Curation and Preparation

#### 1.1 Create Axis Gene Sets
**Script**: `code/three-axes/01_curate_gene_sets.py`

For each axis, curate comprehensive gene sets:

**Oxidative Stress Axis**:
- Download/extract: HALLMARK_OXIDATIVE_PHOSPHORYLATION, NRF2 targets, GO:0006979
- Create: `oxidative_pro_disease.csv` (pro-oxidant genes)
- Create: `oxidative_anti_disease.csv` (antioxidant genes, NRF2 targets)

**Inflammatory Axis**:
- Download/extract: HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_COMPLEMENT, NFκB targets
- Create: `inflammatory_pro_disease.csv` (pro-inflammatory)
- Create: `inflammatory_anti_disease.csv` (anti-inflammatory, resolution)

**Senescence Axis**:
- Download SenMayo from MSigDB (SAUL_SEN_MAYO)
- Add CDKN1A, CDKN2A as gold standards
- Extract SASP subset (IL6, IL8, MMP3, CXCL8, SERPINE1, etc.)
- Create: `senescence_pro_disease.csv` (senescence markers, SASP)
- Create: `senescence_anti_disease.csv` (anti-senescence, proliferative)

**Known AMD Genes**:
- Extract from GWAS report: CFH, C3, ARMS2, HTRA1, APOE, TIMP3, etc.
- Add literature markers: RPE65, BEST1, RLBP1, TTR
- Create: `known_amd_genes.csv`

**Output Directory**: `results/three-axes/gene-sets/`

#### 1.2 Map Gene Sets to Data
**Script**: `code/three-axes/02_map_gene_sets.py`

- Map gene symbols to Ensembl IDs
- Filter to genes present in bulk RNA-seq data
- Filter to genes present in each Perturb-seq dataset
- Create coverage summary tables

**Output**: 
- `results/three-axes/gene-sets/coverage_summary.csv`
- Filtered gene set files with Ensembl IDs

---

### Part 2: Axis Expression Analysis

#### 2.1 H2O2 Induces All Three Axes
**Script**: `code/three-axes/03_h2o2_axis_induction.R`

**Question**: Does H2O2 treatment recapitulate AMD pathobiology by inducing all three axes?

**Method**:
1. For each axis gene set, calculate:
   - Number of genes induced by H2O2 (log2FC > 0.5, FDR < 0.05)
   - Mean log2FC across axis genes
   - Enrichment test (hypergeometric)
2. Validate with canonical markers:
   - Oxidative: NQO1, HMOX1, GCLC (NRF2 targets should be induced as stress response)
   - Inflammatory: IL6, IL1B, CXCL8, CFH (should change with stress)
   - Senescence: CDKN1A, CDKN2A, SERPINE1, IL6 (should be induced)

**Output Figures**:
- `Fig_H2O2_axis_induction_barplot.pdf`: Bar chart showing % genes induced per axis
- `Fig_H2O2_canonical_markers.pdf`: Expression of canonical markers across conditions

**Output Data**:
- `results/three-axes/h2o2-induction/h2o2_axis_summary.csv`
- `results/three-axes/h2o2-induction/h2o2_induced_genes_per_axis.csv`

#### 2.2 PRG4 Reverses All Three Axes
**Script**: `code/three-axes/04_prg4_axis_reversal.R`

**Question**: Does PRG4 treatment reverse the AMD-like phenotype induced by H2O2?

**Method**:
1. For each axis gene set:
   - Identify genes induced by H2O2 AND reversed by PRG4 (opposite direction)
   - Calculate reversal rate: (H2O2-induced genes reversed by PRG4) / (total H2O2-induced)
   - Calculate mean PRG4 rescue effect (log2FC H2O2+PRG4 vs H2O2)
2. Directional analysis:
   - Pro-disease genes: Should be ↑ by H2O2, ↓ by PRG4
   - Anti-disease genes: Should be ↓ by H2O2, ↑ by PRG4

**Output Figures**:
- `Fig_PRG4_reversal_rates.pdf`: Reversal rate per axis (bar chart with error bars)
- `Fig_PRG4_GSEA_Enrichment.pdf`: GSEA enrichment plot of Axis genes in PRG4 ranked list (replacing scatter correlation)

**Output Data**:
- `results/three-axes/prg4-reversal/prg4_reversal_summary.csv`
- `results/three-axes/prg4-reversal/prg4_reversed_genes_per_axis.csv`

#### 2.3 Axis Heatmaps
**Script**: `code/three-axes/05_axis_heatmaps.R`

**Question**: How do all axis genes behave across all experimental conditions?

**Method**:
1. Extract TPM values for all axis genes
2. Convert to gene-wise z-scores
3. Create heatmaps with:
   - Rows: Genes (grouped by axis, ordered by H2O2 effect)
   - Columns: Samples (CTRL, H2O2, PRG4, H2O2+PRG4)
   - Color: BrBG (brown = high)
   - Annotations: Axis membership, PRG4 reversal status

**Output Figures** (one per axis + one combined):
- `Fig_Heatmap_oxidative_axis.pdf`
- `Fig_Heatmap_inflammatory_axis.pdf`
- `Fig_Heatmap_senescence_axis.pdf`
- `Fig_Heatmap_three_axes_combined.pdf` (tri-band visualization)

**Output Data**:
- `results/three-axes/heatmaps/heatmap_data_zscore.csv`

---

### Part 3: Axis Scores and Quantitative Metrics

#### 3.1 GSEA-Based Axis Metrics
**Script**: `code/three-axes/06_axis_scores.R`

**Question**: Is the gene set as a whole significantly upregulated or downregulated by each condition?

**Method**:
1. Perform GSEA for each axis gene set against the full ranked list of genes for each contrast:
   - H2O2 vs CTRL
   - PRG4 vs CTRL
   - H2O2+PRG4 vs H2O2
2. **Metric**: Use the **Normalized Enrichment Score (NES)** as the primary measure of pathway activity change.
3. **Significance**: Use GSEA FDR/p-value.
4. Report NES and FDR for all axes across all contrasts.

**Output Figures**:
- `Fig_Axis_scores_boxplot.pdf`: Boxplot of disease scores by condition (3 panels)
- `Fig_Axis_scores_radar.pdf`: Radar/spider plot showing all 3 axes per condition

**Output Data**:
- `results/three-axes/axis-scores/sample_axis_scores.csv`
- `results/three-axes/axis-scores/axis_score_statistics.csv`

#### 3.2 Enrichment Analysis
**Script**: `code/three-axes/07_axis_enrichment.R`

**Question**: Are PRG4-rescued genes significantly enriched for each axis?

**Method**:
1. Hypergeometric test: PRG4-upregulated genes vs each axis gene set
2. GSEA: Rank PRG4 rescue signature, test enrichment of each axis
3. Compare enrichment across axes

**Output Figures**:
- `Fig_Enrichment_barplot.pdf`: -log10(p-value) for each axis
- `Fig_GSEA_curves.pdf`: GSEA enrichment curves for each axis

**Output Data**:
- `results/three-axes/enrichment/enrichment_results.csv`
- `results/three-axes/enrichment/gsea_results.csv`

---

### Part 4: AMD Cohort Validation

#### 4.1 Axis Expression in AMD Patients
**Script**: `code/three-axes/08_amd_cohort_validation.R`

**Question**: Are the three axes dysregulated in AMD patients, and does PRG4 rescue oppose the AMD signature?

**Method** (for GSE135092 and GSE29801):
1. **Scoring**: Calculate **GSVA scores** (or ssGSEA) for each axis gene set for every patient sample. This provides a robust, rank-based metric per individual.
2. **Analysis**:
   - Compare GSVA scores between AMD and Control (Wilcoxon test).
   - Correlate GSVA scores with clinical metadata (if available).
   - Correlate AMD Status (binary) with GSVA score.

**Output Figures**:
- `Fig_AMD_axis_dysregulation.pdf`: Volcano plots annotated with axis genes
- `Fig_AMD_PRG4_GSEA.pdf`: GSEA enrichment plot of PRG4 rescue signature in ranked AMD gene list (replacing scatter correlation)

**Output Data**:
- `results/three-axes/amd-validation/amd_axis_expression.csv`
- `results/three-axes/amd-validation/amd_prg4_correlation.csv`

#### 4.2 Known AMD Gene Panel
**Script**: `code/three-axes/09_known_amd_genes.R`

**Question**: Do canonical AMD genes follow the expected pattern (dysregulated in AMD, rescued by PRG4)?

**Method**:
1. Select known AMD genes: CFH, C3, RPE65, BEST1, RLBP1, HTRA1, ARMS2, etc.
2. For each gene, show:
   - Expression in bulk RNA-seq (CTRL, H2O2, PRG4, H2O2+PRG4)
   - Expression in AMD cohorts (Control vs AMD)
   - PRG4 rescue direction
3. Create dedicated figure panel

**Output Figures**:
- `Fig_Known_AMD_genes_heatmap.pdf`: Heatmap of known AMD genes across all datasets
- `Fig_Known_AMD_genes_barplot.pdf`: Bar plot of expression per gene per condition
- `Fig_CFH_RPE65_spotlight.pdf`: Detailed panels for top genes

**Output Data**:
- `results/three-axes/known-amd/known_amd_gene_expression.csv`

---

### Part 5: Perturb-seq Integration

#### 5.1 Approach A: Axis Regulators
**Script**: `code/three-axes/10_perturbseq_axis_regulators.py`

**Question**: Which gene knockdowns activate or suppress each axis?

**Method**:
1. For each Perturb-seq dataset (RPE1, K562 Essential, K562 GWPS):
   - For each knockdown, calculate axis score (mean Z-score of axis genes)
   - Rank knockdowns by axis score
   - Identify top 50 "pro-axis" KDs (activate axis) and top 50 "anti-axis" KDs (suppress axis)
2. Cross-cell-type comparison:
   - Identify regulators found in both RPE1 and K562 (high confidence)

**Output Figures**:
- `Fig_Axis_regulators_rankings.pdf`: Top 20 regulators per axis, per cell type
- `Fig_Axis_regulators_concordance.pdf`: Scatter plot comparing RPE1 vs K562

**Output Data**:
- `results/three-axes/perturbseq/axis_regulators_RPE1.csv`
- `results/three-axes/perturbseq/axis_regulators_K562_essential.csv`
- `results/three-axes/perturbseq/axis_regulators_K562_gwps.csv`
- `results/three-axes/perturbseq/axis_regulators_concordant.csv`

#### 5.2 Approach B: Axis-Specific PRG4 Mimetics
**Script**: `code/three-axes/11_perturbseq_prg4_mimetics.py`

**Question**: Which knockdowns mimic PRG4's effect on each axis?

**Method**:
1. Define PRG4 rescue gene set for each axis.
2. For each knockdown, rank all genes by expression change.
3. Calculate **GSEA NES** of the PRG4 rescue set in the KD ranked list.
4. Rank knockdowns by NES (high positive NES = Mimetic).

**Output Figures**:
- `Fig_PRG4_mimetics_per_axis.pdf`: Top 20 mimetics per axis ranked by NES
- `Fig_PRG4_mimetics_pathway.pdf`: Pathway enrichment of mimetics per axis

**Output Data**:
- `results/three-axes/perturbseq/prg4_mimetics_oxidative.csv`
- `results/three-axes/perturbseq/prg4_mimetics_inflammatory.csv`
- `results/three-axes/perturbseq/prg4_mimetics_senescence.csv`

#### 5.3 Approach C: GWAS Gene Knockdown Effects
**Script**: `code/three-axes/12_perturbseq_gwas_genes.py`

**Question**: Do AMD GWAS genes regulate the three axes when knocked down?

**Method**:
1. Identify GWAS genes with knockdowns in Perturb-seq (~6-10 genes)
2. For each GWAS gene KD, calculate effect on each axis
3. Connect: Genetic risk → Axis dysregulation → PRG4 rescue

**Output Figures**:
- `Fig_GWAS_gene_axis_effects.pdf`: Heatmap (GWAS genes x axes, color = effect size)
- `Fig_CFH_KD_axis_profile.pdf`: Detailed profile of CFH knockdown effects

**Output Data**:
- `results/three-axes/perturbseq/gwas_gene_axis_effects.csv`

#### 5.4 Regulator-PRG4 Connection
**Script**: `code/three-axes/13_regulator_prg4_connection.R`

**Question**: Does PRG4 affect expression of the identified axis regulators?

**Method**:
1. Take top regulators from Approach A (e.g., top 20 per axis)
2. Check their expression in bulk RNA-seq:
   - Are they differentially expressed by PRG4?
   - Are they induced by H2O2?
3. Create mechanistic model: PRG4 → Regulator → Axis

**Output Figures**:
- `Fig_Regulator_expression_changes.pdf`: Heatmap of regulator expression
- `Fig_Mechanistic_model.pdf`: Network diagram of PRG4 → Regulators → Axes

**Output Data**:
- `results/three-axes/perturbseq/regulator_prg4_expression.csv`

---

### Part 6: SASP/Senescence Deep Dive

#### 6.1 SenMayo Validation
**Script**: `code/three-axes/14_senmayo_validation.R`

**Question**: Does H2O2 induce the senescence signature, and does PRG4 reverse it?

**Method**:
1. Use full SenMayo gene set (125 genes) + CDKN1A/CDKN2A
2. Create detailed heatmap with gene annotations
3. Calculate SenMayo score per sample
4. Perform GSEA on PRG4 rescue signature

**Output Figures**:
- `Fig_SenMayo_heatmap.pdf`: Full SenMayo gene set heatmap
- `Fig_SenMayo_score.pdf`: SenMayo score by condition
- `Fig_SASP_subset.pdf`: SASP-specific genes highlighted

**Output Data**:
- `results/three-axes/senescence/senmayo_expression.csv`
- `results/three-axes/senescence/senmayo_scores.csv`

#### 6.2 Senescence Markers Spotlight
**Script**: `code/three-axes/15_senescence_markers.R`

**Question**: Do canonical senescence markers (p21, p16, IL6) follow expected patterns?

**Method**:
1. Select gold-standard markers: CDKN1A, CDKN2A, SERPINE1, IL6, MMP3
2. Show individual gene expression across all conditions
3. Compare with GO:Cellular Senescence enrichment

**Output Figures**:
- `Fig_Senescence_markers_barplot.pdf`: Individual gene expression bar plots
- `Fig_p21_p16_spotlight.pdf`: Detailed panels for CDKN1A and CDKN2A

---

### Part 7: Integrated Visualization

#### 7.1 Tri-Axis Summary Figure
**Script**: `code/three-axes/16_summary_figure.R`

**Question**: Can we create a single comprehensive figure showing all three axes?

**Method**:
1. Create multi-panel figure:
   - Panel A: Three-band heatmap (Oxidative | Inflammatory | Senescence)
   - Panel B: Radar plot of axis scores
   - Panel C: Reversal rate summary
   - Panel D: Known AMD gene validation

**Output Figures**:
- `Fig_Three_Axes_Summary.pdf`: Publication-ready multi-panel figure
- `Fig_Three_Axes_Summary.tiff`: High-resolution TIFF (300 dpi, LZW+P)

#### 7.2 Mechanism Summary
**Script**: `code/three-axes/17_mechanism_summary.R`

**Question**: What is the integrated mechanistic model?

**Method**:
1. Compile top regulators from all three axes
2. Cross-reference with PRG4 expression effects
3. Create network diagram showing:
   - PRG4 → [Regulators] → [Axes] → [Rescue phenotype]

**Output Figures**:
- `Fig_Mechanism_network.pdf`: Network diagram
- `Fig_Mechanism_pathway.pdf`: Pathway-centric view

---

## Output Organization

### Directory Structure
```
results/three-axes/
├── gene-sets/                    # Curated gene sets
│   ├── oxidative_pro_disease.csv
│   ├── oxidative_anti_disease.csv
│   ├── inflammatory_pro_disease.csv
│   ├── inflammatory_anti_disease.csv
│   ├── senescence_pro_disease.csv
│   ├── senescence_anti_disease.csv
│   ├── known_amd_genes.csv
│   └── coverage_summary.csv
│
├── h2o2-induction/               # H2O2 axis induction analysis
│   ├── h2o2_axis_summary.csv
│   ├── h2o2_induced_genes_per_axis.csv
│   ├── Fig_H2O2_axis_induction_barplot.pdf
│   └── Fig_H2O2_canonical_markers.pdf
│
├── prg4-reversal/                # PRG4 reversal analysis
│   ├── prg4_reversal_summary.csv
│   ├── prg4_reversed_genes_per_axis.csv
│   ├── Fig_PRG4_reversal_rates.pdf
│   └── Fig_PRG4_bidirectional_effect.pdf
│
├── heatmaps/                     # Expression heatmaps
│   ├── heatmap_data_zscore.csv
│   ├── Fig_Heatmap_oxidative_axis.pdf
│   ├── Fig_Heatmap_inflammatory_axis.pdf
│   ├── Fig_Heatmap_senescence_axis.pdf
│   └── Fig_Heatmap_three_axes_combined.pdf
│
├── axis-scores/                  # Quantitative axis scores
│   ├── sample_axis_scores.csv
│   ├── axis_score_statistics.csv
│   ├── Fig_Axis_scores_boxplot.pdf
│   └── Fig_Axis_scores_radar.pdf
│
├── enrichment/                   # Enrichment analysis
│   ├── enrichment_results.csv
│   ├── gsea_results.csv
│   ├── Fig_Enrichment_barplot.pdf
│   └── Fig_GSEA_curves.pdf
│
├── amd-validation/               # AMD cohort validation
│   ├── amd_axis_expression.csv
│   ├── amd_prg4_correlation.csv
│   ├── Fig_AMD_axis_dysregulation.pdf
│   └── Fig_AMD_PRG4_correlation.pdf
│
├── known-amd/                    # Known AMD genes
│   ├── known_amd_gene_expression.csv
│   ├── Fig_Known_AMD_genes_heatmap.pdf
│   ├── Fig_Known_AMD_genes_barplot.pdf
│   └── Fig_CFH_RPE65_spotlight.pdf
│
├── perturbseq/                   # Perturb-seq integration
│   ├── axis_regulators_RPE1.csv
│   ├── axis_regulators_K562_essential.csv
│   ├── axis_regulators_K562_gwps.csv
│   ├── axis_regulators_concordant.csv
│   ├── prg4_mimetics_oxidative.csv
│   ├── prg4_mimetics_inflammatory.csv
│   ├── prg4_mimetics_senescence.csv
│   ├── gwas_gene_axis_effects.csv
│   ├── regulator_prg4_expression.csv
│   ├── Fig_Axis_regulators_rankings.pdf
│   ├── Fig_Axis_regulators_concordance.pdf
│   ├── Fig_PRG4_mimetics_per_axis.pdf
│   ├── Fig_PRG4_mimetics_pathway.pdf
│   ├── Fig_GWAS_gene_axis_effects.pdf
│   ├── Fig_CFH_KD_axis_profile.pdf
│   ├── Fig_Regulator_expression_changes.pdf
│   └── Fig_Mechanistic_model.pdf
│
├── senescence/                   # Senescence deep dive
│   ├── senmayo_expression.csv
│   ├── senmayo_scores.csv
│   ├── Fig_SenMayo_heatmap.pdf
│   ├── Fig_SenMayo_score.pdf
│   ├── Fig_SASP_subset.pdf
│   ├── Fig_Senescence_markers_barplot.pdf
│   └── Fig_p21_p16_spotlight.pdf
│
├── summary/                      # Summary figures
│   ├── Fig_Three_Axes_Summary.pdf
│   ├── Fig_Three_Axes_Summary.tiff
│   ├── Fig_Mechanism_network.pdf
│   └── Fig_Mechanism_pathway.pdf
│
└── three_axes_detailed_report.md # Comprehensive report
```

### Code Directory
```
code/three-axes/
├── 01_curate_gene_sets.py
├── 02_map_gene_sets.py
├── 03_h2o2_axis_induction.R
├── 04_prg4_axis_reversal.R
├── 05_axis_heatmaps.R
├── 06_axis_scores.R
├── 07_axis_enrichment.R
├── 08_amd_cohort_validation.R
├── 09_known_amd_genes.R
├── 10_perturbseq_axis_regulators.py
├── 11_perturbseq_prg4_mimetics.py
├── 12_perturbseq_gwas_genes.py
├── 13_regulator_prg4_connection.R
├── 14_senmayo_validation.R
├── 15_senescence_markers.R
├── 16_summary_figure.R
├── 17_mechanism_summary.R
└── utils/
    ├── gene_set_utils.py
    ├── axis_scoring.R
    └── plotting_theme.R
```

---

## Report Standards

### Figure Caption Format (Nature-style)

Each figure caption must follow this template:

> **Figure X. [Descriptive Title]**
> 
> **(a)** [Panel description]. To [objective/question], we [method with sufficient detail]. [Key finding from this panel].
>
> **(b)** [Next panel...]
>
> **Data**: [Source datasets, sample sizes]. **Statistics**: [Test used, p-value thresholds, multiple testing correction]. **Gene sets**: [Source and size of gene sets used].

### Visualization Standards

**All figures must follow the standards defined in [`planning-mds/coding_preferences.md`](file:///home/ysuhail/work/Tannin-AMD/planning-mds/coding_preferences.md)**:

| Element | Specification |
|:--------|:--------------|
| **Language** | R with ggplot2 |
| **Font** | Palatino Linotype |
| **Theme** | `theme_classic()` with `axis.text = element_text(color = "black")` |
| **Figure Size** | Small physical dimensions (inches) |
| **File Formats** | PDF + TIFF (LZW+P compression, 300 dpi) |
| **Heatmaps** | Gene-wise z-scores from TPM; no dendrograms unless needed |
| **Color Palettes** | BrBG for diverging (brown = high); YlOrB for one-sided positive |
| **Code Chunks** | Separated by `#%%` |
| **Caching** | Store heavy computation results in dedicated subfolders |

### Example Caption

> **Figure 3. PRG4 Comprehensively Reverses the Senescence Axis Induced by Oxidative Stress**
> 
> **(a)** Heatmap of SenMayo senescence signature genes across experimental conditions. To determine whether H2O2-induced stress recapitulates cellular senescence and whether PRG4 rescues this phenotype, we extracted 125 genes from the SenMayo curated senescence gene set (Saul et al., 2022) and visualized their expression (gene-wise z-scored TPM) across CTRL, H2O2, PRG4, and H2O2+PRG4 conditions (n=3 biological replicates each). H2O2 treatment induced 78/125 (62.4%) senescence genes (log2FC > 0.5, FDR < 0.05), including CDKN1A (+1.8 log2FC), SERPINE1 (+2.1 log2FC), and IL6 (+2.4 log2FC). PRG4 co-treatment reversed 61/78 (78.2%) of these genes toward baseline levels.
> 
> **(b)** SASP sub-signature quantification. The Senescence-Associated Secretory Phenotype (SASP) subset (42 genes including IL6, IL8, MMP3, CXCL8) was scored per sample as the mean z-score. H2O2 treatment increased SASP score by 1.8-fold compared to CTRL (p < 0.001, two-tailed t-test); PRG4 co-treatment reduced the stress-induced SASP score by 45% (p < 0.01, two-tailed t-test).
> 
> **Data**: Internal RPE bulk RNA-seq (n=12 samples total, 3 per condition). **Gene set**: SenMayo (SAUL_SEN_MAYO, MSigDB C2:CGP, 125 genes). **Statistics**: Two-tailed Student's t-test, FDR correction (Benjamini-Hochberg) for differential expression.

---

## Execution Order

### Phase 1: Gene Set Preparation (Scripts 01-02)
```bash
# Run in sequence
python code/three-axes/01_curate_gene_sets.py
python code/three-axes/02_map_gene_sets.py
```

### Phase 2: Core Expression Analysis (Scripts 03-07)
```bash
# Can run in parallel after Phase 1
Rscript code/three-axes/03_h2o2_axis_induction.R
Rscript code/three-axes/04_prg4_axis_reversal.R
Rscript code/three-axes/05_axis_heatmaps.R
Rscript code/three-axes/06_axis_scores.R
Rscript code/three-axes/07_axis_enrichment.R
```

### Phase 3: AMD Validation (Scripts 08-09)
```bash
# Can run in parallel
Rscript code/three-axes/08_amd_cohort_validation.R
Rscript code/three-axes/09_known_amd_genes.R
```

### Phase 4: Perturb-seq Integration (Scripts 10-13)
```bash
# Run in sequence within phase
python code/three-axes/10_perturbseq_axis_regulators.py
python code/three-axes/11_perturbseq_prg4_mimetics.py
python code/three-axes/12_perturbseq_gwas_genes.py
Rscript code/three-axes/13_regulator_prg4_connection.R
```

### Phase 5: Senescence Deep Dive (Scripts 14-15)
```bash
Rscript code/three-axes/14_senmayo_validation.R
Rscript code/three-axes/15_senescence_markers.R
```

### Phase 6: Summary (Scripts 16-17)
```bash
# Run after all prior phases complete
Rscript code/three-axes/16_summary_figure.R
Rscript code/three-axes/17_mechanism_summary.R
```

### Phase 7: Report Generation
Generate `results/three-axes/three_axes_detailed_report.md` with all figure captions and interpretations.

---

## Success Criteria

### Quantitative Metrics
- [ ] H2O2 induces ≥50% of genes in each axis (log2FC > 0.5, FDR < 0.05)
- [ ] PRG4 reverses ≥60% of H2O2-induced genes per axis
- [ ] Axis scores significantly different between conditions (p < 0.05)
- [ ] GSEA enrichment significant for all three axes (FDR < 0.1)
- [ ] ≥50% of known AMD genes show expected pattern

### Qualitative Criteria
- [ ] All figures have complete, publication-ready captions
- [ ] Heatmaps clearly show three-axis structure
- [ ] Perturb-seq analysis identifies biologically plausible regulators
- [ ] Mechanism model is coherent and testable

### Output Deliverables
- [ ] 30+ figures in `results/three-axes/` subdirectories
- [ ] 15+ data CSV files with axis-level results
- [ ] Comprehensive report (`three_axes_detailed_report.md`)
- [ ] Summary figure suitable for publication

---

## Dependencies

### Required Prior Tasks
- **Task 1 (Robustness)**: DEG threshold selection (FDR < 0.05) ✅
- **Task 2 (Coverage)**: Dataset selection (K562 GWPS) ✅
- **Task 7 (Virtual Screen)**: PRG4 mimetic identification ✅
- **Task 8 (Human Cohorts)**: AMD cohort DEG results ✅

### Software Dependencies
- **R**: ggplot2, pheatmap, clusterProfiler, fgsea, ggpubr
- **Python**: scanpy, pandas, numpy, scipy, matplotlib, seaborn

### Data Dependencies
- Gene set files from MSigDB (will be downloaded)
- SenMayo paper supplementary data (will be extracted)

---

## Notes and Considerations

### H2O2 as Senescence Model
H2O2 induces **stress-induced premature senescence (SIPS)**, characterized by:
- p21/p16 upregulation
- SASP activation
- Cell cycle arrest
- SA-β-gal activity

This makes H2O2 an appropriate model for both oxidative stress AND senescence axes.

### SASP vs Inflammation Overlap
SASP (Senescence-Associated Secretory Phenotype) overlaps heavily with inflammatory genes:
- Shared: IL6, IL8, CXCL8, MMP3
- SASP-specific: SERPINE1, GDF15
- Inflammation-specific: IL1B, TNF, NFκB targets

**Resolution**: Treat SASP as a sub-axis of senescence, but acknowledge overlap in report.

### Cell-Type Considerations
- **RPE1 Perturb-seq**: Most relevant cell type, but limited to essential genes
- **K562 Perturb-seq**: Broader coverage, different cell type
- **Strategy**: Report both, highlight concordant findings

---

## Version History

| Version | Date | Changes |
|:--------|:-----|:--------|
| 1.0 | 2026-01-10 | Initial planning document |

---

**Prepared by:** Gemini Agent  
**Last Updated:** January 10, 2026
