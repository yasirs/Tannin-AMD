# Comprehensive Perturb-seq Gene Set Analysis Report

**Analysis Date:** January 17, 2025  
**Output Directory:** `results/perturbseq-analysis/`

---

## Executive Summary

This report documents the analysis of **13 gene sets** across **2 Perturb-seq datasets** (RPE1 Essential, K562 GWPS). The analysis identifies genetic knockdowns that mimic or antagonize age-associated, AMD-associated, and pathobiological (Senescence, Redox, Inflammation) gene expression signatures.

### Key Findings

1. **399 convergent genes** appear as top antagonists across multiple analyses
2. **Significant RPE1-K562 concordance** for Age (21 genes, p<1e-15) and Redox (11 genes, p<5.7e-6)
3. **High gene coverage** for Redox-PRO (94%) and Senescence-ANTI (92%)
4. **Lower coverage for AMD gene sets** (7-45%), reflecting tissue-specific expression

---

## Methods

### Gene Sets Analyzed

| Category | Gene Set | Size | Description |
|----------|----------|------|-------------|
| **Age** | Age-UP | 323 | Upregulated with age (GSE29801) |
| | Age-DOWN | 161 | Downregulated with age (GSE29801) |
| | ARE | 242 | NRF2/ARE targets (MSigDB + literature) |
| **AMD** | AMD-Macula-UP | 126 | Upregulated in AMD RPE Macula (p<0.01) |
| | AMD-Macula-DOWN | 431 | Downregulated in AMD RPE Macula (p<0.01) |
| | AMD-nonMacula-UP | 395 | Upregulated in AMD RPE non-Macula (p<0.01) |
| | AMD-nonMacula-DOWN | 508 | Downregulated in AMD RPE non-Macula (p<0.01) |
| **Senescence** | Senescence-PRO | 250 | Pro-senescence (SenMayo + GO + SASP) |
| | Senescence-ANTI | 327 | Anti-senescence (Cell cycle + E2F) |
| **Redox** | Redox-PRO | 245 | Pro-oxidative stress (HALLMARK + GO) |
| | Redox-ANTI | 16 | Antioxidant response (NRF2 targets) |
| **Inflammation** | Inflammation-PRO | 571 | Pro-inflammatory (HALLMARK + NFkB) |
| | Inflammation-ANTI | 111 | Anti-inflammatory (GO regulatory) |

### Perturb-seq Datasets

| Dataset | Perturbations | Genes | Cell Type |
|---------|---------------|-------|-----------|
| RPE1 Essential | 2,679 | 8,749 | RPE1 (human RPE cells) |
| K562 GWPS | 11,258 | 8,248 | K562 (erythroleukemia) |

### Analytical Approach

**Backward Analysis (Mimetic Scoring):**
For each knockdown, we computed a mimetic score:
$$\text{Mimetic Score} = \text{Enrichment}(\text{UP genes}) - \text{Enrichment}(\text{DOWN genes})$$

- **Positive scores**: Knockdown upregulates UP genes and downregulates DOWN genes (mimetic)
- **Negative scores**: Knockdown has opposite effect (antagonist)

---

## Results

### Gene Set Coverage in Perturb-seq Data

| Gene Set | RPE1 Coverage | K562 Coverage |
|----------|---------------|---------------|
| Redox-PRO | 94.7% (232/245) | 93.9% (230/245) |
| Senescence-ANTI | 91.7% (300/327) | 92.7% (303/327) |
| ARE | 62.4% (151/242) | 52.9% (128/242) |
| Age-DOWN | 49.7% (80/161) | 46.6% (75/161) |
| Inflammation-PRO | 46.9% (268/571) | 35.9% (205/571) |
| Age-UP | 43.7% (141/323) | 36.5% (118/323) |
| AMD-Macula-UP | 45.2% (57/126) | 34.1% (43/126) |
| Senescence-PRO | 43.2% (108/250) | 31.6% (79/250) |
| Inflammation-ANTI | 35.1% (39/111) | 29.7% (33/111) |
| AMD-nonMacula-UP | 32.7% (129/395) | 24.6% (97/395) |
| AMD-nonMacula-DOWN | 14.8% (75/508) | 12.4% (63/508) |
| AMD-Macula-DOWN | 7.4% (32/431) | 7.7% (33/431) |
| Redox-ANTI | 62.5% (10/16) | 56.3% (9/16) |

> [!NOTE]
> Low coverage for AMD-DOWN gene sets reflects tissue-specific expression in RPE that differs from cell line expression profiles.

![Results Coverage Summary](concordance/Fig_coverage_summary.pdf)

---

### Analysis Summary Statistics

| Gene Set Pair | Dataset | UP Genes | DOWN Genes | Mean Score | Std | Range |
|---------------|---------|----------|------------|------------|-----|-------|
| Age | RPE1 | 141 | 80 | -0.019 | 0.114 | [-0.39, 0.31] |
| Age | K562 | 118 | 75 | -0.027 | 0.150 | [-0.46, 0.38] |
| AMD-Macula | RPE1 | 57 | 32 | 0.174 | 0.228 | [-0.51, 0.60] |
| AMD-Macula | K562 | 43 | 33 | 0.021 | 0.220 | [-0.55, 0.67] |
| AMD-nonMacula | RPE1 | 129 | 75 | 0.013 | 0.141 | [-0.39, 0.36] |
| AMD-nonMacula | K562 | 97 | 63 | -0.001 | 0.152 | [-0.40, 0.40] |
| Senescence | RPE1 | 108 | 300 | 0.350 | 0.295 | [-0.58, 0.83] |
| Senescence | K562 | 79 | 303 | 0.039 | 0.147 | [-0.40, 0.67] |
| Redox | RPE1 | 232 | 10 | 0.007 | 0.281 | [-0.65, 0.60] |
| Redox | K562 | 230 | 9 | 0.006 | 0.282 | [-0.80, 0.68] |
| Inflammation | RPE1 | 268 | 39 | 0.014 | 0.122 | [-0.34, 0.38] |
| Inflammation | K562 | 205 | 33 | -0.005 | 0.170 | [-0.43, 0.48] |

---

### RPE1 vs K562 Concordance

| Gene Set | Overlap (top 200) | Jaccard | Hypergeom p-value | Spearman r | Status |
|----------|-------------------|---------|-------------------|------------|--------|
| **Age** | 21 | 0.055 | **1.0e-15** | 0.153 | YES Significant |
| **Redox** | 11 | 0.028 | **5.7e-6** | 0.144 | YES Significant |
| **Inflammation** | 9 | 0.023 | **1.9e-4** | -0.001 | YES Overlap significant |
| **AMD-Macula** | 6 | 0.015 | 0.015 | 0.064 | Marginal |
| Senescence | 5 | 0.013 | 0.051 | 0.264 | Not significant |
| AMD-nonMacula | 2 | 0.005 | 0.597 | 0.011 | Not significant |

> [!IMPORTANT]
> Age, Redox, and Inflammation show significant cross-dataset concordance, suggesting these pathways are conserved across cell types.

![Concordance Heatmap](concordance/Fig_concordance_heatmap.pdf)

#### Concordance Scatter Plots

**Age Concordance**
![Age concordance scatter](concordance/Fig_Age_rpe1_vs_k562_scatter.pdf)

**Redox Concordance**
![Redox concordance scatter](concordance/Fig_Redox_rpe1_vs_k562_scatter.pdf)

**Inflammation Concordance**
![Inflammation concordance scatter](concordance/Fig_Inflammation_rpe1_vs_k562_scatter.pdf)

---

### Convergent Genes

We identified **399 genes** that appear as top antagonists in multiple analyses (>=2 appearances). Top convergent genes include knockdowns that score highly across multiple pathways and/or both datasets.

**Genes appearing in >=6 analyses:**
- See `convergent_genes_all.csv` for complete list

![Convergent Genes Ranked](concordance/Fig_convergent_genes_ranked.pdf)

---

## Output Files

### Per-Analysis Outputs

Each analysis produces the following files in its respective directory:

```
results/perturbseq-analysis/
|-- Axes/
|   |-- Senescence/
|   |   |-- RPE1/backward/
|   |   |   |-- backward_Senescence_mimetics.csv
|   |   |   |-- top_200_mimetics.csv
|   |   |   |-- top_200_antagonists.csv
|   |   |   |-- Fig_Senescence_distribution.pdf
|   |   |   `-- Fig_Senescence_top_genes.pdf
|   |   `-- K562_GWPS/backward/[same structure]
|   |-- Redox/[same structure]
|   `-- Inflammation/[same structure]
|-- AMD/
|   |-- RPE1/backward/[AMD-Macula, AMD-nonMacula results]
|   `-- K562_GWPS/backward/[same]
|-- concordance/
|   |-- rpe1_vs_k562_concordance.csv
|   |-- cross_geneset_overlaps_RPE1.csv
|   |-- cross_geneset_overlaps_K562_GWPS.csv
|   `-- convergent_genes_all.csv
`-- gene-sets/
    |-- registry.csv
    |-- coverage_in_perturbseq.csv
    `-- [individual gene set CSVs]
```

---

## Interpretation

### Age Signature
The Age signature analysis identifies knockdowns that reverse age-associated gene expression changes. The significant RPE1-K562 concordance (21 genes, p<1e-15) suggests conserved aging mechanisms across cell types.

**Age Analysis Results (RPE1)**
![Age Distribution](Age/RPE1/backward/Fig_Age_distribution.pdf)
![Age Top Genes](Age/RPE1/backward/Fig_Age_top_genes.pdf)

### Senescence Axis
High coverage of Senescence-ANTI genes (92%) but lower Senescence-PRO coverage (32-43%) provides good detection of anti-senescence programs. The non-significant cross-dataset overlap may reflect cell-type-specific senescence programs.

**Senescence Analysis Results (RPE1)**
![Senescence Distribution](Axes/Senescence/RPE1/backward/Fig_Senescence_distribution.pdf)
![Senescence Top Genes](Axes/Senescence/RPE1/backward/Fig_Senescence_top_genes.pdf)

### Redox Axis
Excellent coverage of Redox-PRO genes (94%) but limited Redox-ANTI coverage (56-63%) due to the small NRF2 target gene set. The significant cross-dataset concordance (11 genes, p<5.7e-6) indicates conserved oxidative stress response.

**Redox Analysis Results (RPE1)**
![Redox Distribution](Axes/Redox/RPE1/backward/Fig_Redox_distribution.pdf)
![Redox Top Genes](Axes/Redox/RPE1/backward/Fig_Redox_top_genes.pdf)

### Inflammation Axis
Good coverage of Inflammation-PRO genes (36-47%) with significant cross-dataset overlap (9 genes, p<1.9e-4). The near-zero correlation (r=-0.001) suggests that while some antagonists overlap, the overall ranking differs between cell types.

**Inflammation Analysis Results (RPE1)**
![Inflammation Distribution](Axes/Inflammation/RPE1/backward/Fig_Inflammation_distribution.pdf)
![Inflammation Top Genes](Axes/Inflammation/RPE1/backward/Fig_Inflammation_top_genes.pdf)

### AMD Signatures
Low coverage for AMD gene sets reflects tissue-specific expression patterns in RPE that differ from cell line expression profiles. AMD-specific findings should be interpreted cautiously.

**AMD-Macula Analysis Results (RPE1)**
![AMD Macula Distribution](AMD/RPE1/backward/Fig_AMD-Macula_distribution.pdf)
![AMD Macula Top Genes](AMD/RPE1/backward/Fig_AMD-Macula_top_genes.pdf)

---

## Conclusions

1. **Conserved pathways**: Age, Redox, and Inflammation pathways show significant cross-dataset concordance
2. **Cell-type specificity**: AMD and Senescence signatures show more cell-type-specific patterns
3. **Convergent targets**: 399 genes appear across multiple analyses as potential intervention targets
4. **Coverage limitations**: AMD-DOWN gene sets have low coverage, limiting sensitivity

---

## Technical Notes

- **Analysis runtime**: 7.5 minutes
- **Code location**: `code/perturbseq-analysis/`
- **Python environment**: `~/venvs/torch`
- **Key files**: `run_analysis.py` (orchestration), `enrichment_analysis.py` (GSEA scoring)
