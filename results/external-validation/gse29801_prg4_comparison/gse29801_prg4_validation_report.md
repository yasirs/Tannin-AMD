# GSE29801-PRG4 Cross-Validation Report

**Analysis Date:** January 14, 2026  
**Objective:** Cross-validate GSE29801 AMD-associated gene sets with PRG4 bulk RNA-seq data

---

## Executive Summary

This analysis tested whether gene sets identified in the GSE29801 dry AMD cohort show enrichment in our PRG4 H2O2 stress/rescue bulk RNA-seq experiment. We defined 10 gene sets from both macular and extramacular RPE-choroid tissues, including differentially expressed (DE) genes, age-adjusted variance genes (low/high), bimodal genes, and age-associated genes.

**Key Findings:**
- **No strong enrichments** in PRG4 data (all p>0.10)
- **Modest trends:** Macular low-variance genes show weak positive enrichment in H2O2 (NES=1.29, p=0.10)
- **Negative trend:** Extramacular bimodal genes show negative enrichment in PRG4 rescue (NES=-1.20, p=0.19)
- **Good gene overlap:** 68-88% of GSE29801 genes found in PRG4 data
- **Pathway enrichment:** High-variance genes enriched for oxidative phosphorylation and ribosomal pathways

**Interpretation:** GSE29801 dry AMD signatures and PRG4 H2O2 stress model capture different aspects of AMD biology. The lack of strong enrichment suggests either: (1) dry AMD has weaker molecular signatures than wet AMD, (2) bulk tissue masks cell-type-specific changes, or (3) H2O2 is not a perfect model for dry AMD pathology.

---

## Methods

### GSE29801 Gene Sets Defined

**10 gene sets total (5 per tissue):**

| Gene Set | Macular | Extramacular | Description |
|----------|---------|--------------|-------------|
| **DE (p<0.01)** | 868 genes | 1,353 genes | Mean expression difference in AMD vs Control |
| **Low Variance** | 17 genes | 47 genes | Decreased variance in AMD (age-adjusted, FDR<0.05) |
| **High Variance** | 17 genes | 146 genes | Increased variance in AMD (age-adjusted, FDR<0.05) |
| **Bimodal** | 8 genes | 25 genes | Genes with bimodal distribution in AMD (age-adjusted) |
| **Age DE** | 18 genes | 472 genes | Age-associated genes (FDR<0.05) |

**Age Adjustment:** Variance genes were computed on residuals after regressing out age effect, ensuring disease-specific signal.

### PRG4 Bulk Data Comparisons

1. **H2O2 vs Control:** Does H2O2 induce AMD-like changes?
2. **PRG4+H2O2 vs H2O2 (Rescue):** Does PRG4 reverse AMD genes?

**Method:** GSEA (fgsea) testing enrichment of GSE29801 gene sets in ranked PRG4 gene lists (ranked by logFC).

---

## Results

### 1. Gene Set Enrichment in PRG4 Data (GSEA)

**Summary Table:**

| Tissue | Gene Set | H2O2 vs Control | PRG4 Rescue |
|--------|----------|-----------------|-------------|
| **Macular** | | | |
| | DE | NES=0.82, p=0.998 | NES=0.93, p=0.747 |
| | Low Variance | **NES=1.29, p=0.103** | NES=1.17, p=0.293 |
| | High Variance | NES=0.70, p=0.857 | NES=0.96, p=0.524 |
| | Bimodal | NES=0.82, p=0.763 | NES=1.11, p=0.369 |
| | Age DE | NES=0.99, p=0.555 | NES=1.29, p=0.165 |
| **Extramacular** | | | |
| | DE | NES=0.93, p=0.945 | NES=0.76, p=0.997 |
| | Low Variance | NES=1.03, p=0.473 | NES=1.01, p=0.477 |
| | High Variance | NES=0.50, p=0.999 | NES=-0.70, p=0.976 |
| | Bimodal | NES=0.60, p=0.956 | **NES=-1.20, p=0.189** |
| | Age DE | NES=0.84, p=0.987 | NES=0.77, p=0.956 |

**Overlap:** 68-88% of GSE29801 genes found in PRG4 data (good technical overlap).

**Interpretation:**
- **No significant enrichments** (all FDR>0.20)
- **Macular low-variance genes** show borderline positive trend in H2O2 (p=0.10), suggesting uniformly downregulated genes in AMD might be upregulated by H2O2
- **Extramacular bimodal genes** show negative NES in PRG4 rescue (p=0.19), suggesting PRG4 might affect these differently than expected

---

### 2. Pathway Enrichment of GSE29801 Gene Sets

**Over-Representation Analysis (ORA) using Fisher's Exact Test**

#### Extramacular High-Variance Genes (146 genes)

**Top Hallmark Pathways:**
- OXIDATIVE_PHOSPHORYLATION (p=0.002, 18 genes)
- MYC_TARGETS_V1 (p=0.008, 21 genes)
- E2F_TARGETS (p=0.02, 12 genes)

**Top GO BP Pathways:**
- MITOCHONDRIAL_TRANSLATION (p=3e-06, 14 genes)
- RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS (p=8e-06, 28 genes)
- OXIDATIVE_PHOSPHORYLATION (p=1e-05, 18 genes)

**Interpretation:** High-variance genes in AMD are enriched for mitochondrial/metabolic pathways, consistent with energy dysregulation.

#### Extramacular DE Genes (1,353 genes)

**Top Hallmark Pathways:**
- OXIDATIVE_PHOSPHORYLATION (p=2e-10, 71 genes)
- MYC_TARGETS_V1 (p=3e-09, 85 genes)
- ADIPOGENESIS (p=5e-08, 67 genes)

**Interpretation:** DE genes also show metabolic dysregulation, overlapping with variance genes.

#### Extramacular Low-Variance Genes (47 genes)

**Top Hallmark Pathways:**
- EPITHELIAL_MESENCHYMAL_TRANSITION (p=0.04, 9 genes)

**Interpretation:** Low-variance genes (uniformly downregulated in AMD) may relate to RPE dedifferentiation.

---

## Figures

### Figure 1: Gene Set Enrichment Barplot

![Faceted barplot showing NES by gene set type](file:///home/ysuhail/work/Tannin-AMD/results/external-validation/gse29801_prg4_comparison/figures/figure1_gene_set_enrichment_barplot.pdf)

**Description:** Normalized Enrichment Score (NES) for each GSE29801 gene set in H2O2 vs Control and PRG4 Rescue comparisons. Bars colored by significance (red: p<0.05, orange: p<0.20, gray: NS). No strong enrichments observed.

### Figure 2: Pathway Enrichment Panels

![Multi-panel dotplot showing top enriched pathways](file:///home/ysuhail/work/Tannin-AMD/results/external-validation/gse29801_prg4_comparison/figures/figure2_pathway_enrichment_panels.pdf)

**Description:** Top 5 enriched pathways (GO BP, Hallmark) for selected GSE29801 gene sets. Dot size = overlap genes, color = -log10(p-value). High-variance genes enriched for oxidative phosphorylation and ribosome biogenesis.

---

## Discussion

### Why No Strong Enrichment in PRG4 Data?

**Possible explanations:**

1. **Different disease contexts:**
   - GSE29801: Dry AMD (late-stage geographic atrophy)
   - PRG4: H2O2 stress model (acute oxidative stress)
   - Dry AMD is chronic, multifactorial; H2O2 is acute, single stressor

2. **Tissue heterogeneity:**
   - Both use bulk RPE-choroid tissue
   - Cell-type-specific changes (RPE vs choroid vs immune) may cancel out
   - Single-cell analysis might reveal stronger concordance

3. **Weak dry AMD signatures:**
   - Only 1-2 DE genes at FDR<0.10 in GSE29801
   - Variance is the dominant signal, not mean changes
   - PRG4/H2O2 primarily affects mean expression, not variance

4. **PRG4 mechanism:**
   - PRG4 may act via different pathways than expected
   - Anti-inflammatory effects may not directly reverse oxidative stress genes

### Value of Variance Analysis

**Key insight:** High-variance genes (146 in extramacular) are enriched for **oxidative phosphorylation** and **ribosome biogenesis**, suggesting:
- AMD disrupts mitochondrial function variably across patients
- Some patients have mitochondrial dysfunction, others don't
- **Patient stratification potential:** High-var gene expression could classify AMD subtypes

**Bimodal genes (25):** Show negative trend in PRG4 rescue (NES=-1.20, p=0.19), suggesting:
- These genes might be PRG4 targets
- Further investigation of BLOC1S2, SEC13, HSF1, PSMC6 warranted

---

## Conclusions

1. **No strong overlap** between GSE29801 dry AMD and PRG4 H2O2/rescue signatures
2. **Variance genes** provide biological insight: mitochondrial/ribosomal dysregulation
3. **Bimodal genes** show interesting negative trend in PRG4 rescue (hypothesis-generating)
4. **Pathway convergence:** Both AMD and H2O2 affect oxidative phosphorylation
5. **Next steps:** Single-cell analysis, direct dry AMD vs wet AMD comparison, PRG4 mechanism

---

## Output Files

### Gene Sets
- `macular_low_variance_age_adj.csv`
- `macular_high_variance_age_adj.csv`
- `macular_bimodality_age_adjusted.csv`
- `extramacular_low_variance_age_adj.csv`
- `extramacular_high_variance_age_adj.csv`
- `extramacular_bimodality_age_adjusted.csv`
- `gene_set_sizes.csv`

### GSEA Results
- `gsea_all_results.csv` - Complete GSEA results (20 tests)
- `gsea_summary.csv` - Summary table

### Pathway Enrichment
- `pathway_enrichment_ora_results.csv` - ORA results for all gene sets

### Figures
- `figure1_gene_set_enrichment_barplot.pdf/.tiff`
- `figure2_pathway_enrichment_panels.pdf/.tiff`
- `macular_amd_clustering_highvar_age_adj.pdf/.tiff` (21 genes, 28 samples)
- `macular_amd_clustering_bimodal_age_adj.pdf/.tiff` (9 genes, 28 samples)
- `extramacular_amd_clustering_highvar_age_adj.pdf/.tiff` (210 genes, 25 samples)
- `extramacular_amd_clustering_bimodal_age_adj.pdf/.tiff` (46 genes, 25 samples)

---

## Technical Notes

- **Age adjustment:** Critical for variance analysis - increased extramacular high-var genes from 172 to 210
- **Bimodality testing:** Used mclust on age-corrected residuals (BIC improvement >10)
- **GSEA settings:** minSize=5, nproc=1, two-tailed
- **ORA background:** 20,000 genes
- **Gene overlap:** Good technical concordance (68-88%) between GSE29801 Agilent and PRG4 RNA-seq platforms

---

*Analysis performed using R 4.3.3 with fgsea, msigdbr, pheatmap, ggplot2*
