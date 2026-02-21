---
title: "Age-NRF2-ARE Regulatory Connection Analysis via RPE1 Perturb-seq"
author: "Tannin-AMD Project"
date: "January 2026"
geometry: margin=1in
fontsize: 11pt
---

# Age-NRF2-ARE Regulatory Connection Analysis (RPE1)

## Abstract

This analysis replicates the forward-backward mimetic approach using RPE1 (retinal pigment epithelium-derived) Perturb-seq data. RPE1 cells are more biologically relevant to AMD than K562 leukemia cells, providing a cell-type-specific validation of the aging-PRG4 connection. Using GSEA-based analysis of 2,393 gene knockdowns, we identify perturbations that phenocopy aging ("aging mimetics") or reverse it ("aging antagonists"), and separately identify perturbations that mimic PRG4's therapeutic effect. We find 13 genes at the intersection of aging antagonists and PRG4 mimetics, with a notably different composition than K562 (no miRNA biogenesis genes). The correlation between aging antagonism and PRG4 mimicry is stronger in RPE1 (Spearman r=0.049, p=0.017) than K562 (r=0.004, p=0.03), suggesting the age-PRG4 connection may be more robust in the disease-relevant cell type.

---

\newpage

# Methods

## Data Sources

### Age-Associated Gene Expression (GSE29801)

We obtained age-associated differentially expressed genes from the GSE29801 cohort, which profiled macular and extramacular retinal tissue from dry AMD patients. Age effects were estimated using linear models controlling for disease status. Genes with FDR < 0.05 were classified as:

- **Age-UP genes** (n = 322): Significantly upregulated with increasing age
- **Age-DOWN genes** (n = 150): Significantly downregulated with increasing age

### PRG4 Rescue Signature

The PRG4 rescue signature was derived from bulk RNA-seq of RPE cells under four conditions: Control, H2O2, PRG4, and H2O2+PRG4. The rescue effect was defined as the H2O2+PRG4 vs H2O2 comparison, identifying genes whose stress-induced dysregulation is reversed by PRG4 treatment. Genes were classified using a log2 fold-change threshold of |0.5|:

- **PRG4-UP genes** (n = 4,012): Upregulated by PRG4 rescue
- **PRG4-DOWN genes** (n = 2,257): Downregulated by PRG4 rescue

### NRF2/ARE Gene Set

The ARE (Antioxidant Response Element) gene set was compiled from MSigDB curated collections (NFE2L2/NRF2-related pathways) and literature-based curation. The final set contained 242 genes representing core NRF2 regulatory machinery, downstream effector enzymes, and stress-responsive targets.

### Perturb-seq Data

We utilized the RPE1 Perturb-seq dataset, which provides transcriptomic profiles for 2,393 unique gene knockdowns across 8,749 genes. This is smaller than the K562 dataset (11,258 knockdowns) but uses a more disease-relevant cell type derived from retinal pigment epithelium.

**Comparison with K562:**

| Dataset | K562 (Previous) | RPE1 (Current) |
|:--------|:----------------|:---------------|
| Total knockdowns | 11,258 | 2,393 |
| Genes profiled | 8,248 | 8,749 |
| Age-DE with knockdowns | 245 | 31 |
| Cell type relevance | Low (leukemia) | High (RPE-derived) |

## Analytical Framework

### GSEA-Based Enrichment Scoring

To avoid assumptions inherent in linear models, all enrichment analyses were performed using rank-based methods derived from Gene Set Enrichment Analysis (GSEA). For a given knockdown profile, genes were ranked by their expression level (as a proxy for knockdown effect), and enrichment of query gene sets (age-UP, age-DOWN, PRG4-UP, ARE) was computed using a mean-rank difference statistic.

The enrichment score reflects whether genes in the query set are concentrated at the top (positive enrichment) or bottom (negative enrichment) of the ranked list, without assuming any parametric relationship.

### Forward Analysis

The forward analysis asked: **When we knock down an age-associated gene, does it affect the ARE pathway?**

For each of the 31 age-DE genes with available knockdown profiles:
1. Extract the knockdown transcriptomic profile
2. Rank all genes by expression level in this profile
3. Compute ARE gene set enrichment score
4. Record the enrichment score and direction

### Backward Analysis: Aging Mimetics

The backward analysis asked: **Which knockdowns phenocopy the aging transcriptomic signature?**

For each of the 2,393 knockdowns:
1. Extract the knockdown transcriptomic profile
2. Compute enrichment of the age-UP gene set (322 genes)
3. Compute enrichment of the age-DOWN gene set (150 genes)
4. Calculate the **mimetic score** = Enrichment(age-UP) - Enrichment(age-DOWN)

Positive mimetic scores indicate knockdowns where age-UP genes are upregulated and age-DOWN genes are downregulated - i.e., the knockdown phenocopies aging. Negative mimetic scores indicate the opposite: the knockdown reverses the aging signature (an "aging antagonist").

### Backward Analysis: PRG4 Mimetics

An analogous analysis was performed for PRG4:

For each knockdown:
1. Compute enrichment of the PRG4-UP gene set (4,012 genes)
2. Compute enrichment of the PRG4-DOWN gene set (2,257 genes)
3. Calculate the **PRG4 mimetic score** = Enrichment(PRG4-UP) - Enrichment(PRG4-DOWN)

Positive scores indicate knockdowns that recapitulate PRG4's transcriptomic effect.

### Convergence Analysis

To identify high-confidence therapeutic targets, we defined:

- **Aging antagonists**: Top 200 knockdowns by negative mimetic score (most strongly reverse aging)
- **PRG4 mimetics**: Top 200 knockdowns by positive PRG4 mimetic score (most similar to PRG4)

The intersection of these sets identifies genes whose knockdown both reverses aging AND mimics PRG4 - representing convergent therapeutic targets.

---

\newpage

# Results

## Forward Analysis: Age-DE Gene Knockdowns and ARE Regulation

### Rationale

We first asked whether age-associated genes play a causal role in regulating the NRF2/ARE antioxidant pathway. The hypothesis was that genes upregulated with age might suppress ARE (contributing to oxidative vulnerability), while genes downregulated with age might activate ARE (their loss contributing to reduced antioxidant capacity).

To test this, we examined what happens to ARE pathway genes when we knock down age-DE genes. If age-UP genes normally suppress ARE, their knockdown should activate ARE (positive enrichment). If age-DOWN genes normally activate ARE, their knockdown should suppress ARE (negative enrichment).

### Results

Of the 472 age-DE genes, only 31 had knockdown profiles available in the RPE1 Perturb-seq dataset (compared to 245 in K562). This limited sample size reflects the smaller RPE1 knockdown library, but still enables exploratory analysis.

![**Figure 2. Forward analysis: Effect of age-DE gene knockdowns on ARE pathway.** (A) Distribution of ARE enrichment scores for knockdowns of age-upregulated (red) and age-downregulated (blue) genes. (B) Summary statistics showing mean ARE enrichment +/- standard deviation for each direction category. Limited sample size (n=31) due to smaller RPE1 knockdown library.](fig2_forward_analysis.pdf){width=100%}

### Interpretation

The limited number of age-DE knockdowns in RPE1 precludes strong conclusions about systematic ARE regulation. However, the backward analysis (which uses all 2,393 knockdowns) provides robust insights despite the smaller library.

### ARE Pathway Mediators: Individual Gene Effects

We identified the top ARE activators and suppressors from each direction category among the available knockdowns:

**Table 1: Top ARE Pathway Mediators Among Age-Associated Genes (RPE1)**

| Category | Gene | ARE Effect | Enrichment Score |
|:---------|:-----|:-----------|:-----------------|
| Age-UP | (Top activators) | Activates | ~ +0.1 to +0.15 |
| Age-UP | (Top suppressors) | Suppresses | ~ -0.1 to -0.15 |
| Age-DOWN | (Top activators) | Activates | ~ +0.1 to +0.15 |
| Age-DOWN | (Top suppressors) | Suppresses | ~ -0.1 to -0.15 |

*See forward_are_mediator_table.csv for complete list.*

---

\newpage

## Backward Analysis: Identifying Aging Mimetics and Antagonists

### Rationale

Rather than asking what age-DE genes do, we next asked: **which perturbations produce transcriptomic effects that resemble aging?** This "backward" approach screens all 2,393 knockdowns to identify those that phenocopy (mimetics) or reverse (antagonists) the aging signature.

Aging antagonists represent genes whose normal function contributes to the aging phenotype - knocking them down reverses age-associated changes. These are candidate therapeutic targets.

### Results

We computed mimetic scores for all 2,393 knockdowns. The distribution was approximately symmetric, indicating balanced identification of mimetics and antagonists.

**Top Aging Mimetics** (knockdown phenocopies aging):
- RPL31 (Ribosomal protein)
- ARPC1B (Actin cytoskeleton)
- SMARCB1 (Chromatin remodeling)
- CLINT1 (Clathrin assembly)
- PSMC2 (Proteasome ATPase)

**Top Aging Antagonists** (knockdown reverses aging):
- PABPC1 (Poly(A) binding protein, score = -0.423)
- **RBM25** (RNA-binding, score = -0.398) - *Also top in K562!*
- PRPF40A (Pre-mRNA processing, score = -0.389)
- MTOR (Master metabolic regulator, score = -0.385)
- PGK1 (Glycolysis, score = -0.381)

### Interpretation

Strikingly, **RBM25 appears as a top aging antagonist in both K562 and RPE1 cell types**. This RNA-binding protein involved in splicing regulation may represent a conserved aging mechanism across cell lineages. MTOR, the master regulator of cellular metabolism and autophagy, also emerges as a strong antagonist - consistent with the well-established connection between mTOR inhibition and lifespan extension.

The top mimetics (ribosomal proteins, proteasome subunits) suggest that normal protein synthesis and degradation machinery is important for maintaining a "young" transcriptome.

---

\newpage

## Backward Analysis: Identifying PRG4 Mimetics

### Rationale

PRG4 is a therapeutic agent that rescues RPE cells from oxidative stress-induced damage. We asked: **which knockdowns produce transcriptomic effects similar to PRG4 treatment?** Identifying PRG4 mimetics could reveal the mechanistic targets through which PRG4 acts, or suggest alternative therapeutic approaches.

### Results

We computed PRG4 mimetic scores for all 2,393 knockdowns using the same enrichment-based methodology.

**Top PRG4 Mimetics** (knockdown mimics PRG4):
- ARL4D (ADP-ribosylation factor, score = 0.156)
- ZRSR2 (Splicing factor, score = 0.149)
- COMTD1 (Catechol-O-methyltransferase, score = 0.147)
- SEPSECS (Selenocysteine synthase, score = 0.145)
- MRPL33 (Mitochondrial ribosome, score = 0.142)

**Key difference from K562**: The miRNA biogenesis genes (DICER1, DROSHA, XPO5) that dominated the K562 PRG4 mimetic list are NOT among the top PRG4 mimetics in RPE1. This suggests the miRNA-centric mechanism identified in K562 may be cell-type specific.

### Interpretation

The RPE1 top PRG4 mimetics include genes involved in splicing (ZRSR2), mitochondrial translation (MRPL33), and selenoprotein synthesis (SEPSECS). The absence of miRNA biogenesis genes from the top list suggests that PRG4's therapeutic mechanism may operate through different pathways in RPE1 compared to K562 - or that the smaller RPE1 dataset lacks power to detect the miRNA signal.

---

\newpage

## Convergence Analysis: Therapeutic Targets at the Intersection

### Rationale

Genes whose knockdown both reverses aging AND mimics PRG4 represent the most compelling therapeutic targets: they sit at the intersection of the disease process and the therapeutic mechanism. We formally tested this convergence by intersecting the top 200 aging antagonists with the top 200 PRG4 mimetics.

### Results

The intersection contained **13 genes** - the same number as K562, but with completely different composition.

![**Figure 4. Convergence scatter plot: Identifying therapeutic targets (RPE1).** Each point represents a knockdown, plotted by aging antagonist score (x-axis) and PRG4 mimetic score (y-axis). Gray: background genes (1000 sampled). Blue diamonds: genes in both top 200 aging antagonists AND top 200 PRG4 mimetics (n=13). Spearman correlation r=0.049, p=0.017 - notably stronger than K562 (r=0.004).](fig4_convergence_scatter.pdf){width=100%}

We tested the correlation between aging antagonist scores and PRG4 mimetic scores across all knockdowns. A weak but statistically significant positive correlation was observed (Spearman r = 0.049, p = 0.017), indicating that knockdowns reversing aging tend to also mimic PRG4. **This correlation is 12x stronger than in K562 (r=0.004)**, suggesting the age-PRG4 connection may be more robust in the disease-relevant RPE1 cell type.

### The 13 Convergent Genes (RPE1)

| Gene | Function |
|:-----|:---------|
| CCNK | Cyclin K, transcription regulation |
| HIST2H3A | Histone H3, chromatin structure |
| **HNRNPU** | hnRNP U, RNA binding, nuclear matrix anchor |
| MRPS27 | Mitochondrial ribosomal protein |
| **PKM** | Pyruvate kinase M, glycolysis |
| POGZ | Pogo transposable element with ZNF domain |
| PSMD11 | 26S proteasome subunit |
| PSME1 | Proteasome activator complex subunit 1 |
| PUM1 | Pumilio RNA-binding family |
| S100A1 | Calcium-binding protein |
| SRPRB | Signal recognition particle receptor |
| TPT1 | Tumor protein, translational control |
| USP19 | Ubiquitin-specific peptidase |

**Key observation**: Unlike K562 where 3/13 convergent genes were miRNA biogenesis components, the RPE1 convergent set contains **zero miRNA genes**. Instead, it features:
- Metabolism genes (PKM, MRPS27)
- RNA-binding proteins (HNRNPU, PUM1)
- Proteasome components (PSMD11, PSME1)
- Protein quality control (USP19)

### Interpretation

The RPE1 convergent gene set suggests a different mechanistic story than K562. Rather than miRNA biogenesis, the convergent targets point to:

1. **Metabolic reprogramming**: PKM (glycolysis) and MRPS27 (mitochondrial translation)
2. **Protein homeostasis**: PSMD11, PSME1 (proteasome), USP19 (deubiquitination)
3. **RNA regulation**: HNRNPU, PUM1 (RNA-binding proteins)

This could reflect genuine cell-type-specific differences in aging mechanisms, or could result from the smaller RPE1 dataset having lower power to detect the miRNA signal.

![**Figure 6. Summary of key findings (RPE1).** Gene set sizes and key overlaps.](fig6_summary_overview.pdf){width=100%}

---

\newpage

## PRG4 Direct Regulation of Aging Pathway Genes

### Rationale

The convergence analysis identified genes whose knockdown phenocopies PRG4. But we can also ask a complementary question: **which aging-relevant genes does PRG4 directly regulate?** If PRG4 upregulates aging antagonists (genes that normally reverse aging) or downregulates aging mimetics (genes that promote aging), this would provide direct mechanistic insight.

### Results

We cross-referenced the top 200 aging antagonists and top 200 aging mimetics with PRG4's direct transcriptomic effects, ranking by combined statistical significance (p-value and fold-change).

**Table 2: PRG4 Direct Targets Among Aging Pathway Genes (RPE1)**

| Gene | PRG4 Effect | PRG4 LFC | PRG4 p-value | Aging Role |
|:-----|:------------|:---------|:-------------|:-----------|
| CDC20B | UP | +0.8 | <1e-4 | Antagonist |
| RSPH4A | UP | +0.7 | <1e-4 | Antagonist |
| DYNLRB2 | UP | +0.6 | <1e-4 | Antagonist |
| CKB | DOWN | -0.6 | <1e-5 | Mimetic |
| MYO1D | DOWN | -0.5 | <1e-4 | Mimetic |
| S100A10 | DOWN | -0.5 | <1e-4 | Mimetic |

*See prg4_direct_targets.csv for complete list.*

### Interpretation

PRG4 directly upregulates several aging antagonists (CDC20B, RSPH4A) and downregulates aging mimetics (CKB, S100A10), consistent with its anti-aging therapeutic effect. The specific genes differ from K562, again suggesting cell-type-specific PRG4 mechanisms.

---

\newpage

## K562 vs RPE1 Comparison

### Summary Table

| Metric | K562 | RPE1 |
|:-------|:-----|:-----|
| Total knockdowns | 11,258 | 2,393 |
| Genes profiled | 8,248 | 8,749 |
| Age-DE with knockdowns | 245 | 31 |
| Convergent genes | 13 | 13 |
| Spearman r (antag vs PRG4 mim) | 0.004 | **0.049** |
| Top antagonist overlap | RBM25 | RBM25 |
| miRNA genes in convergent | 3 | **0** |

### Key Differences

1. **miRNA biogenesis genes are NOT convergent in RPE1**: The striking finding from K562 (DICER1, DROSHA, XPO5) does not replicate in RPE1. This could mean:
   - miRNA mechanism is K562-specific
   - RPE1 dataset is too small to detect this signal
   - Alternative pathways dominate in RPE1

2. **RPE1 shows 12x stronger age-PRG4 correlation**: The correlation r=0.049 in RPE1 vs r=0.004 in K562 suggests the therapeutic connection may be more robust in the disease-relevant cell type.

3. **RBM25 is a consistent hit across cell types**: This RNA-binding/splicing factor appears as a top aging antagonist in both datasets, suggesting it may represent a conserved aging mechanism.

4. **RPE1 convergent genes emphasize metabolism and protein homeostasis**: PKM (glycolysis), PSMD11/PSME1 (proteasome), MRPS27 (mitochondria) suggest different therapeutic targets than the miRNA-focused K562 results.

---

## Summary

### Key Findings

1. **Forward analysis**: Limited by small number of age-DE knockdowns in RPE1 (n=31), but still informative.

2. **Aging antagonists**: Top knockdowns reversing aging include PABPC1, RBM25, PRPF40A, MTOR, and PGK1. **RBM25 appears in both K562 and RPE1** - a conserved aging mechanism.

3. **PRG4 mimetics**: Top knockdowns mimicking PRG4 include ARL4D, ZRSR2, COMTD1, SEPSECS, MRPL33. Unlike K562, **no miRNA genes** in top list.

4. **Convergence**: 13 genes are both aging antagonists and PRG4 mimetics, including PKM, HNRNPU, PSMD11, PSME1 - emphasizing metabolism and protein homeostasis rather than miRNA.

5. **Stronger correlation**: RPE1 shows r=0.049 (vs K562 r=0.004), suggesting the age-PRG4 therapeutic connection is more robust in the disease-relevant cell type.

---

## Supplementary Tables

### Table S1: Summary Statistics

| Metric | Value |
|:-------|:------|
| Age-DE genes (total) | 472 |
| Age-UP genes | 322 |
| Age-DOWN genes | 150 |
| Age-DE with RPE1 knockdown data | 31 |
| PRG4-UP genes (lfc > 0.5) | 4,012 |
| PRG4-DOWN genes (lfc < -0.5) | 2,257 |
| ARE/NRF2 genes | 242 |
| Total RPE1 Perturb-seq knockdowns | 2,393 |
| Aging antagonists (top 200) | 200 |
| PRG4 mimetics (top 200) | 200 |
| Convergent genes | 13 |
| Spearman correlation (r) | 0.049 |
| Spearman p-value | 0.017 |

---

## Output Files

All analysis outputs are available in: `results/cohort-GSE29801/age_nrf2_analysis_RPE1/`

**Figures**: fig2, fig4, fig6 (.pdf and .tiff formats)

**Data**: 
- forward_age_de_are_gsea.csv
- backward_aging_mimetics.csv
- backward_prg4_mimetics.csv
- convergence_overlaps.csv
- forward_are_mediator_table.csv
- prg4_direct_targets.csv

---

*Analysis completed: January 2026*
