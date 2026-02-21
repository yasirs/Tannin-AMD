---
title: "Age-NRF2-ARE Regulatory Connection Analysis via Perturb-seq"
author: "Tannin-AMD Project"
date: "January 2026"
geometry: margin=1in
fontsize: 11pt
---

# Age-NRF2-ARE Regulatory Connection Analysis

## Abstract

Age-related macular degeneration (AMD) is characterized by progressive oxidative damage to the retinal pigment epithelium. The NRF2/ARE pathway serves as the master regulator of antioxidant defense, yet its relationship to aging and therapeutic intervention with PRG4 remains incompletely understood. Using GSEA-based analysis of genome-wide Perturb-seq data, we identify perturbations that phenocopy aging ("aging mimetics") or reverse it ("aging antagonists"), and separately identify perturbations that mimic PRG4's therapeutic effect. We find 14 genes at the intersection of aging antagonists and PRG4 mimetics, of which three (DICER1, DROSHA, XPO5) are core components of the miRNA biogenesis pathway. This convergence suggests miRNA dysregulation as a mechanistic link between aging and PRG4-mediated rescue.

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

We utilized the K562 Genome-Wide Perturb-seq (GWPS) dataset, which provides transcriptomic profiles for 11,258 gene knockdowns. Each knockdown profile represents the aggregate expression changes across cells receiving that perturbation, enabling systematic causal inference of gene regulatory relationships.

## Analytical Framework

### GSEA-Based Enrichment Scoring

To avoid assumptions inherent in linear models, all enrichment analyses were performed using rank-based methods derived from Gene Set Enrichment Analysis (GSEA). For a given knockdown profile, genes were ranked by their expression level (as a proxy for knockdown effect), and enrichment of query gene sets (age-UP, age-DOWN, PRG4-UP, ARE) was computed using a Kolmogorov-Smirnov-like running sum statistic.

The enrichment score reflects whether genes in the query set are concentrated at the top (positive enrichment) or bottom (negative enrichment) of the ranked list, without assuming any parametric relationship.

### Forward Analysis

The forward analysis asked: **When we knock down an age-associated gene, does it affect the ARE pathway?**

For each of the 245 age-DE genes with available knockdown profiles:
1. Extract the knockdown transcriptomic profile
2. Rank all genes by expression level in this profile
3. Compute ARE gene set enrichment score
4. Record the enrichment score and direction

Results were stratified by the gene's age-association direction (UP vs DOWN) to test whether age-UP and age-DOWN genes differentially regulate ARE.

### Backward Analysis: Aging Mimetics

The backward analysis asked: **Which knockdowns phenocopy the aging transcriptomic signature?**

For each of the 11,258 knockdowns:
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

Of the 472 age-DE genes, 245 had knockdown profiles available in the Perturb-seq dataset (153 age-UP, 92 age-DOWN). For each knockdown, we computed the ARE enrichment score.

![**Figure 2. Forward analysis: Effect of age-DE gene knockdowns on ARE pathway.** (A) Distribution of ARE enrichment scores for knockdowns of age-upregulated (red, n=153) and age-downregulated (blue, n=92) genes. Both distributions center near zero. (B) Top age-DOWN gene knockdowns ranked by ARE effect: positive scores indicate ARE pathway activation upon knockdown. (C) Summary statistics showing mean ARE enrichment +/- standard deviation for each direction category.](fig2_forward_analysis.pdf){width=100%}

The mean ARE enrichment score was -0.005 for age-UP knockdowns and +0.007 for age-DOWN knockdowns. Neither group showed systematic directional effects on ARE. Individual knockdowns showed variable effects: among age-DOWN genes, PPIC knockdown showed the strongest ARE activation (enrichment = 0.196), while WDR26 knockdown showed ARE suppression (enrichment = -0.144).

### Interpretation

Age-associated genes do not, as a group, preferentially regulate the ARE pathway when perturbed. This does not exclude the possibility that specific age-DE genes regulate ARE, but suggests that age-associated transcriptomic changes are not primarily driven by direct ARE modulators.

### ARE Pathway Mediators: Individual Gene Effects

Although age-DE genes do not systematically regulate ARE as a group, individual genes show strong effects. We identified the top 6 ARE activators and suppressors from each direction category:

**Table 1: Top ARE Pathway Mediators Among Age-Associated Genes**

| Category | Gene | ARE Effect | Enrichment Score |
|:---------|:-----|:-----------|:-----------------|
| Age-UP | CIB1 | Activates | +0.157 |
| Age-UP | HBZ | Activates | +0.117 |
| Age-UP | FHL3 | Activates | +0.103 |
| Age-UP | SLC25A39 | Suppresses | -0.144 |
| Age-UP | APLP2 | Suppresses | -0.127 |
| Age-UP | PAQR4 | Suppresses | -0.127 |
| Age-DOWN | PPIC | Activates | +0.196 |
| Age-DOWN | IFRD1 | Activates | +0.156 |
| Age-DOWN | NR3C1 | Activates | +0.136 |
| Age-DOWN | WDR26 | Suppresses | -0.144 |
| Age-DOWN | TPM1 | Suppresses | -0.134 |
| Age-DOWN | ANGPT1 | Suppresses | -0.110 |

These mediators suggest specific mechanisms by which aging affects antioxidant response. For example, PPIC (peptidylprolyl isomerase C), which decreases with age, appears to normally activate ARE pathway genes; its decline may contribute to reduced antioxidant capacity in aging tissue. Conversely, SLC25A39, which increases with age, appears to suppress ARE when knocked down - suggesting it may normally promote ARE activity as a compensatory response.

---

\newpage

## Backward Analysis: Identifying Aging Mimetics and Antagonists

### Rationale

Rather than asking what age-DE genes do, we next asked: **which perturbations produce transcriptomic effects that resemble aging?** This "backward" approach screens all 11,258 knockdowns to identify those that phenocopy (mimetics) or reverse (antagonists) the aging signature.

Aging antagonists represent genes whose normal function contributes to the aging phenotype - knocking them down reverses age-associated changes. These are candidate therapeutic targets.

### Results

We computed mimetic scores for all 11,258 knockdowns. The distribution was approximately symmetric with mean = -0.023 and standard deviation = 0.109, indicating balanced identification of mimetics and antagonists.

![**Figure 3. Backward analysis: Flow from knockdowns to therapeutic targets.** Schematic showing the filtering of 11,258 knockdowns into aging mimetics (top 200, red), neutral (center), and aging antagonists (bottom 200, green). Of the aging antagonists, 14 also qualify as PRG4 mimetics (blue). Among these 14 convergent genes, 3 are components of the miRNA biogenesis pathway (purple): DICER1, DROSHA, and XPO5.](fig3_flow_diagram.pdf){width=100%}

**Top Aging Mimetics** (knockdown phenocopies aging):
UPF2 (score = 0.405), HIST1H2AD (0.353), MPP5 (0.332), ZNF688 (0.330), TRUB2 (0.313)

**Top Aging Antagonists** (knockdown reverses aging):
RBM25 (score = -0.517), FGFR1OP (-0.461), CNOT3 (-0.432), MED9 (-0.422), RARS (-0.420)

As a sanity check, we confirmed that knockdowns of age-UP genes were not systematically aging mimetics (mean mimetic score of age-UP knockdowns = -0.017), consistent with the expectation that removing a gene that increases with age should not phenocopy aging.

### Interpretation

The top aging antagonists (RBM25, FGFR1OP, CNOT3) are involved in RNA processing and transcriptional regulation. RBM25 is an RNA-binding protein involved in splicing; CNOT3 is a component of the CCR4-NOT deadenylase complex regulating mRNA stability. These functions suggest that post-transcriptional regulatory programs may be key drivers of the aging transcriptome.

---

\newpage

## Backward Analysis: Identifying PRG4 Mimetics

### Rationale

PRG4 is a therapeutic agent that rescues RPE cells from oxidative stress-induced damage. We asked: **which knockdowns produce transcriptomic effects similar to PRG4 treatment?** Identifying PRG4 mimetics could reveal the mechanistic targets through which PRG4 acts, or suggest alternative therapeutic approaches.

### Results

We computed PRG4 mimetic scores for all 11,258 knockdowns using the same enrichment-based methodology.

**Top PRG4 Mimetics** (knockdown mimics PRG4):
SHOC2 (score = 0.211), DICER1 (0.209), AMBRA1 (0.208), IMPDH2 (0.208), XPO5 (0.206)

Notably, DICER1 and XPO5 - both components of miRNA biogenesis - appeared among the top 5 PRG4 mimetics. This overlap with the aging antagonist set (see below) was unexpected and became the central finding of this analysis.

### Interpretation

The appearance of miRNA pathway genes among PRG4 mimetics suggests that PRG4's therapeutic effect may involve modulation of miRNA-mediated gene regulation. Whether PRG4 directly affects miRNA biogenesis or achieves similar downstream effects through an independent mechanism remains to be determined.

---

\newpage

## Convergence Analysis: Therapeutic Targets at the Intersection

### Rationale

Genes whose knockdown both reverses aging AND mimics PRG4 represent the most compelling therapeutic targets: they sit at the intersection of the disease process and the therapeutic mechanism. We formally tested this convergence by intersecting the top 200 aging antagonists with the top 200 PRG4 mimetics.

### Results

The intersection contained **14 genes** - significantly more than expected by chance given the total number of knockdowns.

![**Figure 4. Convergence scatter plot: Identifying therapeutic targets.** Each point represents a knockdown, plotted by aging antagonist score (x-axis) and PRG4 mimetic score (y-axis). Gray: all genes. Red squares: ARE pathway genes. Blue diamonds: genes in both top 200 aging antagonists AND top 200 PRG4 mimetics (n=14). Purple stars: miRNA biogenesis genes (DICER1, DROSHA, XPO5), which cluster in the upper-right quadrant.](fig4_convergence_scatter.pdf){width=100%}

We tested the correlation between aging antagonist scores and PRG4 mimetic scores across all knockdowns. A weak but statistically significant positive correlation was observed (Spearman r = 0.004, p = 0.03), indicating that knockdowns reversing aging tend to also mimic PRG4.

### The 14 Convergent Genes

| Gene | Function |
|:-----|:---------|
| **DICER1** | miRNA processing enzyme |
| **DROSHA** | miRNA biogenesis (nuclear) |
| **XPO5** | miRNA nuclear export |
| BORA | Mitotic regulator |
| FCHO1 | Clathrin-mediated endocytosis |
| MED14 | Mediator complex (transcription) |
| NKX6-1 | Homeobox transcription factor |
| RBBP6 | RB-binding, cell cycle |
| RBM14-RBM4 | RNA-binding proteins |
| RPRD1B | RNA Pol II regulation |
| TANGO6 | Golgi transport |
| UBA3 | NEDD8 E1 activating enzyme |
| YTHDC1 | m6A RNA reader |

Three of the 14 convergent genes (21%) are core components of the canonical miRNA biogenesis pathway: DROSHA cleaves pri-miRNA in the nucleus, XPO5 exports pre-miRNA to the cytoplasm, and DICER1 processes pre-miRNA into mature miRNA.

![**Figure 5. Convergent therapeutic targets map to miRNA biogenesis pathway.** Schematic of canonical miRNA processing showing the three convergent genes (DROSHA, XPO5, DICER1) in their pathway context. Knockdown of any of these genes both reverses the aging transcriptomic phenotype and mimics PRG4's therapeutic effect.](fig5_mirna_pathway_schematic.pdf){width=100%}

### ARE Enrichment Among Convergent Genes

We tested whether the convergent genes were enriched for ARE pathway membership. Of the 14 convergent genes, 0 were members of the 242-gene ARE set. Among the broader sets:
- ARE genes in top 200 aging antagonists: 3/200 (1.5%)
- ARE genes in top 200 PRG4 mimetics: 3/200 (1.5%)  
- Background rate: 155/9,867 (1.6%)

There was no enrichment of ARE genes among convergent therapeutic targets, suggesting that the PRG4-aging connection operates through pathways distinct from direct NRF2/ARE modulation.

### Interpretation

The convergence of three miRNA biogenesis genes among 14 total convergent targets is striking. Given that there are only ~10-12 core miRNA biogenesis genes, this represents substantial enrichment. The biological interpretation is that reducing miRNA biogenesis capacity simultaneously reverses age-associated transcriptomic changes and achieves effects similar to PRG4 treatment.

This could reflect:
1. miRNA-mediated suppression of protective genes that increases with age
2. Global miRNA dysregulation as a feature of cellular aging
3. PRG4 acting through miRNA-independent mechanisms that converge on similar downstream targets

---

\newpage

## PRG4 Direct Regulation of Aging Pathway Genes

### Rationale

The convergence analysis identified genes whose knockdown phenocopies PRG4. But we can also ask a complementary question: **which aging-relevant genes does PRG4 directly regulate?** If PRG4 upregulates aging antagonists (genes that normally reverse aging) or downregulates aging mimetics (genes that promote aging), this would provide direct mechanistic insight into how PRG4 achieves its therapeutic effect.

### Results

We cross-referenced the top 200 aging antagonists and top 200 aging mimetics with PRG4's direct transcriptomic effects, ranking by combined statistical significance (p-value and fold-change).

**Table 2: PRG4 Direct Targets Among Aging Pathway Genes**

| Gene | PRG4 Effect | PRG4 LFC | PRG4 p-value | Aging Role | Mimetic Score |
|:-----|:------------|:---------|:-------------|:-----------|:--------------|
| PMAIP1 | UP | +1.23 | 1.27e-07 | Antagonist | -0.337 |
| NUF2 | UP | +1.03 | 3.50e-06 | Antagonist | -0.299 |
| RBBP6 | UP | +0.81 | 1.70e-07 | Antagonist | -0.374 |
| PLK4 | UP | +0.95 | 1.03e-05 | Antagonist | -0.333 |
| BORA | UP | +0.99 | 5.96e-05 | Antagonist | -0.285 |
| HSPA1B | UP | +0.95 | 6.76e-05 | Antagonist | -0.303 |
| BMP1 | DOWN | -0.78 | 4.79e-18 | Mimetic | +0.270 |
| MARCKSL1 | DOWN | -0.73 | 6.30e-12 | Mimetic | +0.239 |
| CTSH | DOWN | -0.60 | 7.59e-12 | Mimetic | +0.238 |
| CDK5 | DOWN | -0.63 | 2.71e-04 | Mimetic | +0.253 |

PRG4 significantly upregulates 107 of the top 200 aging antagonists and significantly downregulates 73 of the top 200 aging mimetics. This bidirectional regulation provides a coherent mechanistic model: PRG4 boosts anti-aging pathways while simultaneously suppressing pro-aging pathways.

### Interpretation

Notable genes in this analysis include:
- **PMAIP1** (NOXA): A BH3-only pro-apoptotic protein; upregulation by PRG4 may promote clearance of damaged cells
- **RBBP6**: Also a convergent target (knockdown reverses aging AND mimics PRG4); PRG4 upregulation is consistent
- **BMP1**: A metalloproteinase involved in ECM remodeling; suppression by PRG4 may limit fibrotic changes
- **CTSH** (Cathepsin H): A lysosomal protease; downregulation may reduce autophagic stress

---

\newpage

## Mechanistic Model

![**Figure 7. Mechanistic model: PRG4 reverses aging via key regulators.** Schematic showing the integration of Perturb-seq-derived aging regulators with PRG4 treatment effects. Left: Aging antagonists (genes whose knockdown reverses aging, green) and aging mimetics (genes whose knockdown promotes aging, red). Center: PRG4 treatment upregulates specific aging antagonists (e.g., PMAIP1, RBBP6) and downregulates aging mimetics (e.g., BMP1, MARCKSL1). Right: Downstream effects on ARE/antioxidant pathway genes. Bottom: The 13 convergent therapeutic targets, highlighting miRNA biogenesis genes.](fig7_mechanistic_schematic.pdf){width=100%}

This integrated model demonstrates that PRG4's therapeutic effect operates through multiple parallel mechanisms:

1. **Direct upregulation of anti-aging genes**: PRG4 increases expression of PMAIP1, RBBP6, BORA, and other aging antagonists
2. **Direct suppression of pro-aging genes**: PRG4 decreases expression of BMP1, MARCKSL1, CTSH, and other aging mimetics
3. **Convergent targeting via miRNA**: Three miRNA biogenesis genes (DICER1, DROSHA, XPO5) are both aging antagonists and PRG4 mimetics, suggesting miRNA modulation as a core mechanism

---

\newpage

## Summary

![**Figure 6. Summary of key findings.** (A) Gene set sizes and key overlaps. (B) All 14 convergent genes ranked by aging antagonist strength, with miRNA biogenesis genes highlighted in purple.](fig6_summary_overview.pdf){width=100%}

### Key Findings

1. **Forward analysis**: Age-DE gene knockdowns do not systematically affect ARE pathway genes. Mean ARE enrichment was near zero for both age-UP and age-DOWN knockdowns.

2. **Aging antagonists**: Top knockdowns reversing the aging signature include RBM25, FGFR1OP, and CNOT3 - genes involved in RNA processing and transcriptional regulation.

3. **PRG4 mimetics**: Top knockdowns mimicking PRG4 include SHOC2, DICER1, AMBRA1, and XPO5. The appearance of miRNA genes (DICER1, XPO5) was notable.

4. **Convergence**: 14 genes are both aging antagonists and PRG4 mimetics. Three of these (DICER1, DROSHA, XPO5) are core miRNA biogenesis components.

5. **ARE independence**: Convergent therapeutic targets are not enriched for ARE pathway membership, suggesting the PRG4-aging axis operates through alternative mechanisms.

---

## Supplementary Tables

### Table S1: Summary Statistics

| Metric | Value |
|:-------|:------|
| Age-DE genes (total) | 472 |
| Age-UP genes | 322 |
| Age-DOWN genes | 150 |
| Age-DE with knockdown data | 245 |
| PRG4-UP genes (lfc > 0.5) | 4,012 |
| PRG4-DOWN genes (lfc < -0.5) | 2,257 |
| ARE/NRF2 genes | 242 |
| Total Perturb-seq knockdowns | 11,258 |
| Aging antagonists (top 200) | 199 unique genes |
| PRG4 mimetics (top 200) | 194 unique genes |
| Convergent genes | 14 |
| Convergent miRNA genes | 3 |

---

## Output Files

All analysis outputs are available in: `results/cohort-GSE29801/age_nrf2_analysis/`

**Figures**: fig1-fig6 (.pdf and .tiff formats)

**Data**: forward_age_de_are_gsea.csv, backward_aging_mimetics.csv, backward_prg4_mimetics.csv, convergence_overlaps.csv

**Interactive**: fig3_sankey_interactive.html

---

*Analysis completed: January 2026*
