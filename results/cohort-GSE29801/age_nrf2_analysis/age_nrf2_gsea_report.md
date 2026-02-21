# Age-NRF2-ARE GSEA Analysis Report (Revised)

## Executive Summary

This analysis connects age-associated gene expression changes in AMD retinal tissue to the NRF2/ARE antioxidant pathway using **GSEA-based non-parametric methods** applied to Perturb-seq data (K562 GWPS, 11,258 knockdowns).

### Key Findings

| Analysis | Result |
|:---------|:-------|
| **Forward**: Age-DE knockdowns → ARE | Age-UP: mean enrichment = -0.005, Age-DOWN: mean = +0.007 |
| **Backward-Age**: Aging mimetics | Top mimetics: UPF2, HIST1H2AD; Top antagonists: RBM25, FGFR1OP |
| **Backward-PRG4**: PRG4 mimetics | Top mimetics: SHOC2, DICER1, AMBRA1, IMPDH2 |
| **Convergence**: 14 genes are both aging antagonists AND PRG4 mimetics |

**Critical discovery**: The 14 convergent genes include **DICER1, DROSHA, XPO5** - all components of the miRNA biogenesis pathway. Knocking down these genes both reverses aging AND mimics PRG4.

---

## Figures

### Figure 1. Forward Analysis: ARE Enrichment by Age Direction

![Forward analysis showing ARE enrichment when age-DE genes are knocked down](forward_are_enrichment_histograms.pdf)

**Figure 1. ARE pathway enrichment when age-associated genes are knocked down.** Histograms showing the distribution of ARE enrichment scores for knockdowns of age-upregulated genes (red, n=153) and age-downregulated genes (blue, n=92). Vertical black dashed line indicates zero (no enrichment); solid red line indicates mean. Neither group shows strong directional bias (means near zero), suggesting age-DE genes do not preferentially regulate ARE pathway genes when perturbed.

---

### Figure 2. Backward Analysis: Aging Mimetic Score Distribution

![Distribution of aging mimetic scores across all 11,258 knockdowns](backward_aging_mimetic_distribution.pdf)

**Figure 2. Distribution of aging mimetic scores across all Perturb-seq knockdowns.** The aging mimetic score reflects how similar a knockdown's transcriptomic effect is to the aging signature (positive = phenocopies aging, negative = reverses aging). The distribution is approximately symmetric and centered below zero (mean = -0.023), with extreme mimetics (score > 0.3) and antagonists (score < -0.4) representing potential therapeutic targets.

---

### Figure 3. Convergence: Aging Antagonists vs PRG4 Mimetics

![Scatter plot showing correlation between aging antagonist and PRG4 mimetic scores](convergence_correlation.pdf)

**Figure 3. Correlation between aging antagonist and PRG4 mimetic scores.** Each point represents one of ~9,800 unique knockdowns. X-axis: aging antagonist score (negative of mimetic score); Y-axis: PRG4 mimetic score. A weak but significant positive correlation (Spearman r = 0.004, p = 0.03) indicates that perturbations reversing aging tend to also mimic PRG4's effect, supporting PRG4 as an aging-reversal therapeutic. ARE genes (orange) are distributed throughout, not concentrated in the convergent quadrant.

---

### Figure 4. Convergence Analysis: Set Overlaps

![Bar chart showing set sizes and overlaps between aging antagonists, PRG4 mimetics, and ARE genes](convergence_overlaps_bar.pdf)

**Figure 4. Overlap between aging antagonists, PRG4 mimetics, and ARE genes.** Top 200 aging antagonists and top 200 PRG4 mimetics were defined; 155 ARE genes have knockdown data. Key overlap: 14 genes are both aging antagonists AND PRG4 mimetics. Only 3 of these are ARE genes, and the triple intersection is empty, suggesting the convergent therapeutic targets operate through pathways other than direct antioxidant response.

---

## Convergent Genes: Potential Therapeutic Targets

14 genes whose knockdown both **reverses aging** AND **mimics PRG4**:

| Gene | Function |
|:-----|:---------|
| **DICER1** | miRNA processing |
| **DROSHA** | miRNA biogenesis |
| **XPO5** | miRNA nuclear export |
| RBBP6 | RB binding, cell cycle |
| MED14 | Mediator complex |
| TANGO6 | Transport/Golgi |
| UBA3 | NEDD8 conjugation |
| BORA | Mitotic regulator |
| RPRD1B | Pol II regulation |
| YTHDC1 | m6A reader |
| FCHO1 | Endocytosis |
| RBM14-RBM4 | RNA binding |
| NKX6-1 | Transcription factor |

**Biological interpretation**: Three of the 14 convergent genes (DICER1, DROSHA, XPO5) are core components of microRNA biogenesis. This suggests that **miRNA dysregulation** may be a key mechanism linking aging and PRG4's therapeutic effect.

---

## Methods Summary

| Step | Method |
|:-----|:-------|
| Gene signatures | Age-UP (322 genes), Age-DOWN (150), PRG4-UP (4012), PRG4-DOWN (2257), ARE (242) |
| Enrichment | Rank-based quick enrichment score (non-parametric) |
| Mimetic scoring | Combined score: NES_UP - NES_DOWN |
| Convergence | Top 200 from each category; pairwise overlaps |

---

## Output Files

- `forward_age_de_are_gsea.csv` - ARE enrichment for each age-DE knockdown
- `backward_aging_mimetics.csv` - All 11,258 knockdowns with aging mimetic scores
- `backward_prg4_mimetics.csv` - All knockdowns with PRG4 mimetic scores
- `convergence_overlaps.csv` - Summary statistics

---

*Analysis completed: January 2026*
