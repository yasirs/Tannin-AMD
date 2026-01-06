## 1. Project Context (RPE Bulk RNA-seq)
Before linking to Perturb-seq, we established the transcriptomic landscape of RPE cells under Stress and PRG4 rescue.

````carousel
![PCA of RPE Bulk Samples](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/PCA_RPE_Bulk.png)
<!-- slide -->
![Volcano: H2O2 vs CTRL](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/Volcano_H2O2_vs_CTRL.png)
<!-- slide -->
![Volcano: PRG4 Rescue](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/Volcano_PRG4_Rescue.png)
````

## 2. Bridge Methodology & Technical Details

### 1.1 Signature Definition (Bulk RNA-seq)
*   **Source Data**: `RPE_gene pvals.xlsx` processed from RPE cell treatments.
*   **Discovery Threshold**: Genes were defined as differentially expressed (DEGs) if they met a nominal significance of **p < 0.05**.
*   **Signature Vectors**:
    *   **H2O2 Stress**: 5,734 DEGs (Comparison: `H2O2_vs_CTRL`).
    *   **PRG4 Baseline**: 4,618 DEGs (Comparison: `PRG4_vs_CTRL`).
    *   **PRG4 Rescue**: 5,625 DEGs (Comparison: `H2O2PRG4_vs_H2O2`).
*   **Data Values**: The `log2FoldChange` values for these significant genes were used as the reference vectors for correlation.

### 1.2 Data Harmonization
*   **Common Genes**: Bulk signatures were matched against the RPE1 Perturb-seq library using Ensembl Gene IDs.
*   **H2O2 Stress Overlay**: 2,816 common genes.
*   **PRG4 Baseline Overlay**: 2,868 common genes.
*   **PRG4 Rescue Overlay**: 3,671 common genes.
*   Only these overlapping gene sets were used for subsequent correlation calculations.

### 1.3 Correlation and Ranking
*   **Metric**: **Spearman Rank Correlation (ρ)**. This robust, non-parametric metric was chosen to account for the different data distributions between bulk bulk log2FC and single-cell Z-normalized expression.
*   **Calculation**: For each of the 2,679 knockdowns (KDs) in the Perturb-seq dataset, the Spearman correlation was calculated between the KD's expression profile and the Bulk signature vector across the common gene sets.
*   **Ranking Logic**:
    *   **Positive Correlation (High ρ)**: Indicates the KD produces a Transcriptomic state similar to the signature (e.g., AMD-mimics or PRG4-proxies).
    *   **Negative Correlation (Low ρ)**: Indicates the KD produces the opposite effect of the signature (e.g., PRG4-antagonists).

### 3.2 Cross-Signature Comparison
This Dot Plot compares how our top candidates behave across all three experimental signatures. It highlights specific genes (like **ARL4D**) that appear as strong proxies in both Baseline and Rescue conditions.

![Cross-Signature Dot Plot](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/Cross_Signature_Comparison_DotPlot.png)

## 4. Paper & PPT Presentation Suggestions

For presenting these results (e.g., in a group meeting or manuscript draft), we suggest grouping the findings as follows:

| Story Phase | Core Visual | Purpose |
| :--- | :--- | :--- |
| **1. The AMD Phenotype** | [Volcano: H2O2 vs CTRL](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/Volcano_H2O2_vs_CTRL.png) | Establish that stressed RPE cells have a massive inflammatory signature. |
| **2. The PRG4 Effect** | [Volcano: PRG4 Rescue](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/Volcano_PRG4_Rescue.png) | Demonstrate that direct addition of PRG4 significantly reverses the stress phenotype. |
| **3. Bridge Methodology** | [Bridge Evidence (SNIP1)](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/H2O2_Stress_SNIP1_bridge_scatter.png) | Show how we can "bridge" from our bulk data to identify genetic drivers in Perturb-seq. |
| **4. Top Candidates** | [Cross-Signature DotPlot](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/Cross_Signature_Comparison_DotPlot.png) | Summarize the top genetic mimics and proxies identified. |

## 5. Data Locations
*   **Code**: [correlation_analysis.py](file:///home/ysuhail/work/Tannin-AMD/code/bridge-analysis/correlation_analysis.py)
*   **Results**:
    *   [H2O2_Stress_bridge_correlations.csv](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/H2O2_Stress_bridge_correlations.csv)
    *   [PRG4_Baseline_bridge_correlations.csv](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/PRG4_Baseline_bridge_correlations.csv)
    *   [PRG4_Rescue_bridge_correlations.csv](file:///home/ysuhail/work/Tannin-AMD/results/bridge-results/PRG4_Rescue_bridge_correlations.csv)

## 4. Next Steps
*   Perform functional enrichment (GO/KEGG) on the top 50 correlates for each signature.
*   Investigate the role of DNA repair pathways (`RFC`, `RAD51C`) in PRG4-mediated protection.
