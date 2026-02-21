# Mimetic Forward-Backward Perturb-seq Analysis Template

This document provides a detailed, reproducible template for running the forward-backward mimetic analysis using Perturb-seq data. It can be adapted to different gene sets (e.g., disease-associated genes from a different cohort) and different Perturb-seq datasets (e.g., K562 vs RPE1).

---

## Overview

### Analysis Goals

1. **Forward Analysis**: Test whether knockdown of condition-associated genes (e.g., age-DE, disease-DE) affects a target pathway (e.g., ARE/antioxidant response)
2. **Backward Analysis - Mimetics**: Identify knockdowns that phenocopy the condition signature
3. **Backward Analysis - Antagonists**: Identify knockdowns that reverse the condition signature
4. **PRG4 Mimetics**: Identify knockdowns that recapitulate PRG4's therapeutic effect
5. **Convergence**: Find genes at the intersection of antagonists and PRG4 mimetics
6. **Direct Targets**: Identify condition-relevant genes directly regulated by PRG4

### Key Outputs

| Output | Description |
|:-------|:------------|
| `forward_{{ANALYSIS_NAME}}_pathway_gsea.csv` | Forward analysis results |
| `backward_{{ANALYSIS_NAME}}_mimetics.csv` | All knockdowns with mimetic scores |
| `backward_prg4_mimetics.csv` | All knockdowns with PRG4 mimetic scores |
| `convergence_overlaps.csv` | Summary of set overlaps |
| `{{ANALYSIS_NAME}}_mediator_table.csv` | Top pathway mediators |
| `prg4_direct_targets.csv` | PRG4 targets among mimetics/antagonists |
| Figures 1-7 | Publication-ready figures (PDF + TIFF) |
| Comprehensive report | Markdown + PDF + DOCX |

---

## Parameters

### Required Inputs

```yaml
# ==========================================
# ANALYSIS PARAMETERS - MODIFY THESE
# ==========================================

ANALYSIS_NAME: "age_nrf2"  # Used in output file naming

# Condition Gene Set (the signature you want to test)
CONDITION_UP_GENES:
  path: "path/to/condition_up_genes.csv"
  format: "CSV with column 'gene' containing HGNC symbols"
  description: "Genes upregulated in condition (e.g., aging, disease)"
  
CONDITION_DOWN_GENES:
  path: "path/to/condition_down_genes.csv"
  format: "CSV with column 'gene' containing HGNC symbols"
  description: "Genes downregulated in condition"

# Target Pathway Gene Set (optional, for forward analysis)
TARGET_PATHWAY:
  path: "path/to/pathway_genes.csv"
  format: "CSV with column 'gene' containing HGNC symbols"
  name: "ARE"  # For labeling

# Perturb-seq Dataset
PERTURBSEQ_DATA:
  path: "path/to/perturbseq.h5ad"
  format: "AnnData h5ad file with normalized expression"
  obs_column: "target_gene"  # Column containing knockdown target
  description: "K562 GWPS or RPE1 Perturb-seq"

# PRG4 Signature (FIXED - do not modify)
PRG4_DATA:
  path: "data/RPE_DE_results.csv"
  up_column: "H2O2PRG4_vs_H2O2.log2FoldChange"
  pval_column: "H2O2PRG4_vs_H2O2.pvalue"
  gene_column: "hgnc_symbol"
  lfc_threshold: 0.5  # |LFC| > 0.5 for PRG4-UP/DOWN

# Output Directory
OUTPUT_DIR: "results/{{COHORT_NAME}}/{{ANALYSIS_NAME}}_analysis"
```

### Input File Formats

#### Condition Gene Set CSV
```csv
gene,direction,log2fc,pvalue
TP53,UP,1.23,0.001
BRCA1,DOWN,-0.87,0.0001
...
```
- Required columns: `gene`
- Optional columns: `direction`, `log2fc`, `pvalue` (for ranking/filtering)

#### Perturb-seq h5ad File
```python
# Expected structure:
adata.obs  # DataFrame with perturbation metadata
  - index: cell/sample IDs
  - 'target_gene': HGNC symbol of knocked-down gene
  - Other QC columns (optional)

adata.var  # DataFrame with gene metadata
  - index: gene names (HGNC symbols)
  
adata.X    # Expression matrix (cells × genes)
  # Should be normalized (e.g., log-transformed CPM)
```

---

## Step-by-Step Implementation

### Step 0: Setup and Data Loading

```python
#!/usr/bin/env python3
"""
Mimetic Forward-Backward Analysis
Template implementation - modify parameters above
"""

import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from PIL import Image

# ==========================================
# LOAD PARAMETERS (modify these)
# ==========================================
ANALYSIS_NAME = "age_nrf2"
OUTPUT_DIR = Path(f"results/cohort-GSE29801/{ANALYSIS_NAME}_analysis")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load gene sets
condition_up = set(pd.read_csv("path/to/condition_up.csv")['gene'].dropna())
condition_down = set(pd.read_csv("path/to/condition_down.csv")['gene'].dropna())
target_pathway = set(pd.read_csv("path/to/pathway.csv")['gene'].dropna())

print(f"Condition-UP genes: {len(condition_up)}")
print(f"Condition-DOWN genes: {len(condition_down)}")
print(f"Target pathway genes: {len(target_pathway)}")

# Load Perturb-seq data
adata = ad.read_h5ad("path/to/perturbseq.h5ad")
perturbseq_genes = set(adata.var_names)
print(f"Perturb-seq knockdowns: {adata.n_obs}")
print(f"Perturb-seq genes profiled: {adata.n_vars}")

# Load PRG4 data (FIXED)
prg4_full = pd.read_csv("data/RPE_DE_results.csv")
prg4_full = prg4_full[prg4_full['hgnc_symbol'].notna()]
prg4_up = set(prg4_full[prg4_full['H2O2PRG4_vs_H2O2.log2FoldChange'] > 0.5]['hgnc_symbol'])
prg4_down = set(prg4_full[prg4_full['H2O2PRG4_vs_H2O2.log2FoldChange'] < -0.5]['hgnc_symbol'])
print(f"PRG4-UP: {len(prg4_up)}, PRG4-DOWN: {len(prg4_down)}")
```

---

### Step 1: Forward Analysis

**Goal**: For each condition-associated gene that has a knockdown, compute enrichment of the target pathway.

```python
def quick_enrichment_score(ranked_genes, gene_set):
    """
    Non-parametric enrichment score based on mean rank difference.
    
    Args:
        ranked_genes: pd.Series with gene names as index, expression values as values
        gene_set: set of gene names to test for enrichment
    
    Returns:
        enrichment_score: float (-1 to +1), positive = genes enriched at top
        n_overlap: int, number of gene_set genes found in ranked list
    """
    gene_set_in_list = [g for g in gene_set if g in ranked_genes.index]
    if len(gene_set_in_list) < 3:
        return np.nan, len(gene_set_in_list)
    
    n_genes = len(ranked_genes)
    gene_ranks = ranked_genes.rank(ascending=False)  # 1 = highest expression
    
    set_ranks = gene_ranks[gene_set_in_list].values
    expected_mean = (n_genes + 1) / 2
    observed_mean = np.mean(set_ranks)
    
    # Normalize to [-1, 1]
    enrichment = (expected_mean - observed_mean) / expected_mean
    return enrichment, len(gene_set_in_list)

# Get condition genes with knockdown profiles
condition_genes = condition_up | condition_down
condition_with_kd = condition_genes & set(adata.obs['target_gene'].unique())
print(f"Condition genes with knockdowns: {len(condition_with_kd)}")

# Run forward analysis
forward_results = []
for gene in condition_with_kd:
    # Get knockdown profile (average across all perturbations of this gene)
    mask = adata.obs['target_gene'] == gene
    if mask.sum() == 0:
        continue
    
    profile = pd.Series(
        adata[mask].X.mean(axis=0).A1 if hasattr(adata.X, 'A1') else adata[mask].X.mean(axis=0),
        index=adata.var_names
    )
    
    # Compute pathway enrichment
    enrichment, n_genes = quick_enrichment_score(profile, target_pathway)
    
    forward_results.append({
        'target_gene': gene,
        'direction': 'UP' if gene in condition_up else 'DOWN',
        'pathway_enrichment': enrichment,
        'n_pathway_genes': n_genes
    })

forward_df = pd.DataFrame(forward_results)
forward_df.to_csv(OUTPUT_DIR / f"forward_{ANALYSIS_NAME}_pathway_gsea.csv", index=False)
print(f"Forward analysis: {len(forward_df)} knockdowns tested")
```

---

### Step 2: Backward Analysis - Condition Mimetics

**Goal**: For each knockdown, compute how similar its effect is to the condition signature.

```python
# Filter data
adata_clean = adata[adata.obs['target_gene'] != 'non-targeting'].copy()

# Deduplicate to one entry per gene (take first/best)
gene_to_idx = {}
for idx, row in adata_clean.obs.iterrows():
    gene = row['target_gene']
    if gene not in gene_to_idx:
        gene_to_idx[gene] = idx

unique_genes = list(gene_to_idx.keys())
unique_idx = list(gene_to_idx.values())
print(f"Unique knockdown genes: {len(unique_genes)}")

# Compute mimetic scores
mimetic_results = []
for gene in unique_genes:
    mask = adata_clean.obs['target_gene'] == gene
    profile = pd.Series(
        adata_clean[mask].X.mean(axis=0).A1 if hasattr(adata_clean.X, 'A1') else adata_clean[mask].X.mean(axis=0),
        index=adata_clean.var_names
    )
    
    # Enrichment of condition-UP genes
    enrich_up, n_up = quick_enrichment_score(profile, condition_up)
    
    # Enrichment of condition-DOWN genes  
    enrich_down, n_down = quick_enrichment_score(profile, condition_down)
    
    # Mimetic score: positive = phenocopies condition
    # Knockdown causes UP genes to go up AND DOWN genes to go down
    mimetic_score = (enrich_up if not np.isnan(enrich_up) else 0) - \
                    (enrich_down if not np.isnan(enrich_down) else 0)
    
    mimetic_results.append({
        'target_gene': gene,
        'enrichment_UP': enrich_up,
        'enrichment_DOWN': enrich_down,
        'mimetic_score': mimetic_score,
        'n_up_genes': n_up,
        'n_down_genes': n_down,
        'is_pathway_gene': gene in target_pathway
    })

mimetic_df = pd.DataFrame(mimetic_results)
mimetic_df = mimetic_df.sort_values('mimetic_score').reset_index(drop=True)
mimetic_df.to_csv(OUTPUT_DIR / f"backward_{ANALYSIS_NAME}_mimetics.csv", index=False)

# Define top mimetics and antagonists
top_mimetics = set(mimetic_df.tail(200)['target_gene'])  # High score = phenocopies condition
top_antagonists = set(mimetic_df.head(200)['target_gene'])  # Low score = reverses condition

print(f"Top 200 mimetics (phenocopy condition)")
print(f"Top 200 antagonists (reverse condition)")
```

---

### Step 3: Backward Analysis - PRG4 Mimetics

**Goal**: For each knockdown, compute how similar its effect is to PRG4 treatment.

```python
prg4_mimetic_results = []
for gene in unique_genes:
    mask = adata_clean.obs['target_gene'] == gene
    profile = pd.Series(
        adata_clean[mask].X.mean(axis=0).A1 if hasattr(adata_clean.X, 'A1') else adata_clean[mask].X.mean(axis=0),
        index=adata_clean.var_names
    )
    
    # Enrichment of PRG4-UP genes
    enrich_prg4_up, n_up = quick_enrichment_score(profile, prg4_up)
    
    # Enrichment of PRG4-DOWN genes
    enrich_prg4_down, n_down = quick_enrichment_score(profile, prg4_down)
    
    # PRG4 mimetic score
    prg4_mimetic_score = (enrich_prg4_up if not np.isnan(enrich_prg4_up) else 0) - \
                          (enrich_prg4_down if not np.isnan(enrich_prg4_down) else 0)
    
    prg4_mimetic_results.append({
        'target_gene': gene,
        'enrichment_PRG4_UP': enrich_prg4_up,
        'enrichment_PRG4_DOWN': enrich_prg4_down,
        'prg4_mimetic_score': prg4_mimetic_score
    })

prg4_df = pd.DataFrame(prg4_mimetic_results)
prg4_df = prg4_df.sort_values('prg4_mimetic_score', ascending=False).reset_index(drop=True)
prg4_df.to_csv(OUTPUT_DIR / "backward_prg4_mimetics.csv", index=False)

top_prg4_mimetics = set(prg4_df.head(200)['target_gene'])
print(f"Top 200 PRG4 mimetics")
```

---

### Step 4: Convergence Analysis

**Goal**: Find intersection of condition antagonists and PRG4 mimetics.

```python
from scipy import stats

# Compute overlaps
convergent_genes = top_antagonists & top_prg4_mimetics
pathway_in_antagonists = top_antagonists & target_pathway
pathway_in_prg4_mimetics = top_prg4_mimetics & target_pathway

print(f"Convergent genes (antagonist ∩ PRG4 mimetic): {len(convergent_genes)}")
print(f"Convergent genes: {sorted(convergent_genes)}")

# Merge for correlation
merged = mimetic_df.merge(prg4_df[['target_gene', 'prg4_mimetic_score']], on='target_gene')
merged = merged.dropna(subset=['mimetic_score', 'prg4_mimetic_score'])

# Spearman correlation
r, p = stats.spearmanr(-merged['mimetic_score'], merged['prg4_mimetic_score'])
print(f"Correlation (antagonist vs PRG4 mimetic): r={r:.4f}, p={p:.4e}")

# Save convergence results
convergence = {
    'n_antagonists': len(top_antagonists),
    'n_prg4_mimetics': len(top_prg4_mimetics),
    'n_convergent': len(convergent_genes),
    'convergent_genes': ','.join(sorted(convergent_genes)),
    'spearman_r': r,
    'spearman_p': p
}
pd.DataFrame([convergence]).to_csv(OUTPUT_DIR / "convergence_overlaps.csv", index=False)
```

---

### Step 5: Forward Mediator Table

**Goal**: Identify top pathway activators/suppressors from each direction.

```python
mediator_table = []

for direction in ['UP', 'DOWN']:
    subset = forward_df[forward_df['direction'] == direction]
    
    # Top activators (positive enrichment)
    activators = subset.nlargest(6, 'pathway_enrichment')
    for _, row in activators.iterrows():
        mediator_table.append({
            'category': f'Condition-{direction}',
            'gene': row['target_gene'],
            'effect': 'Activates',
            'enrichment': row['pathway_enrichment']
        })
    
    # Top suppressors (negative enrichment)
    suppressors = subset.nsmallest(6, 'pathway_enrichment')
    for _, row in suppressors.iterrows():
        mediator_table.append({
            'category': f'Condition-{direction}',
            'gene': row['target_gene'],
            'effect': 'Suppresses',
            'enrichment': row['pathway_enrichment']
        })

mediator_df = pd.DataFrame(mediator_table)
mediator_df.to_csv(OUTPUT_DIR / f"forward_{ANALYSIS_NAME}_mediator_table.csv", index=False)
```

---

### Step 6: PRG4 Direct Targets

**Goal**: Identify condition-relevant genes directly regulated by PRG4.

```python
# PRG4 full data with p-values
prg4_ranked = prg4_full[['hgnc_symbol', 'H2O2PRG4_vs_H2O2.log2FoldChange', 'H2O2PRG4_vs_H2O2.pvalue']].copy()
prg4_ranked.columns = ['gene', 'prg4_lfc', 'prg4_pval']
prg4_ranked = prg4_ranked.dropna()
prg4_ranked['rank_score'] = prg4_ranked['prg4_lfc'].abs() * (-np.log10(prg4_ranked['prg4_pval'].clip(1e-300)))

# PRG4-UP among antagonists
prg4_up_antagonists = prg4_ranked[(prg4_ranked['prg4_lfc'] > 0) & 
                                   (prg4_ranked['gene'].isin(top_antagonists))]
prg4_up_antagonists = prg4_up_antagonists.nlargest(10, 'rank_score')

# PRG4-DOWN among mimetics
prg4_down_mimetics = prg4_ranked[(prg4_ranked['prg4_lfc'] < 0) & 
                                  (prg4_ranked['gene'].isin(top_mimetics))]
prg4_down_mimetics = prg4_down_mimetics.nlargest(10, 'rank_score')

# Combine
prg4_targets = []
for _, row in prg4_up_antagonists.iterrows():
    prg4_targets.append({
        'gene': row['gene'],
        'prg4_direction': 'UP',
        'prg4_lfc': row['prg4_lfc'],
        'prg4_pval': row['prg4_pval'],
        'condition_role': 'Antagonist',
        'interpretation': 'PRG4 boosts anti-condition pathway'
    })

for _, row in prg4_down_mimetics.iterrows():
    prg4_targets.append({
        'gene': row['gene'],
        'prg4_direction': 'DOWN',
        'prg4_lfc': row['prg4_lfc'],
        'prg4_pval': row['prg4_pval'],
        'condition_role': 'Mimetic',
        'interpretation': 'PRG4 suppresses pro-condition pathway'
    })

prg4_targets_df = pd.DataFrame(prg4_targets)
prg4_targets_df.to_csv(OUTPUT_DIR / "prg4_direct_targets.csv", index=False)
```

---

## Figure Generation

### Figure 1: Analysis Overview Schematic

Create a flow diagram showing:
- Left: Input gene sets (Condition-UP, Condition-DOWN, Target Pathway)
- Center: Perturb-seq analysis hub
- Right: Output categories (Mimetics, Antagonists, PRG4 Mimetics)
- Bottom: Convergent targets

### Figure 2: Forward Analysis

Three-panel figure:
- A: Histograms of pathway enrichment by direction (Condition-UP vs Condition-DOWN)
- B: Bar chart of top activators/suppressors
- C: Summary statistics (mean ± SD per direction)

### Figure 3: Flow Diagram

Schematic showing knockdown filtering:
- All knockdowns → Mimetics (top 200) / Neutral / Antagonists (top 200)
- Antagonists → PRG4 mimetic overlap → Convergent genes

### Figure 4: Convergence Scatter

- X-axis: Antagonist score (negative of mimetic score)
- Y-axis: PRG4 mimetic score
- Gray: Background (subsample to 1000 points)
- Red squares: Target pathway genes
- Blue diamonds: Convergent genes
- Purple stars: Key genes of interest (label with arrows)

**IMPORTANT**: Subsample background points to ~1000 to reduce file size.

### Figure 5: Pathway Schematic

If convergent genes map to a known pathway (e.g., miRNA biogenesis), create a schematic showing:
- Pathway steps with gene names
- Highlighting which are convergent

### Figure 6: Summary Overview

Two-panel bar chart:
- A: Set sizes and overlaps
- B: Convergent genes ranked by antagonist score

### Figure 7: Mechanistic Model

Comprehensive schematic with gene names:
- Left boxes: Top antagonists and mimetics (6 genes each)
- Center: PRG4 treatment effects (which genes it regulates)
- Right: Target pathway genes
- Arrows showing regulatory relationships

### Figure Saving Helper

```python
def save_figure(fig, name, output_dir, dpi=300):
    """Save figure in PDF and LZW-compressed TIFF"""
    from PIL import Image
    
    pdf_path = output_dir / f"{name}.pdf"
    tiff_path = output_dir / f"{name}.tiff"
    
    fig.savefig(pdf_path, dpi=dpi, bbox_inches='tight')
    
    # Save as PNG then convert to TIFF with LZW compression
    png_temp = output_dir / f"{name}_temp.png"
    fig.savefig(png_temp, dpi=dpi, bbox_inches='tight')
    
    img = Image.open(png_temp)
    img.save(tiff_path, compression='tiff_lzw')
    png_temp.unlink()
    
    print(f"Saved: {name}.pdf, {name}.tiff")
```

---

## Report Generation

### Markdown Template

The report should follow this structure:

```markdown
---
title: "{{ANALYSIS_NAME}} Regulatory Connection Analysis via Perturb-seq"
author: "Tannin-AMD Project"
date: "{{DATE}}"
geometry: margin=1in
fontsize: 11pt
---

# {{ANALYSIS_NAME}} Regulatory Connection Analysis

## Abstract
[2-3 sentences summarizing the analysis and key findings]

## Methods
[Data sources, analytical framework, enrichment scoring method]

## Results

### Forward Analysis
[Rationale → Figure 2 → Results → Table 1 (mediators) → Interpretation]

### Backward Analysis: Condition Mimetics and Antagonists
[Rationale → Figure 3 → Top genes → Interpretation]

### Backward Analysis: PRG4 Mimetics
[Rationale → Top genes → Interpretation]

### Convergence Analysis
[Rationale → Figure 4 → Convergent gene table → Correlation stats → Figure 5 (if applicable)]

### PRG4 Direct Targets
[Rationale → Table 2 → Interpretation]

### Mechanistic Model
[Figure 7 → Interpretation]

## Summary
[Figure 6 → 5 bullet points of key findings]

## Supplementary Tables
[Summary statistics table]
```

### Compilation

```bash
cd {{OUTPUT_DIR}}
pandoc {{ANALYSIS_NAME}}_comprehensive_report.md -o {{ANALYSIS_NAME}}_comprehensive_report.pdf --pdf-engine=pdflatex
pandoc {{ANALYSIS_NAME}}_comprehensive_report.md -o {{ANALYSIS_NAME}}_comprehensive_report.docx
```

**Note**: Replace Unicode characters (∩, ✓, ≈) with ASCII equivalents before PDF compilation.

---

## Validation Checklist

Before considering the analysis complete, verify:

- [ ] Forward analysis CSV has correct number of genes
- [ ] Mimetic score distribution is approximately symmetric around 0
- [ ] Non-targeting controls are excluded from mimetic/antagonist ranking
- [ ] Only one representative per gene (deduplicated)
- [ ] Convergent genes list does not include "non-targeting"
- [ ] All figures render correctly in PDF
- [ ] Figure file sizes are reasonable (<500 KB for scatter plots)
- [ ] Report compiles without LaTeX errors

---

## Example Invocations

### Example 1: Age-DE genes (GSE29801) with K562 Perturb-seq

```python
ANALYSIS_NAME = "age_nrf2"
condition_up = "results/cohort-GSE29801/age_nrf2_analysis/age_de_extramacular.csv"  # Filter to direction=UP
condition_down = "results/cohort-GSE29801/age_nrf2_analysis/age_de_extramacular.csv"  # Filter to direction=DOWN
target_pathway = "results/cohort-GSE29801/age_nrf2_analysis/ARE_gene_set.csv"
perturbseq = "data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad"
```

### Example 2: AMD-DE genes (GSE135092) with RPE1 Perturb-seq

```python
ANALYSIS_NAME = "amd_nrf2"
condition_up = "results/cohort-GSE135092/amd_de_genes.csv"  # Filter to direction=UP
condition_down = "results/cohort-GSE135092/amd_de_genes.csv"  # Filter to direction=DOWN
target_pathway = "results/cohort-GSE29801/age_nrf2_analysis/ARE_gene_set.csv"  # Same pathway
perturbseq = "data/external/perturbseq/RPE1_perturbseq.h5ad"  # Different cell type
```

---

## Troubleshooting

| Issue | Solution |
|:------|:---------|
| Large figure file size | Subsample background points; use `rasterized=True` for scatter |
| Unicode LaTeX errors | Replace ∩ with "AND", ✓ with "[check]", ≈ with "~" |
| Zero convergent genes | Lower threshold (top 300-500 instead of 200) |
| Empty target pathway | Check gene name format matches Perturb-seq var_names |
| Non-targeting in results | Filter `adata.obs['target_gene'] != 'non-targeting'` |

---

*Template version: 1.0 - January 2026*
