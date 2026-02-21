#!/usr/bin/env python3
"""
Age-NRF2-ARE Analysis using RPE1 Perturb-seq Data
Same methodology as K562 analysis, separate output files
"""

import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import FancyBboxPatch
from scipy import stats
from PIL import Image

# ==========================================
# PATHS - RPE1 VERSION
# ==========================================
BASE_DIR = Path("/home/ysuhail/work/Tannin-AMD")
K562_RESULTS_DIR = BASE_DIR / "results/cohort-GSE29801/age_nrf2_analysis"
OUTPUT_DIR = BASE_DIR / "results/cohort-GSE29801/age_nrf2_analysis_RPE1"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PERTURBSEQ_PATH = BASE_DIR / "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad"

print("=" * 60)
print("AGE-NRF2-ARE ANALYSIS - RPE1 PERTURB-SEQ VERSION")
print("=" * 60)

# ==========================================
# STEP 0: LOAD DATA
# ==========================================
print("\n[Step 0] Loading data...")

# Load gene sets from K562 analysis (same gene sets, different Perturb-seq)
age_de = pd.read_csv(K562_RESULTS_DIR / "age_de_extramacular.csv")
age_up = set(age_de[age_de['logFC'] > 0]['gene_name'].dropna())
age_down = set(age_de[age_de['logFC'] < 0]['gene_name'].dropna())
are_genes = set(pd.read_csv(K562_RESULTS_DIR / "ARE_gene_set.csv")['gene'].dropna())

print(f"  Age-UP genes: {len(age_up)}")
print(f"  Age-DOWN genes: {len(age_down)}")
print(f"  ARE genes: {len(are_genes)}")

# Load PRG4 data (FIXED)
prg4_full = pd.read_csv(BASE_DIR / "data/RPE_DE_results.csv")
prg4_full = prg4_full[prg4_full['hgnc_symbol'].notna() & (prg4_full['hgnc_symbol'] != '')]
prg4_up = set(prg4_full[prg4_full['H2O2PRG4_vs_H2O2.log2FoldChange'] > 0.5]['hgnc_symbol'])
prg4_down = set(prg4_full[prg4_full['H2O2PRG4_vs_H2O2.log2FoldChange'] < -0.5]['hgnc_symbol'])
print(f"  PRG4-UP: {len(prg4_up)}, PRG4-DOWN: {len(prg4_down)}")

# Load RPE1 Perturb-seq
print(f"\n  Loading RPE1 Perturb-seq from: {PERTURBSEQ_PATH}")
adata = ad.read_h5ad(PERTURBSEQ_PATH)
print(f"  RPE1 knockdowns: {adata.n_obs}")
print(f"  RPE1 genes profiled: {adata.n_vars}")
print(f"  obs columns: {adata.obs.columns.tolist()}")

# Identify target gene column
target_col = None
for col in ['target_gene', 'gene', 'perturbation', 'guide_identity']:
    if col in adata.obs.columns:
        target_col = col
        break

if target_col is None:
    # Try to infer from index or first column
    print("  Warning: No target_gene column found, checking alternatives...")
    print(f"  Sample obs head:\n{adata.obs.head()}")

# ==========================================
# HELPER FUNCTIONS
# ==========================================
def quick_enrichment_score(ranked_genes, gene_set):
    """Non-parametric enrichment score"""
    gene_set_in_list = [g for g in gene_set if g in ranked_genes.index]
    if len(gene_set_in_list) < 3:
        return np.nan, len(gene_set_in_list)
    
    n_genes = len(ranked_genes)
    gene_ranks = ranked_genes.rank(ascending=False)
    set_ranks = gene_ranks[gene_set_in_list].values
    expected_mean = (n_genes + 1) / 2
    observed_mean = np.mean(set_ranks)
    enrichment = (expected_mean - observed_mean) / expected_mean
    return enrichment, len(gene_set_in_list)

def save_figure(fig, name, dpi=300):
    """Save figure in PDF and LZW-compressed TIFF"""
    pdf_path = OUTPUT_DIR / f"{name}.pdf"
    fig.savefig(pdf_path, dpi=dpi, bbox_inches='tight')
    
    png_temp = OUTPUT_DIR / f"{name}_temp.png"
    fig.savefig(png_temp, dpi=dpi, bbox_inches='tight')
    img = Image.open(png_temp)
    img.save(OUTPUT_DIR / f"{name}.tiff", compression='tiff_lzw')
    png_temp.unlink()
    print(f"  Saved: {name}.pdf, {name}.tiff")

# ==========================================
# STEP 1: FORWARD ANALYSIS
# ==========================================
print("\n[Step 1] Forward Analysis: Age-DE knockdowns → ARE effect...")

if target_col:
    # Filter to age-DE genes with knockdowns
    all_targets = set(adata.obs[target_col].unique())
    age_de_genes = age_up | age_down
    age_with_kd = age_de_genes & all_targets
    print(f"  Age-DE genes with RPE1 knockdowns: {len(age_with_kd)}")
    
    forward_results = []
    for gene in age_with_kd:
        mask = adata.obs[target_col] == gene
        if mask.sum() == 0:
            continue
        
        # Get expression profile
        profile_values = adata[mask].X
        if hasattr(profile_values, 'toarray'):
            profile_values = profile_values.toarray()
        profile = pd.Series(profile_values.mean(axis=0), index=adata.var_names)
        
        enrichment, n_genes = quick_enrichment_score(profile, are_genes)
        
        forward_results.append({
            'target_gene': gene,
            'age_direction': 'UP' if gene in age_up else 'DOWN',
            'ARE_enrichment': enrichment,
            'n_ARE_genes': n_genes
        })
    
    forward_df = pd.DataFrame(forward_results)
    forward_df.to_csv(OUTPUT_DIR / "forward_age_de_are_gsea.csv", index=False)
    print(f"  Saved: forward_age_de_are_gsea.csv ({len(forward_df)} genes)")
    
    # Summary stats
    if len(forward_df) > 0:
        age_up_mean = forward_df[forward_df['age_direction'] == 'UP']['ARE_enrichment'].mean()
        age_down_mean = forward_df[forward_df['age_direction'] == 'DOWN']['ARE_enrichment'].mean()
        print(f"  Age-UP mean ARE enrichment: {age_up_mean:.4f}")
        print(f"  Age-DOWN mean ARE enrichment: {age_down_mean:.4f}")
else:
    print("  ERROR: Cannot identify target gene column in RPE1 data")
    forward_df = pd.DataFrame()

# ==========================================
# STEP 2: BACKWARD ANALYSIS - AGING MIMETICS
# ==========================================
print("\n[Step 2] Backward Analysis: Aging mimetics/antagonists...")

if target_col:
    # Remove non-targeting and deduplicate
    adata_clean = adata[~adata.obs[target_col].str.contains('non-targeting|control', case=False, na=False)].copy()
    unique_genes = adata_clean.obs[target_col].unique()
    print(f"  Unique knockdown genes (after filtering): {len(unique_genes)}")
    
    mimetic_results = []
    for gene in unique_genes:
        mask = adata_clean.obs[target_col] == gene
        if mask.sum() == 0:
            continue
        
        profile_values = adata_clean[mask].X
        if hasattr(profile_values, 'toarray'):
            profile_values = profile_values.toarray()
        profile = pd.Series(profile_values.mean(axis=0), index=adata_clean.var_names)
        
        enrich_up, n_up = quick_enrichment_score(profile, age_up)
        enrich_down, n_down = quick_enrichment_score(profile, age_down)
        
        mimetic_score = (enrich_up if not np.isnan(enrich_up) else 0) - \
                        (enrich_down if not np.isnan(enrich_down) else 0)
        
        mimetic_results.append({
            'target_gene': gene,
            'enrichment_age_UP': enrich_up,
            'enrichment_age_DOWN': enrich_down,
            'mimetic_score': mimetic_score,
            'n_up_genes': n_up,
            'n_down_genes': n_down,
            'is_ARE_gene': gene in are_genes
        })
    
    aging_mim_df = pd.DataFrame(mimetic_results)
    aging_mim_df = aging_mim_df.sort_values('mimetic_score').reset_index(drop=True)
    aging_mim_df.to_csv(OUTPUT_DIR / "backward_aging_mimetics.csv", index=False)
    print(f"  Saved: backward_aging_mimetics.csv ({len(aging_mim_df)} genes)")
    
    # Define top mimetics and antagonists
    top_mimetics = set(aging_mim_df.tail(200)['target_gene'])
    top_antagonists = set(aging_mim_df.head(200)['target_gene'])
    print(f"  Top 200 mimetics (phenocopy aging)")
    print(f"  Top 200 antagonists (reverse aging)")
else:
    aging_mim_df = pd.DataFrame()
    top_mimetics = set()
    top_antagonists = set()

# ==========================================
# STEP 3: BACKWARD ANALYSIS - PRG4 MIMETICS
# ==========================================
print("\n[Step 3] Backward Analysis: PRG4 mimetics...")

if target_col and len(unique_genes) > 0:
    prg4_mimetic_results = []
    for gene in unique_genes:
        mask = adata_clean.obs[target_col] == gene
        if mask.sum() == 0:
            continue
        
        profile_values = adata_clean[mask].X
        if hasattr(profile_values, 'toarray'):
            profile_values = profile_values.toarray()
        profile = pd.Series(profile_values.mean(axis=0), index=adata_clean.var_names)
        
        enrich_prg4_up, n_up = quick_enrichment_score(profile, prg4_up)
        enrich_prg4_down, n_down = quick_enrichment_score(profile, prg4_down)
        
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
    print(f"  Saved: backward_prg4_mimetics.csv ({len(prg4_df)} genes)")
    
    top_prg4_mimetics = set(prg4_df.head(200)['target_gene'])
else:
    prg4_df = pd.DataFrame()
    top_prg4_mimetics = set()

# ==========================================
# STEP 4: CONVERGENCE ANALYSIS
# ==========================================
print("\n[Step 4] Convergence Analysis...")

if len(aging_mim_df) > 0 and len(prg4_df) > 0:
    convergent_genes = top_antagonists & top_prg4_mimetics
    print(f"  Convergent genes (antagonist AND PRG4 mimetic): {len(convergent_genes)}")
    print(f"  Genes: {sorted(convergent_genes)}")
    
    # Correlation
    merged = aging_mim_df.merge(prg4_df[['target_gene', 'prg4_mimetic_score']], on='target_gene')
    merged = merged.dropna(subset=['mimetic_score', 'prg4_mimetic_score'])
    
    if len(merged) > 10:
        r, p = stats.spearmanr(-merged['mimetic_score'], merged['prg4_mimetic_score'])
        print(f"  Spearman correlation: r={r:.4f}, p={p:.4e}")
    else:
        r, p = np.nan, np.nan
    
    # Save convergence
    convergence = {
        'n_antagonists': len(top_antagonists),
        'n_prg4_mimetics': len(top_prg4_mimetics),
        'n_convergent': len(convergent_genes),
        'convergent_genes': ','.join(sorted(convergent_genes)),
        'spearman_r': r,
        'spearman_p': p
    }
    pd.DataFrame([convergence]).to_csv(OUTPUT_DIR / "convergence_overlaps.csv", index=False)
else:
    convergent_genes = set()

# ==========================================
# STEP 5: MEDIATOR TABLE
# ==========================================
print("\n[Step 5] Creating mediator table...")

if len(forward_df) > 0:
    mediator_table = []
    for direction in ['UP', 'DOWN']:
        subset = forward_df[forward_df['age_direction'] == direction]
        if len(subset) > 0:
            # Activators
            activators = subset.nlargest(min(6, len(subset)), 'ARE_enrichment')
            for _, row in activators.iterrows():
                mediator_table.append({
                    'category': f'Age-{direction}',
                    'gene': row['target_gene'],
                    'effect': 'Activates',
                    'enrichment': row['ARE_enrichment']
                })
            # Suppressors
            suppressors = subset.nsmallest(min(6, len(subset)), 'ARE_enrichment')
            for _, row in suppressors.iterrows():
                mediator_table.append({
                    'category': f'Age-{direction}',
                    'gene': row['target_gene'],
                    'effect': 'Suppresses',
                    'enrichment': row['ARE_enrichment']
                })
    
    mediator_df = pd.DataFrame(mediator_table)
    mediator_df.to_csv(OUTPUT_DIR / "forward_are_mediator_table.csv", index=False)
    print(f"  Saved: forward_are_mediator_table.csv ({len(mediator_df)} entries)")

# ==========================================
# STEP 6: PRG4 DIRECT TARGETS
# ==========================================
print("\n[Step 6] PRG4 Direct Targets...")

if len(aging_mim_df) > 0:
    prg4_ranked = prg4_full[['hgnc_symbol', 'H2O2PRG4_vs_H2O2.log2FoldChange', 'H2O2PRG4_vs_H2O2.pvalue']].copy()
    prg4_ranked.columns = ['gene', 'prg4_lfc', 'prg4_pval']
    prg4_ranked = prg4_ranked.dropna()
    prg4_ranked['rank_score'] = prg4_ranked['prg4_lfc'].abs() * (-np.log10(prg4_ranked['prg4_pval'].clip(1e-300)))
    
    # PRG4-UP among antagonists
    prg4_up_antag = prg4_ranked[(prg4_ranked['prg4_lfc'] > 0) & (prg4_ranked['gene'].isin(top_antagonists))]
    prg4_up_antag = prg4_up_antag.nlargest(10, 'rank_score')
    
    # PRG4-DOWN among mimetics
    prg4_down_mim = prg4_ranked[(prg4_ranked['prg4_lfc'] < 0) & (prg4_ranked['gene'].isin(top_mimetics))]
    prg4_down_mim = prg4_down_mim.nlargest(10, 'rank_score')
    
    prg4_targets = []
    for _, row in prg4_up_antag.iterrows():
        prg4_targets.append({
            'gene': row['gene'],
            'prg4_direction': 'UP',
            'prg4_lfc': row['prg4_lfc'],
            'prg4_pval': row['prg4_pval'],
            'condition_role': 'Antagonist'
        })
    for _, row in prg4_down_mim.iterrows():
        prg4_targets.append({
            'gene': row['gene'],
            'prg4_direction': 'DOWN',
            'prg4_lfc': row['prg4_lfc'],
            'prg4_pval': row['prg4_pval'],
            'condition_role': 'Mimetic'
        })
    
    prg4_targets_df = pd.DataFrame(prg4_targets)
    prg4_targets_df.to_csv(OUTPUT_DIR / "prg4_direct_targets.csv", index=False)
    print(f"  Saved: prg4_direct_targets.csv ({len(prg4_targets_df)} entries)")

# ==========================================
# FIGURES
# ==========================================
print("\n[Step 7] Generating figures...")

# Figure 2: Forward Analysis
if len(forward_df) > 0:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    
    for i, direction in enumerate(['UP', 'DOWN']):
        subset = forward_df[forward_df['age_direction'] == direction]
        color = '#E74C3C' if direction == 'UP' else '#3498DB'
        axes[0].hist(subset['ARE_enrichment'].dropna(), bins=20, alpha=0.6, 
                     label=f'Age-{direction} (n={len(subset)})', color=color, edgecolor='black')
    axes[0].axvline(0, color='black', linestyle='--', linewidth=1.5)
    axes[0].set_xlabel('ARE Enrichment Score')
    axes[0].set_ylabel('Count')
    axes[0].set_title('A. Distribution')
    axes[0].legend()
    
    # Summary
    summary = forward_df.groupby('age_direction')['ARE_enrichment'].agg(['mean', 'std', 'count'])
    if len(summary) > 0:
        x_pos = range(len(summary))
        colors = ['#E74C3C' if d == 'UP' else '#3498DB' for d in summary.index]
        axes[2].bar(x_pos, summary['mean'], yerr=summary['std'], color=colors, edgecolor='black', capsize=5)
        axes[2].set_xticks(x_pos)
        axes[2].set_xticklabels([f'Age-{d}\n(n={int(summary.loc[d, "count"])})' for d in summary.index])
        axes[2].axhline(0, color='black', linestyle='--')
        axes[2].set_ylabel('Mean ARE Enrichment')
        axes[2].set_title('C. Summary')
    
    # Top genes
    top_n = min(8, len(forward_df))
    top_genes = pd.concat([forward_df.nlargest(top_n//2, 'ARE_enrichment'),
                           forward_df.nsmallest(top_n//2, 'ARE_enrichment')])
    colors = ['#27AE60' if x > 0 else '#C0392B' for x in top_genes['ARE_enrichment']]
    y_pos = range(len(top_genes))
    axes[1].barh(y_pos, top_genes['ARE_enrichment'], color=colors, edgecolor='black')
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels(top_genes['target_gene'])
    axes[1].axvline(0, color='black', linestyle='-')
    axes[1].set_xlabel('ARE Enrichment')
    axes[1].set_title('B. Top Activators/Suppressors')
    axes[1].invert_yaxis()
    
    plt.suptitle('Forward Analysis (RPE1 Perturb-seq)', fontsize=12, fontweight='bold')
    plt.tight_layout()
    save_figure(fig, "fig2_forward_analysis")
    plt.close()

# Figure 4: Convergence Scatter
if len(aging_mim_df) > 0 and len(prg4_df) > 0:
    merged = aging_mim_df.merge(prg4_df[['target_gene', 'prg4_mimetic_score']], on='target_gene')
    merged = merged.dropna(subset=['mimetic_score', 'prg4_mimetic_score'])
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Subsample background
    np.random.seed(42)
    background_mask = ~(merged['target_gene'].isin(are_genes) | merged['target_gene'].isin(convergent_genes))
    n_sample = min(1000, background_mask.sum())
    if n_sample > 0:
        background = merged[background_mask].sample(n=n_sample)
        ax.scatter(-background['mimetic_score'], background['prg4_mimetic_score'],
                   alpha=0.2, s=10, c='#95A5A6', label=f'Other (n={n_sample})')
    
    # ARE genes
    are_mask = merged['target_gene'].isin(are_genes)
    ax.scatter(-merged.loc[are_mask, 'mimetic_score'], merged.loc[are_mask, 'prg4_mimetic_score'],
               alpha=0.7, s=40, c='#E74C3C', marker='s', label=f'ARE genes (n={are_mask.sum()})')
    
    # Convergent
    conv_mask = merged['target_gene'].isin(convergent_genes)
    conv_data = merged[conv_mask]
    ax.scatter(-conv_data['mimetic_score'], conv_data['prg4_mimetic_score'],
               alpha=1, s=100, c='#3498DB', marker='D', edgecolors='black', linewidth=1.5,
               label=f'Convergent (n={len(conv_data)})')
    
    ax.axhline(0, color='black', linestyle=':', alpha=0.5)
    ax.axvline(0, color='black', linestyle=':', alpha=0.5)
    ax.set_xlabel('Aging Antagonist Score', fontsize=12)
    ax.set_ylabel('PRG4 Mimetic Score', fontsize=12)
    ax.set_title('Convergence Analysis (RPE1 Perturb-seq)', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left')
    
    plt.tight_layout()
    save_figure(fig, "fig4_convergence_scatter")
    plt.close()

# Figure 6: Summary
if len(convergent_genes) >= 0:
    fig, ax = plt.subplots(figsize=(8, 5))
    
    labels = ['Antagonists', 'PRG4 Mim.', 'ARE genes', 'Convergent']
    values = [len(top_antagonists), len(top_prg4_mimetics), 
              len(are_genes & set(aging_mim_df['target_gene'])) if len(aging_mim_df) > 0 else 0,
              len(convergent_genes)]
    colors = ['#27AE60', '#3498DB', '#E74C3C', '#9B59B6']
    
    bars = ax.bar(range(len(labels)), values, color=colors, edgecolor='black')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylabel('Number of Genes')
    ax.set_title('Summary (RPE1 Perturb-seq)', fontsize=12, fontweight='bold')
    
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, str(val),
                ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    save_figure(fig, "fig6_summary_overview")
    plt.close()

# ==========================================
# SUMMARY
# ==========================================
print("\n" + "=" * 60)
print("RPE1 ANALYSIS COMPLETE")
print("=" * 60)

print(f"\nOutput directory: {OUTPUT_DIR}")
print(f"\nKey results:")
print(f"  Forward analysis genes: {len(forward_df) if len(forward_df) > 0 else 'N/A'}")
print(f"  Total knockdowns analyzed: {len(aging_mim_df) if len(aging_mim_df) > 0 else 'N/A'}")
print(f"  Convergent genes: {len(convergent_genes)}")
if len(convergent_genes) > 0:
    print(f"  Convergent gene list: {sorted(convergent_genes)}")
