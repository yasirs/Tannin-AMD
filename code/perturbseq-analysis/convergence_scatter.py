"""
Convergence Scatter Plot Module

Reproduces Figure 4 from Age-NRF2 report:
Scatter plot showing antagonist scores from one analysis vs another,
with highlighted overlapping genes.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))

from config import OUTPUT_DIRS, ANALYSIS_PARAMS, VIZ_SETTINGS, ensure_dirs
from visualization import setup_style, save_figure

setup_style()

#%%
def plot_rpe1_vs_k562_scatter(
    rpe1_df: pd.DataFrame,
    k562_df: pd.DataFrame,
    gene_set_name: str,
    n_top: int = 200,
    output_path: Path = None
):
    """
    Create Figure 4-style scatter: RPE1 scores (x) vs K562 scores (y).
    
    Highlights genes in top N antagonists from both datasets.
    """
    # Merge on target gene
    merged = rpe1_df.merge(
        k562_df[["target_gene", "mimetic_score"]],
        on="target_gene",
        suffixes=("_rpe1", "_k562"),
        how="inner"
    )
    
    if len(merged) < 10:
        print(f"  Only {len(merged)} overlapping genes - insufficient for scatter")
        return None, None
    
    # Get top antagonists from each
    rpe1_top = set(rpe1_df.nsmallest(n_top, "mimetic_score")["target_gene"].str.upper())
    k562_top = set(k562_df.nsmallest(n_top, "mimetic_score")["target_gene"].str.upper())
    
    # Identify overlap
    overlap_genes = rpe1_top & k562_top
    
    # Create significance mask
    merged["is_rpe1_top"] = merged["target_gene"].str.upper().isin(rpe1_top)
    merged["is_k562_top"] = merged["target_gene"].str.upper().isin(k562_top)
    merged["is_overlap"] = merged["target_gene"].str.upper().isin(overlap_genes)
    
    # Plot
    fig, ax = plt.subplots(figsize=(5, 5))
    
    x = merged["mimetic_score_rpe1"].values
    y = merged["mimetic_score_k562"].values
    
    # All points (gray)
    ax.scatter(x, y, c=VIZ_SETTINGS["color_nonsignificant"], s=5, alpha=0.3, rasterized=True)
    
    # RPE1-only top genes (light red)
    rpe1_only = merged["is_rpe1_top"] & ~merged["is_overlap"]
    ax.scatter(x[rpe1_only], y[rpe1_only], c='#FBB4AE', s=15, alpha=0.6, label='RPE1 top only')
    
    # K562-only top genes (light blue)
    k562_only = merged["is_k562_top"] & ~merged["is_overlap"]
    ax.scatter(x[k562_only], y[k562_only], c='#B3CDE3', s=15, alpha=0.6, label='K562 top only')
    
    # Overlap genes (purple, highlighted)
    overlap_mask = merged["is_overlap"]
    ax.scatter(x[overlap_mask], y[overlap_mask], c='#984EA3', s=40, alpha=0.9, 
              edgecolors='black', linewidths=0.5, label=f'Both (n={overlap_mask.sum()})', zorder=10)
    
    # Label top overlap genes
    overlap_df = merged[overlap_mask].copy()
    if len(overlap_df) > 0:
        # Top 5 by combined rank
        overlap_df["combined_score"] = overlap_df["mimetic_score_rpe1"] + overlap_df["mimetic_score_k562"]
        top_to_label = overlap_df.nsmallest(5, "combined_score")
        
        for _, row in top_to_label.iterrows():
            ax.annotate(row["target_gene"], 
                       (row["mimetic_score_rpe1"], row["mimetic_score_k562"]),
                       fontsize=7, xytext=(5, 5), textcoords='offset points')
    
    # Reference lines
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # Correlation
    from scipy import stats
    corr, pval = stats.spearmanr(x, y)
    ax.set_title(f'{gene_set_name}: RPE1 vs K562\nr={corr:.3f}, p={pval:.2e}')
    
    ax.set_xlabel('Mimetic Score (RPE1)')
    ax.set_ylabel('Mimetic Score (K562)')
    ax.legend(loc='upper left', fontsize=7)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def plot_cross_axis_scatter(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    name1: str,
    name2: str,
    n_top: int = 200,
    output_path: Path = None
):
    """
    Scatter of mimetic scores from two different gene sets.
    
    E.g., Age antagonist score vs Senescence antagonist score.
    """
    # Merge on target gene
    merged = df1.merge(
        df2[["target_gene", "mimetic_score"]],
        on="target_gene",
        suffixes=(f"_{name1}", f"_{name2}"),
        how="inner"
    )
    
    if len(merged) < 10:
        print(f"  Only {len(merged)} overlapping genes for {name1} vs {name2}")
        return None, None
    
    # Get top antagonists from each
    top1 = set(df1.nsmallest(n_top, "mimetic_score")["target_gene"].str.upper())
    top2 = set(df2.nsmallest(n_top, "mimetic_score")["target_gene"].str.upper())
    overlap = top1 & top2
    
    merged["is_overlap"] = merged["target_gene"].str.upper().isin(overlap)
    
    # Plot
    fig, ax = plt.subplots(figsize=(5, 5))
    
    score1_col = f"mimetic_score_{name1}"
    score2_col = f"mimetic_score_{name2}"
    
    x = merged[score1_col].values
    y = merged[score2_col].values
    
    # All points (gray)
    ax.scatter(x, y, c=VIZ_SETTINGS["color_nonsignificant"], s=5, alpha=0.3, rasterized=True)
    
    # Overlap points (colored)
    overlap_mask = merged["is_overlap"]
    ax.scatter(x[overlap_mask], y[overlap_mask], c='#E41A1C', s=30, alpha=0.8, 
              edgecolors='black', linewidths=0.3, label=f'Both top {n_top} (n={overlap_mask.sum()})', zorder=10)
    
    # Reference lines
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # Correlation
    from scipy import stats
    corr, pval = stats.spearmanr(x, y)
    ax.set_title(f'{name1} vs {name2}\nr={corr:.3f}, p={pval:.2e}')
    
    ax.set_xlabel(f'Mimetic Score ({name1})')
    ax.set_ylabel(f'Mimetic Score ({name2})')
    ax.legend(loc='upper left', fontsize=8)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def run_all_scatter_plots():
    """
    Generate all convergence scatter plots.
    """
    print("="*60)
    print("CONVERGENCE SCATTER PLOTS")
    print("="*60)
    
    ensure_dirs()
    
    # Find all backward results
    results_dir = OUTPUT_DIRS["Age"].parent
    
    # Collect all backward analysis results
    backward_results = {"RPE1": {}, "K562_GWPS": {}}
    
    for dataset in ["RPE1", "K562_GWPS"]:
        # Search for backward mimetic files
        for csv_file in results_dir.rglob(f"*/{dataset}/backward/backward_*_mimetics.csv"):
            # Extract gene set name from filename
            gene_set = csv_file.stem.replace("backward_", "").replace("_mimetics", "")
            backward_results[dataset][gene_set] = pd.read_csv(csv_file)
            print(f"  Loaded: {gene_set} from {dataset}")
    
    # Create output directory
    conc_dir = OUTPUT_DIRS["concordance"]
    conc_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. RPE1 vs K562 scatter for each gene set
    print("\n[1] Creating RPE1 vs K562 scatter plots...")
    
    common_sets = set(backward_results["RPE1"].keys()) & set(backward_results["K562_GWPS"].keys())
    
    for gene_set in common_sets:
        rpe1_df = backward_results["RPE1"][gene_set]
        k562_df = backward_results["K562_GWPS"][gene_set]
        
        print(f"\n  {gene_set}:")
        plot_rpe1_vs_k562_scatter(
            rpe1_df, k562_df, gene_set,
            output_path=conc_dir / f"Fig_{gene_set}_rpe1_vs_k562_scatter"
        )
    
    # 2. Cross-axis scatter plots (within each dataset)
    print("\n[2] Creating cross-axis scatter plots...")
    
    axis_pairs = [
        ("Age", "Senescence"),
        ("Senescence", "Redox"),
        ("Redox", "Inflammation"),
        ("Age", "Inflammation"),
    ]
    
    for dataset in ["RPE1", "K562_GWPS"]:
        print(f"\n  {dataset}:")
        
        for name1, name2 in axis_pairs:
            if name1 in backward_results[dataset] and name2 in backward_results[dataset]:
                df1 = backward_results[dataset][name1]
                df2 = backward_results[dataset][name2]
                
                print(f"    {name1} vs {name2}")
                plot_cross_axis_scatter(
                    df1, df2, name1, name2,
                    output_path=conc_dir / f"Fig_{name1}_vs_{name2}_{dataset}_scatter"
                )
    
    print("\n" + "="*60)
    print("Scatter plots complete!")
    print("="*60)

#%%
if __name__ == "__main__":
    run_all_scatter_plots()
