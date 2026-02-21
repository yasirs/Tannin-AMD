"""
Single Gene Set Scatter Plots

Creates scatter plots with single gene set ES on each axis:
- X-axis: Gene Set A enrichment score
- Y-axis: Gene Set B enrichment score
Each point = one knockdown

Uses fgsea results which have ES for each gene set independently.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))

from config import OUTPUT_DIRS, VIZ_SETTINGS
from visualization import setup_style, save_figure

setup_style()

#%%
def load_fgsea_pivoted(dataset_name: str) -> pd.DataFrame:
    """
    Load fgsea results and pivot to wide format.
    
    Returns DataFrame with columns:
    - target_gene
    - perturbation_id
    - {pathway}_ES, {pathway}_pval, {pathway}_padj for each pathway
    """
    fgsea_file = OUTPUT_DIRS["concordance"].parent / "fgsea" / f"fgsea_{dataset_name}_results.csv"
    
    if not fgsea_file.exists():
        raise FileNotFoundError(f"fgsea results not found: {fgsea_file}")
    
    fgsea = pd.read_csv(fgsea_file)
    
    # Pivot to wide format
    es_pivot = fgsea.pivot_table(
        index=['perturbation_id', 'target_gene'],
        columns='pathway',
        values='ES',
        aggfunc='first'
    ).reset_index()
    es_pivot.columns = ['perturbation_id', 'target_gene'] + [f'{c}_ES' for c in es_pivot.columns[2:]]
    
    pval_pivot = fgsea.pivot_table(
        index=['perturbation_id', 'target_gene'],
        columns='pathway',
        values='pval',
        aggfunc='first'
    ).reset_index()
    pval_pivot.columns = ['perturbation_id', 'target_gene'] + [f'{c}_pval' for c in pval_pivot.columns[2:]]
    
    padj_pivot = fgsea.pivot_table(
        index=['perturbation_id', 'target_gene'],
        columns='pathway',
        values='padj',
        aggfunc='first'
    ).reset_index()
    padj_pivot.columns = ['perturbation_id', 'target_gene'] + [f'{c}_padj' for c in padj_pivot.columns[2:]]
    
    # Merge all
    result = es_pivot.merge(pval_pivot, on=['perturbation_id', 'target_gene'])
    result = result.merge(padj_pivot, on=['perturbation_id', 'target_gene'])
    
    return result

#%%
def plot_single_geneset_scatter(
    df: pd.DataFrame,
    geneset_x: str,
    geneset_y: str,
    dataset_name: str,
    pval_threshold: float = 0.05,
    n_label: int = 5,
    output_path: Path = None
):
    """
    Create scatter plot: geneset_x ES (x) vs geneset_y ES (y).
    
    Highlights significant knockdowns based on p-values.
    """
    x_col = f"{geneset_x}_ES"
    y_col = f"{geneset_y}_ES"
    x_pval = f"{geneset_x}_pval"
    y_pval = f"{geneset_y}_pval"
    
    if x_col not in df.columns or y_col not in df.columns:
        print(f"  Missing columns: {x_col} or {y_col}")
        return None, None
    
    # Remove NaN
    valid_mask = df[x_col].notna() & df[y_col].notna()
    df_valid = df[valid_mask].copy()
    
    x = df_valid[x_col].values
    y = df_valid[y_col].values
    
    # Significance masks
    sig_x = df_valid[x_pval] < pval_threshold if x_pval in df_valid.columns else np.zeros(len(df_valid), dtype=bool)
    sig_y = df_valid[y_pval] < pval_threshold if y_pval in df_valid.columns else np.zeros(len(df_valid), dtype=bool)
    sig_both = sig_x & sig_y
    
    # Create figure
    fig, ax = plt.subplots(figsize=(5, 5))
    
    # All points (gray)
    ax.scatter(x, y, c=VIZ_SETTINGS["color_nonsignificant"], s=5, alpha=0.3, rasterized=True)
    
    # Significant in X only (light blue)
    mask_x = sig_x & ~sig_y
    if mask_x.sum() > 0:
        ax.scatter(x[mask_x], y[mask_x], c='#B3CDE3', s=15, alpha=0.5, 
                  label=f'Sig {geneset_x} only')
    
    # Significant in Y only (light red)
    mask_y = sig_y & ~sig_x
    if mask_y.sum() > 0:
        ax.scatter(x[mask_y], y[mask_y], c='#FBB4AE', s=15, alpha=0.5,
                  label=f'Sig {geneset_y} only')
    
    # Significant in both (purple/highlighted)
    if sig_both.sum() > 0:
        ax.scatter(x[sig_both], y[sig_both], c='#984EA3', s=25, alpha=0.8,
                  edgecolors='black', linewidths=0.3,
                  label=f'Sig both (n={sig_both.sum()})', zorder=10)
    
    # Label top knockdowns significant in both
    if sig_both.sum() > 0:
        both_df = df_valid[sig_both].copy()
        both_df['combined'] = abs(both_df[x_col]) + abs(both_df[y_col])
        top_to_label = both_df.nlargest(n_label, 'combined')
        
        for _, row in top_to_label.iterrows():
            ax.annotate(row['target_gene'],
                       (row[x_col], row[y_col]),
                       fontsize=6, xytext=(5, 5), textcoords='offset points')
    
    # Reference lines
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # Correlation
    valid_both = np.isfinite(x) & np.isfinite(y)
    corr, pval = stats.spearmanr(x[valid_both], y[valid_both])
    
    ax.set_xlabel(f'{geneset_x} Enrichment Score')
    ax.set_ylabel(f'{geneset_y} Enrichment Score')
    ax.set_title(f'{geneset_x} vs {geneset_y} ({dataset_name})\nr={corr:.3f}, p={pval:.2e}')
    ax.legend(loc='upper left', fontsize=6)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def run_all_single_geneset_plots():
    """
    Generate all single gene set scatter plots.
    """
    print("="*60)
    print("SINGLE GENE SET SCATTER PLOTS (fgsea-based)")
    print("="*60)
    
    # Output directory
    out_dir = OUTPUT_DIRS["concordance"] / "single_geneset_scatter"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Gene set pairs to plot
    pairs = [
        # ARE vs each axis
        ("ARE", "Redox-PRO"),
        ("ARE", "Redox-ANTI"),
        ("ARE", "Senescence-PRO"),
        ("ARE", "Senescence-ANTI"),
        ("ARE", "Inflammation-PRO"),
        ("ARE", "Inflammation-ANTI"),
        
        # Cross-axis comparisons
        ("Redox-PRO", "Senescence-PRO"),
        ("Redox-PRO", "Inflammation-PRO"),
        ("Senescence-PRO", "Inflammation-PRO"),
        
        # Age vs ARE
        ("Age-UP", "ARE"),
        ("Age-DOWN", "ARE"),
        
        # PRO vs ANTI within axis
        ("Redox-PRO", "Redox-ANTI"),
        ("Senescence-PRO", "Senescence-ANTI"),
        ("Inflammation-PRO", "Inflammation-ANTI"),
    ]
    
    for dataset in ["RPE1"]:  # Add K562_GWPS when available
        print(f"\n[{dataset}]")
        
        try:
            df = load_fgsea_pivoted(dataset)
            print(f"  Loaded {len(df)} knockdowns")
        except FileNotFoundError as e:
            print(f"  {e}")
            continue
        
        for gs_x, gs_y in pairs:
            print(f"  {gs_x} vs {gs_y}...")
            
            plot_single_geneset_scatter(
                df, gs_x, gs_y, dataset,
                output_path=out_dir / f"Fig_{gs_x}_vs_{gs_y}_{dataset}"
            )
    
    print("\n" + "="*60)
    print(f"Single gene set scatter plots complete!")
    print(f"Output directory: {out_dir}")
    print("="*60)

#%%
if __name__ == "__main__":
    run_all_single_geneset_plots()
