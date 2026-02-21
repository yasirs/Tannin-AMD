"""
ARE vs Gene Set Scatter Plot Module

For each gene set × dataset combination, create a scatter plot:
- X-axis: ARE gene set mimetic score
- Y-axis: Query gene set mimetic score

Each point represents a knockdown. Highlights "good" scores (top antagonists).
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

from config import OUTPUT_DIRS, VIZ_SETTINGS, ensure_dirs, get_dataset_output_dir
from perturbseq_loader import load_perturbseq_data, get_gene_indices
from gene_set_registry import load_all_gene_sets
from enrichment_analysis import compute_mimetic_score
from visualization import setup_style, save_figure

setup_style()

#%%
def compute_are_mimetic_scores(data: dict, are_genes: list) -> pd.DataFrame:
    """
    Compute ARE mimetic scores for all perturbations.
    
    ARE is a "mixed" gene set (not clearly UP/DOWN), so we use
    enrichment of ARE genes as a proxy for antioxidant pathway activity.
    
    Higher expression of ARE genes = positive "ARE score"
    """
    from perturbseq_loader import get_gene_indices
    
    expression = data["expression"]
    gene_names = data["gene_names"]
    target_genes = data["target_genes"]
    perturbation_ids = data["perturbation_ids"]
    
    # Get ARE gene indices
    are_indices, matched_are, _ = get_gene_indices(gene_names, are_genes)
    
    print(f"  ARE genes mapped: {len(matched_are)}/{len(are_genes)}")
    
    n_perts = len(perturbation_ids)
    results = []
    
    for i in range(n_perts):
        # Skip non-targeting controls
        target = str(target_genes[i]).lower()
        if "non-targeting" in target or "nontargeting" in target or target == "control":
            continue
        
        expr = expression[i, :]
        
        # For ARE, we compute simple enrichment (how much ARE genes are up vs all genes)
        # Use mean z-score approach
        are_expr = expr[are_indices]
        are_score = np.mean(are_expr) - np.mean(expr)  # Simple mean difference
        
        results.append({
            "perturbation_id": perturbation_ids[i],
            "target_gene": target_genes[i],
            "are_score": are_score,
        })
    
    return pd.DataFrame(results)

#%%
def plot_are_vs_geneset_scatter(
    geneset_df: pd.DataFrame,
    are_df: pd.DataFrame,
    gene_set_name: str,
    dataset_name: str,
    n_highlight: int = 50,
    output_path: Path = None
):
    """
    Create scatter plot: ARE score (x) vs Gene Set mimetic score (y).
    
    Highlights top antagonists (negative mimetic scores).
    """
    # Merge on target_gene
    merged = geneset_df.merge(
        are_df[["target_gene", "are_score"]],
        on="target_gene",
        how="inner"
    )
    
    if len(merged) < 10:
        print(f"  Only {len(merged)} overlapping genes - skipping")
        return None, None
    
    x = merged["are_score"].values
    y = merged["mimetic_score"].values
    
    # Identify top antagonists (most negative mimetic scores)
    top_antagonist_genes = set(
        merged.nsmallest(n_highlight, "mimetic_score")["target_gene"].str.upper()
    )
    merged["is_top_antagonist"] = merged["target_gene"].str.upper().isin(top_antagonist_genes)
    
    # Also mark top 50 high ARE + low mimetic (therapeutic quadrant)
    merged["therapeutic_score"] = -merged["mimetic_score"] + merged["are_score"] * 0.5
    top_therapeutic = set(
        merged.nlargest(n_highlight, "therapeutic_score")["target_gene"].str.upper()
    )
    merged["is_therapeutic"] = merged["target_gene"].str.upper().isin(top_therapeutic)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(5, 5))
    
    # All points (gray)
    ax.scatter(x, y, c=VIZ_SETTINGS["color_nonsignificant"], s=5, alpha=0.3, rasterized=True)
    
    # Top antagonists (green - negative mimetic score = reverses signature)
    antagonist_mask = merged["is_top_antagonist"]
    ax.scatter(x[antagonist_mask], y[antagonist_mask], 
              c='#4DAF4A', s=20, alpha=0.7, 
              label=f'Top {n_highlight} antagonists', zorder=5)
    
    # Top therapeutic (purple - high ARE + low mimetic)
    therapeutic_mask = merged["is_therapeutic"] & ~merged["is_top_antagonist"]
    ax.scatter(x[therapeutic_mask], y[therapeutic_mask], 
              c='#984EA3', s=25, alpha=0.8, 
              edgecolors='black', linewidths=0.3,
              label='Therapeutic quadrant', zorder=10)
    
    # Label top 5 therapeutic candidates
    top_to_label = merged[merged["is_therapeutic"]].nlargest(5, "therapeutic_score")
    for _, row in top_to_label.iterrows():
        ax.annotate(row["target_gene"], 
                   (row["are_score"], row["mimetic_score"]),
                   fontsize=7, xytext=(5, 5), textcoords='offset points')
    
    # Reference lines
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # Shade therapeutic quadrant (high ARE, low mimetic)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.fill_between([0, xlim[1]], [ylim[0], ylim[0]], [0, 0], 
                   alpha=0.05, color='green', zorder=0)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    # Correlation
    corr, pval = stats.spearmanr(x, y)
    
    ax.set_xlabel('ARE Pathway Score (higher = antioxidant active)')
    ax.set_ylabel(f'{gene_set_name} Mimetic Score\n(negative = signature reversed)')
    ax.set_title(f'{gene_set_name} vs ARE ({dataset_name})\nr={corr:.3f}, p={pval:.2e}')
    ax.legend(loc='upper left', fontsize=7)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def run_all_are_vs_geneset_plots():
    """
    Generate ARE vs Gene Set scatter plots for all combinations.
    """
    print("="*60)
    print("ARE vs GENE SET SCATTER PLOTS")
    print("="*60)
    
    ensure_dirs()
    
    # Load gene sets
    gene_sets = load_all_gene_sets()
    are_genes = gene_sets.get("ARE", [])
    
    if len(are_genes) == 0:
        print("ERROR: ARE gene set not found")
        return
    
    print(f"\nARE gene set: {len(are_genes)} genes")
    
    # Load Perturb-seq data
    rpe1_data = load_perturbseq_data("RPE1", use_cache=True)
    k562_data = load_perturbseq_data("K562_GWPS", use_cache=True)
    
    # Compute ARE scores for each dataset
    print("\nComputing ARE scores...")
    rpe1_are = compute_are_mimetic_scores(rpe1_data, are_genes)
    k562_are = compute_are_mimetic_scores(k562_data, are_genes)
    
    print(f"  RPE1: {len(rpe1_are)} perturbations")
    print(f"  K562: {len(k562_are)} perturbations")
    
    # Gene sets to plot against ARE
    gene_set_pairs = [
        ("Age", "Age"),
        ("AMD-Macula", "AMD-Macula"),
        ("AMD-nonMacula", "AMD-nonMacula"),
        ("Senescence", "Senescence"),
        ("Redox", "Redox"),
        ("Inflammation", "Inflammation"),
    ]
    
    # Output directory
    out_dir = OUTPUT_DIRS["concordance"] / "are_vs_geneset"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    for pair_name, display_name in gene_set_pairs:
        print(f"\n[{display_name}]")
        
        for dataset_name, data, are_df in [
            ("RPE1", rpe1_data, rpe1_are),
            ("K562_GWPS", k562_data, k562_are)
        ]:
            print(f"  {dataset_name}:")
            
            # Load backward analysis results for this gene set
            backward_dir = get_dataset_output_dir(pair_name, dataset_name, "backward")
            backward_file = backward_dir / f"backward_{pair_name}_mimetics.csv"
            
            if not backward_file.exists():
                print(f"    Backward file not found: {backward_file}")
                continue
            
            geneset_df = pd.read_csv(backward_file)
            
            # Create scatter plot
            plot_are_vs_geneset_scatter(
                geneset_df, are_df, display_name, dataset_name,
                output_path=out_dir / f"Fig_{pair_name}_vs_ARE_{dataset_name}"
            )
    
    print("\n" + "="*60)
    print("ARE vs Gene Set scatter plots complete!")
    print(f"Output directory: {out_dir}")
    print("="*60)

#%%
if __name__ == "__main__":
    run_all_are_vs_geneset_plots()
