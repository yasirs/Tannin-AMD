"""
Complete Forward Analysis Module

Reproduces Figure 2 from Age-NRF2 report:
For each gene in the query set that has knockdown data,
compute enrichment of a TARGET gene set (e.g., ARE).
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
from perturbseq_loader import load_perturbseq_data, get_gene_indices, get_perturbation_indices
from enrichment_analysis import compute_enrichment_for_perturbation
from visualization import setup_style, save_figure

setup_style()

#%%
def run_complete_forward_analysis(
    perturbseq_data: dict,
    up_genes: list,
    down_genes: list,
    target_genes: list,
    gene_set_name: str,
    target_name: str = "ARE"
) -> pd.DataFrame:
    """
    Forward Analysis: What happens to TARGET genes when we knock down UP/DOWN genes?
    
    Parameters
    ----------
    perturbseq_data : dict
        Loaded Perturb-seq data
    up_genes : list
        Genes upregulated in the signature (e.g., Age-UP)
    down_genes : list
        Genes downregulated in the signature (e.g., Age-DOWN)
    target_genes : list
        Target gene set to measure enrichment (e.g., ARE)
    gene_set_name : str
        Name for labeling (e.g., "Age")
    target_name : str
        Name of target set (e.g., "ARE")
    
    Returns
    -------
    pd.DataFrame with enrichment scores for each knockdown
    """
    expression = perturbseq_data["expression"]
    gene_names = perturbseq_data["gene_names"]
    target_genes_ps = perturbseq_data["target_genes"]
    perturbation_ids = perturbseq_data["perturbation_ids"]
    
    # Get target gene indices in expression matrix
    target_indices, matched_targets, _ = get_gene_indices(gene_names, target_genes)
    
    if len(target_indices) < 5:
        print(f"  WARNING: Only {len(target_indices)} target genes mapped - results may be unreliable")
    
    # Get perturbations for UP genes
    up_pert_indices, matched_up = get_perturbation_indices(target_genes_ps, up_genes)
    
    # Get perturbations for DOWN genes
    down_pert_indices, matched_down = get_perturbation_indices(target_genes_ps, down_genes)
    
    print(f"  Forward Analysis: {gene_set_name} → {target_name}")
    print(f"    UP gene knockdowns available: {len(matched_up)}/{len(up_genes)}")
    print(f"    DOWN gene knockdowns available: {len(matched_down)}/{len(down_genes)}")
    print(f"    Target genes mapped: {len(target_indices)}/{len(target_genes)}")
    
    results = []
    
    # Process UP gene knockdowns
    for idx in up_pert_indices:
        expr = expression[idx, :]
        es = compute_enrichment_for_perturbation(expr, target_indices)
        
        results.append({
            "perturbation_id": perturbation_ids[idx],
            "knockdown_gene": target_genes_ps[idx],
            "direction": "UP",
            f"{target_name}_enrichment": es,
        })
    
    # Process DOWN gene knockdowns
    for idx in down_pert_indices:
        expr = expression[idx, :]
        es = compute_enrichment_for_perturbation(expr, target_indices)
        
        results.append({
            "perturbation_id": perturbation_ids[idx],
            "knockdown_gene": target_genes_ps[idx],
            "direction": "DOWN",
            f"{target_name}_enrichment": es,
        })
    
    df = pd.DataFrame(results)
    
    if len(df) > 0:
        # Summary stats
        up_scores = df[df["direction"] == "UP"][f"{target_name}_enrichment"]
        down_scores = df[df["direction"] == "DOWN"][f"{target_name}_enrichment"]
        
        print(f"    UP knockdowns: mean {target_name} enrichment = {up_scores.mean():.4f} ± {up_scores.std():.4f}")
        print(f"    DOWN knockdowns: mean {target_name} enrichment = {down_scores.mean():.4f} ± {down_scores.std():.4f}")
    else:
        print(f"    No knockdown data available for query genes")
    
    return df

#%%
def plot_forward_distribution(
    df: pd.DataFrame,
    target_name: str = "ARE",
    gene_set_name: str = "",
    output_path: Path = None
):
    """
    Plot Figure 2-style distribution: UP vs DOWN knockdowns effect on target.
    """
    if len(df) == 0:
        print(f"    No data to plot for {gene_set_name}")
        return None, None
    
    score_col = f"{target_name}_enrichment"
    
    fig, ax = plt.subplots(figsize=(5, 4))
    
    # Get UP and DOWN scores
    up_mask = df["direction"] == "UP"
    down_mask = df["direction"] == "DOWN"
    
    up_scores = df.loc[up_mask, score_col].values
    down_scores = df.loc[down_mask, score_col].values
    
    # Plot histograms
    bins = np.linspace(-0.5, 0.5, 31)
    
    if len(up_scores) > 0:
        ax.hist(up_scores, bins=bins, alpha=0.6, label=f'UP knockdowns (n={len(up_scores)})', 
               color='#E41A1C', density=True)
    
    if len(down_scores) > 0:
        ax.hist(down_scores, bins=bins, alpha=0.6, label=f'DOWN knockdowns (n={len(down_scores)})', 
               color='#377EB8', density=True)
    
    # Reference line at zero
    ax.axvline(0, color='gray', linestyle='--', linewidth=1)
    
    # Mean lines
    if len(up_scores) > 0:
        ax.axvline(up_scores.mean(), color='#E41A1C', linestyle='-', linewidth=2, alpha=0.7)
    if len(down_scores) > 0:
        ax.axvline(down_scores.mean(), color='#377EB8', linestyle='-', linewidth=2, alpha=0.7)
    
    ax.set_xlabel(f'{target_name} Enrichment Score')
    ax.set_ylabel('Density')
    ax.set_title(f'{gene_set_name} Knockdowns → {target_name} Effect')
    ax.legend(loc='upper right', fontsize=8)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def get_top_mediators(
    df: pd.DataFrame,
    target_name: str = "ARE",
    n_top: int = 6
) -> pd.DataFrame:
    """
    Get top activators and suppressors from forward analysis.
    
    Returns Table 1-style summary.
    """
    if len(df) == 0:
        return pd.DataFrame()
    
    score_col = f"{target_name}_enrichment"
    
    results = []
    
    for direction in ["UP", "DOWN"]:
        subset = df[df["direction"] == direction].copy()
        if len(subset) == 0:
            continue
        
        # Top activators (positive enrichment)
        top_act = subset.nlargest(n_top, score_col)
        for _, row in top_act.iterrows():
            results.append({
                "category": f"Query-{direction}",
                "gene": row["knockdown_gene"],
                f"{target_name}_effect": "Activates" if row[score_col] > 0 else "Suppresses",
                "enrichment_score": row[score_col],
            })
        
        # Top suppressors (negative enrichment)
        top_sup = subset.nsmallest(n_top, score_col)
        for _, row in top_sup.iterrows():
            if row["knockdown_gene"] not in [r["gene"] for r in results]:  # Avoid duplicates
                results.append({
                    "category": f"Query-{direction}",
                    "gene": row["knockdown_gene"],
                    f"{target_name}_effect": "Activates" if row[score_col] > 0 else "Suppresses",
                    "enrichment_score": row[score_col],
                })
    
    return pd.DataFrame(results)

#%%
def run_all_forward_analyses(use_cache: bool = True):
    """
    Run forward analysis for all gene set pairs → ARE.
    """
    from gene_set_registry import load_all_gene_sets
    from run_analysis import get_paired_gene_sets
    
    print("="*60)
    print("FORWARD ANALYSIS: Query KDs → ARE Enrichment")
    print("="*60)
    
    ensure_dirs()
    
    # Load gene sets
    gene_sets = load_all_gene_sets()
    
    # Get ARE genes as target
    are_genes = gene_sets.get("ARE", [])
    if len(are_genes) == 0:
        print("ERROR: ARE gene set not found")
        return
    
    print(f"\nTarget: ARE gene set ({len(are_genes)} genes)")
    
    # Get paired gene sets
    pairs = get_paired_gene_sets(gene_sets)
    
    # Load Perturb-seq data
    rpe1_data = load_perturbseq_data("RPE1", use_cache=use_cache)
    k562_data = load_perturbseq_data("K562_GWPS", use_cache=use_cache)
    
    all_results = {}
    
    for pair_name, pair_info in pairs.items():
        up_genes = pair_info["up"]
        down_genes = pair_info["down"]
        
        for dataset_name, data in [("RPE1", rpe1_data), ("K562_GWPS", k562_data)]:
            print(f"\n{'='*50}")
            print(f"{pair_name} on {dataset_name}")
            print("="*50)
            
            # Get output directory
            from config import get_dataset_output_dir
            out_dir = get_dataset_output_dir(pair_name, dataset_name, "forward")
            
            # Run forward analysis
            forward_df = run_complete_forward_analysis(
                data,
                up_genes,
                down_genes,
                are_genes,
                gene_set_name=pair_name,
                target_name="ARE"
            )
            
            # Save results
            if len(forward_df) > 0:
                forward_df.to_csv(out_dir / f"forward_{pair_name}_ARE_enrichment.csv", index=False)
                
                # Plot distribution
                plot_forward_distribution(
                    forward_df,
                    target_name="ARE",
                    gene_set_name=pair_name,
                    output_path=out_dir / f"Fig_forward_{pair_name}_ARE_distribution"
                )
                
                # Get top mediators
                mediators_df = get_top_mediators(forward_df, target_name="ARE")
                if len(mediators_df) > 0:
                    mediators_df.to_csv(out_dir / f"forward_{pair_name}_top_mediators.csv", index=False)
                
                all_results[f"{pair_name}_{dataset_name}"] = forward_df
            else:
                print(f"  No results for {pair_name} on {dataset_name}")
    
    print("\n" + "="*60)
    print("Forward analysis complete!")
    print("="*60)
    
    return all_results

#%%
if __name__ == "__main__":
    run_all_forward_analyses()
