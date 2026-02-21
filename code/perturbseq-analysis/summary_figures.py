"""
Summary Overview Figure Module

Reproduces Figure 6 from Age-NRF2 report:
- Gene set sizes and coverage
- Key overlaps visualization
- Ranked convergent genes
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

from config import OUTPUT_DIRS, VIZ_SETTINGS, ensure_dirs
from visualization import setup_style, save_figure

setup_style()

#%%
def plot_coverage_summary(output_path: Path = None):
    """
    Bar chart of gene set coverage in Perturb-seq data.
    """
    coverage_file = OUTPUT_DIRS["gene_sets"] / "coverage_in_perturbseq.csv"
    
    if not coverage_file.exists():
        print(f"Coverage file not found: {coverage_file}")
        return None, None
    
    coverage_df = pd.read_csv(coverage_file)
    
    # Pivot for easier plotting
    pivot_df = coverage_df.pivot(index="gene_set_name", columns="dataset", values="gene_coverage_pct")
    
    # Sort by average coverage
    pivot_df["avg"] = pivot_df.mean(axis=1)
    pivot_df = pivot_df.sort_values("avg", ascending=True)
    pivot_df = pivot_df.drop("avg", axis=1)
    
    fig, ax = plt.subplots(figsize=(6, 5))
    
    x = np.arange(len(pivot_df))
    width = 0.35
    
    rpe1_vals = pivot_df.get("RPE1", [0] * len(pivot_df))
    k562_vals = pivot_df.get("K562_GWPS", [0] * len(pivot_df))
    
    bars1 = ax.barh(x - width/2, rpe1_vals, width, label='RPE1', color='#E41A1C', alpha=0.8)
    bars2 = ax.barh(x + width/2, k562_vals, width, label='K562', color='#377EB8', alpha=0.8)
    
    ax.set_yticks(x)
    ax.set_yticklabels(pivot_df.index, fontsize=8)
    ax.set_xlabel('Gene Coverage (%)')
    ax.set_title('Gene Set Coverage in Perturb-seq Data')
    ax.legend(loc='lower right')
    
    # Add value labels
    for bar in bars1:
        width_val = bar.get_width()
        if width_val > 5:
            ax.text(width_val - 3, bar.get_y() + bar.get_height()/2,
                   f'{width_val:.0f}%', va='center', ha='right', fontsize=6, color='white')
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def plot_concordance_heatmap(output_path: Path = None):
    """
    Heatmap of RPE1 vs K562 concordance for all gene sets.
    """
    conc_file = OUTPUT_DIRS["concordance"] / "rpe1_vs_k562_concordance.csv"
    
    if not conc_file.exists():
        print(f"Concordance file not found: {conc_file}")
        return None, None
    
    conc_df = pd.read_csv(conc_file)
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    
    # Panel A: Overlap counts
    ax1 = axes[0]
    gene_sets = conc_df["gene_set"].values
    overlaps = conc_df["overlap"].values
    
    colors = ['#4DAF4A' if p < 0.05 else '#CCCCCC' for p in conc_df["hypergeom_pval"]]
    ax1.barh(range(len(gene_sets)), overlaps, color=colors)
    ax1.set_yticks(range(len(gene_sets)))
    ax1.set_yticklabels(gene_sets)
    ax1.set_xlabel('Overlap (top 200 antagonists)')
    ax1.set_title('RPE1-K562 Antagonist Overlap')
    
    # Add significance markers
    for i, (overlap, pval) in enumerate(zip(overlaps, conc_df["hypergeom_pval"])):
        if pval < 0.001:
            marker = "***"
        elif pval < 0.01:
            marker = "**"
        elif pval < 0.05:
            marker = "*"
        else:
            marker = ""
        ax1.text(overlap + 0.5, i, f"{overlap}{marker}", va='center', fontsize=8)
    
    # Panel B: Correlation
    ax2 = axes[1]
    corrs = conc_df["spearman_corr"].values
    
    colors2 = ['#377EB8' if p < 0.05 else '#CCCCCC' for p in conc_df["corr_pval"]]
    ax2.barh(range(len(gene_sets)), corrs, color=colors2)
    ax2.set_yticks(range(len(gene_sets)))
    ax2.set_yticklabels(gene_sets)
    ax2.set_xlabel('Spearman Correlation')
    ax2.set_title('Score Correlation (all shared genes)')
    ax2.axvline(0, color='gray', linestyle='--', linewidth=0.5)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, axes

#%%
def plot_convergent_genes_ranked(n_top: int = 20, output_path: Path = None):
    """
    Bar chart ranking convergent genes by number of appearances.
    """
    conv_file = OUTPUT_DIRS["concordance"] / "convergent_genes_all.csv"
    
    if not conv_file.exists():
        print(f"Convergent genes file not found: {conv_file}")
        return None, None
    
    conv_df = pd.read_csv(conv_file)
    
    if len(conv_df) == 0:
        print("No convergent genes found")
        return None, None
    
    # Top N by appearances
    top_df = conv_df.nlargest(n_top, "n_appearances")
    
    fig, ax = plt.subplots(figsize=(6, 5))
    
    y = range(len(top_df))
    ax.barh(y, top_df["n_appearances"].values, color='#984EA3', alpha=0.8)
    
    ax.set_yticks(y)
    ax.set_yticklabels(top_df["gene"].values, fontsize=8)
    ax.set_xlabel('Number of Analyses')
    ax.set_title(f'Top {n_top} Convergent Genes\n(appearing in multiple analyses)')
    
    ax.invert_yaxis()  # Top at top
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def create_summary_table():
    """
    Create Table S1-style summary statistics.
    """
    summary = {
        "Metric": [],
        "Value": []
    }
    
    # Gene set sizes from registry
    registry_file = OUTPUT_DIRS["gene_sets"] / "registry.csv"
    if registry_file.exists():
        registry = pd.read_csv(registry_file)
        for _, row in registry.iterrows():
            summary["Metric"].append(f"{row['name']} genes")
            summary["Value"].append(row["size"])
    
    # Perturb-seq stats
    summary["Metric"].extend(["RPE1 knockdowns", "K562 GWPS knockdowns"])
    summary["Value"].extend([2679, 11258])
    
    # Convergent genes
    conv_file = OUTPUT_DIRS["concordance"] / "convergent_genes_all.csv"
    if conv_file.exists():
        conv_df = pd.read_csv(conv_file)
        summary["Metric"].append("Convergent genes (≥2 analyses)")
        summary["Value"].append(len(conv_df))
        
        # Genes in ≥4 analyses
        high_conf = len(conv_df[conv_df["n_appearances"] >= 4])
        summary["Metric"].append("High-confidence convergent (≥4 analyses)")
        summary["Value"].append(high_conf)
    
    df = pd.DataFrame(summary)
    
    # Save
    df.to_csv(OUTPUT_DIRS["concordance"] / "summary_statistics.csv", index=False)
    
    return df

#%%
def run_summary_generation():
    """
    Generate all summary figures and tables.
    """
    print("="*60)
    print("SUMMARY OVERVIEW GENERATION")
    print("="*60)
    
    ensure_dirs()
    
    out_dir = OUTPUT_DIRS["concordance"]
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Coverage summary
    print("\n[1] Coverage summary...")
    plot_coverage_summary(output_path=out_dir / "Fig_coverage_summary")
    
    # 2. Concordance heatmap
    print("\n[2] Concordance heatmap...")
    plot_concordance_heatmap(output_path=out_dir / "Fig_concordance_heatmap")
    
    # 3. Convergent genes ranked
    print("\n[3] Convergent genes ranked...")
    plot_convergent_genes_ranked(n_top=20, output_path=out_dir / "Fig_convergent_genes_ranked")
    
    # 4. Summary table
    print("\n[4] Summary statistics table...")
    summary_df = create_summary_table()
    print(summary_df.to_string(index=False))
    
    print("\n" + "="*60)
    print("Summary generation complete!")
    print("="*60)

#%%
if __name__ == "__main__":
    run_summary_generation()
