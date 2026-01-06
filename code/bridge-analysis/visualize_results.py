import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

# Paths
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
DATA_DIR = os.path.join(BASE_DIR, "results/bridge-results")
PERTURB_FILE = os.path.join(BASE_DIR, "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad")
BULK_DIR = os.path.join(BASE_DIR, "data/RPE_cells/code")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/bridge-results")

os.makedirs(OUTPUT_DIR, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")

def plot_top_correlations(file_path, label):
    df = pd.read_csv(file_path)
    # Get top 10 and bottom 10
    top_10 = df.head(10).copy()
    bottom_10 = df.tail(10).copy()
    plot_df = pd.concat([top_10, bottom_10])
    
    plt.figure(figsize=(10, 8))
    colors = ["#d73027" if r > 0 else "#4575b4" for r in plot_df["spearman_rho"]]
    sns.barplot(data=plot_df, x="spearman_rho", y="target_gene", palette=colors)
    plt.title(f"Top Correlates: {label}")
    plt.xlabel("Spearman Rho (ρ)")
    plt.ylabel("Knockdown Target")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{label}_rankings.png"), dpi=300)
    plt.close()

def plot_bridge_evidence(bulk_file, bulk_col, perturb_kd, label, common_genes_file):
    # This function creates a scatter plot of bulk log2FC vs Perturb-seq Z-score
    # to show the 'bridge' in action.
    
    # Load bulk
    df_bulk = pd.read_csv(bulk_file)
    df_bulk = df_bulk.set_index("ensembl_gene_id")
    
    # Load perturb-seq
    ad = sc.read_h5ad(PERTURB_FILE)
    
    # Extract KD profile
    if perturb_kd not in ad.obs_names:
        # Try to find it by target gene if exact KD name is unknown
        matches = [o for o in ad.obs_names if f"_{perturb_kd}_" in o]
        if not matches:
             print(f"Could not find KD for {perturb_kd}")
             return
        perturb_kd = matches[0]
        
    # Get common genes
    common_genes = list(set(df_bulk.index) & set(ad.var_names))
    
    # Subset
    bulk_vec = df_bulk.loc[common_genes, bulk_col]
    ad_sub = ad[perturb_kd, common_genes]
    kd_vec = ad_sub.X.toarray().flatten() if hasattr(ad_sub.X, "toarray") else ad_sub.X.flatten()
    
    # Create plot
    plt.figure(figsize=(8, 8))
    sns.regplot(x=bulk_vec, y=kd_vec, scatter_kws={'alpha':0.3, 's':10, 'color':'gray'}, line_kws={'color':'red'})
    
    r, p = spearmanr(bulk_vec, kd_vec)
    plt.title(f"Bridge Evidence: {label}\n(Bulk {bulk_col} vs {perturb_kd})\nρ = {r:.3f}")
    plt.xlabel(f"Bulk log2FoldChange ({bulk_col})")
    plt.ylabel(f"Perturb-seq Z-score ({perturb_kd})")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{label}_bridge_scatter.png"), dpi=300)
    plt.close()

if __name__ == "__main__":
    # 1. Bar plots for all three analyses
    files_to_plot = [
        ("H2O2_Stress_bridge_correlations.csv", "H2O2 Stress"),
        ("PRG4_Baseline_bridge_correlations.csv", "PRG4 Baseline"),
        ("PRG4_Rescue_bridge_correlations.csv", "PRG4 Rescue")
    ]
    
    for fname, label in files_to_plot:
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            print(f"Plotting rankings for {label}...")
            plot_top_correlations(fpath, label)

    # 2. Scatter plots for key "mimics"
    # Example: SNIP1 for H2O2 Stress
    print("Plotting bridge evidence for SNIP1 (Stress Mimic)...")
    plot_bridge_evidence(
        os.path.join(BULK_DIR, "H2O2_DEGs.csv"),
        "H2O2_vs_CTRL.log2FoldChange",
        "SNIP1",
        "H2O2_Stress_SNIP1",
        None
    )
    
    # Example: ARL4D for PRG4 Rescue
    print("Plotting bridge evidence for ARL4D (PRG4 Rescue Proxy)...")
    plot_bridge_evidence(
        os.path.join(BULK_DIR, "PRG4_Rescue_DEGs.csv"),
        "H2O2PRG4_vs_H2O2.log2FoldChange",
        "ARL4D",
        "PRG4_Rescue_ARL4D",
        None
    )

    print(f"All visualizations saved to {OUTPUT_DIR}")
