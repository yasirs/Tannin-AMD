import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import linregress
import seaborn as sns

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
FILE_RPE1 = os.path.join(BASE_DIR, "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad")
FILE_K562 = os.path.join(BASE_DIR, "data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad")
CONCORDANCE_FILE = os.path.join(BASE_DIR, "results/concordance-analysis/gene_concordance.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/transfer-model")

# ==========================================
# FUNCTIONS
# ==========================================

def get_gene_profile(ad, gene_id):
    # Find observations for this gene KD
    # Obs names: [ID]_[Symbol]_[Site]_[Ensembl]
    mask = ad.obs_names.str.contains(gene_id)
    if mask.sum() == 0:
        return None
    
    # Average profile
    X_sub = ad[mask].X
    if hasattr(X_sub, "toarray"):
        X_sub = X_sub.toarray()
    
    return np.mean(X_sub, axis=0)

def main():
    print("Loading Concordance Data...")
    conc_df = pd.read_csv(CONCORDANCE_FILE)
    
    # Select Top 500 Concordant Genes for Model Building
    # We want to train the transfer function on genes we KNOW work similarly
    top_genes = conc_df.head(500)["target_ensembl_id"].tolist()
    
    print(f"Loading Perturb-seq Data for {len(top_genes)} training genes...")
    ad_rpe = sc.read_h5ad(FILE_RPE1)
    ad_k562 = sc.read_h5ad(FILE_K562)
    
    # Common measured genes (Features)
    common_feats = ad_rpe.var_names.intersection(ad_k562.var_names)
    print(f"  Common Measured Features: {len(common_feats)}")
    
    x_vals = [] # K562 LFCs
    y_vals = [] # RPE1 LFCs
    
    print("Extracting profiles...")
    for gene in top_genes:
        prof_rpe = get_gene_profile(ad_rpe, gene)
        prof_k562 = get_gene_profile(ad_k562, gene)
        
        if prof_rpe is not None and prof_k562 is not None:
            # Subset to common features
            # Need indices
            # Or make pandas series
            s_rpe = pd.Series(prof_rpe, index=ad_rpe.var_names)
            s_k562 = pd.Series(prof_k562, index=ad_k562.var_names)
            
            x_vals.extend(s_k562.loc[common_feats].values)
            y_vals.extend(s_rpe.loc[common_feats].values)
            
    print(f"Total data points: {len(x_vals)}")
    
    # Downsample for plotting if massive
    if len(x_vals) > 100000:
        idx = np.random.choice(len(x_vals), 50000, replace=False)
        x_plot = np.array(x_vals)[idx]
        y_plot = np.array(y_vals)[idx]
    else:
        x_plot = x_vals
        y_plot = y_vals
        
    # Regression
    slope, intercept, r_value, p_value, std_err = linregress(x_vals, y_vals)
    
    print("\nTransfer Model Stats:")
    print(f"  Slope (Beta): {slope:.4f}")
    print(f"  Intercept: {intercept:.4f}")
    print(f"  R-squared: {r_value**2:.4f}")
    
    # Save Stats
    with open(os.path.join(OUTPUT_DIR, "model_stats.txt"), "w") as f:
        f.write(f"Transfer Model (RPE = beta * K562 + alpha)\n")
        f.write(f"Based on Top 500 Concordant Genes\n")
        f.write(f"Slope (Beta): {slope:.4f}\n")
        f.write(f"Intercept: {intercept:.4f}\n")
        f.write(f"Correlation (r): {r_value:.4f}\n")
    
    # Plot
    plt.figure(figsize=(7, 7))
    # Hexbin for density
    plt.hexbin(x_plot, y_plot, gridsize=50, cmap='Blues', bins='log')
    # Regression line
    x_range = np.linspace(min(x_plot), max(x_plot), 100)
    plt.plot(x_range, intercept + slope * x_range, 'r--', label=f'Fit: y={slope:.2f}x')
    plt.plot(x_range, x_range, 'k:', alpha=0.5, label='Identity (y=x)')
    
    plt.title("Transfer Model: K562 vs RPE1\n(Top 500 Concordant Genes)")
    plt.xlabel("K562 LFC (Z-score)")
    plt.ylabel("RPE1 LFC (Z-score)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "transfer_scatter_density.png"))
    plt.savefig(os.path.join(OUTPUT_DIR, "transfer_scatter_density.pdf")) # For user

if __name__ == "__main__":
    main()
