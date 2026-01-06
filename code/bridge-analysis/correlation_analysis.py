import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import spearmanr

# Paths
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
BULK_DIR = os.path.join(BASE_DIR, "data/RPE_cells/code")
PERTURB_FILE = os.path.join(BASE_DIR, "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/bridge-results")

# Load Bulk DEGs
print("Loading bulk DEG signatures...")
df_h2o2 = pd.read_csv(os.path.join(BULK_DIR, "H2O2_DEGs.csv"))
df_rescue = pd.read_csv(os.path.join(BULK_DIR, "PRG4_Rescue_DEGs.csv"))
df_prg4_ctrl = pd.read_csv(os.path.join(BULK_DIR, "PRG4_Baseline_DEGs.csv"))

# Create a clean signature dictionary {ensembl_id: log2fc}
sig_h2o2 = df_h2o2.set_index("ensembl_gene_id")["H2O2_vs_CTRL.log2FoldChange"].dropna().to_dict()
sig_rescue = df_rescue.set_index("ensembl_gene_id")["H2O2PRG4_vs_H2O2.log2FoldChange"].dropna().to_dict()
sig_prg4_ctrl = df_prg4_ctrl.set_index("ensembl_gene_id")["PRG4_vs_CTRL.log2FoldChange"].dropna().to_dict()

# Load Perturb-seq
print(f"Loading Perturb-seq data from {PERTURB_FILE}...")
ad = sc.read_h5ad(PERTURB_FILE)

# Harmonize genes
# ad.var index is ensembl_id
common_genes_h2o2 = list(set(sig_h2o2.keys()) & set(ad.var_names))
common_genes_rescue = list(set(sig_rescue.keys()) & set(ad.var_names))
common_genes_prg4_ctrl = list(set(sig_prg4_ctrl.keys()) & set(ad.var_names))

print(f"Common genes for H2O2 Signature: {len(common_genes_h2o2)}")
print(f"Common genes for PRG4 Rescue Signature: {len(common_genes_rescue)}")
print(f"Common genes for PRG4 Baseline Signature: {len(common_genes_prg4_ctrl)}")

def run_correlation(sig_dict, common_genes, label):
    if not common_genes:
        print(f"No common genes found for {label}. Skipping.")
        return None
    print(f"Running correlations for {label}...")
    
    # Extract signature values in specific order
    sig_vec = np.array([sig_dict[g] for g in common_genes])
    
    # Subset AnnData to common genes
    ad_sub = ad[:, common_genes]
    
    # Get expression matrix (dense for easier correlation)
    X = ad_sub.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    
    # Perturb-seq profiles (KDs) are in rows (obs)
    perturbations = ad_sub.obs_names.tolist()
    
    corrs = []
    pvals = []
    
    for i in range(len(perturbations)):
        kd_vec = X[i, :]
        r, p = spearmanr(sig_vec, kd_vec)
        corrs.append(r)
        pvals.append(p)
    
    res_df = pd.DataFrame({
        "perturbation": perturbations,
        "spearman_rho": corrs,
        "pvalue": pvals
    })
    
    # Add target gene name if available in index or obs
    res_df["target_gene"] = res_df["perturbation"].apply(lambda x: x.split("_")[1] if "_" in x else x)
    
    res_df = res_df.sort_values("spearman_rho", ascending=False)
    
    out_path = os.path.join(OUTPUT_DIR, f"{label}_bridge_correlations.csv")
    res_df.to_csv(out_path, index=False)
    print(f"Saved results to {out_path}")
    return res_df

# Running Analysis
if len(common_genes_h2o2) > 0:
    res_h2o2 = run_correlation(sig_h2o2, common_genes_h2o2, "H2O2_Stress")
    print("\nTop 5 AMD-Mimics (Correlated with H2O2 stress):")
    print(res_h2o2.head(5)[["target_gene", "spearman_rho"]])

if len(common_genes_rescue) > 0:
    res_rescue = run_correlation(sig_rescue, common_genes_rescue, "PRG4_Rescue")
    print("\nTop 5 PRG4-Proxies (Correlated with PRG4 rescue):")
    print(res_rescue.head(5)[["target_gene", "spearman_rho"]])

if len(common_genes_prg4_ctrl) > 0:
    res_prg4_ctrl = run_correlation(sig_prg4_ctrl, common_genes_prg4_ctrl, "PRG4_Baseline")
    print("\nTop 5 PRG4-Baseline Mimics (Correlated with PRG4 vs CTRL):")
    print(res_prg4_ctrl.head(5)[["target_gene", "spearman_rho"]])
