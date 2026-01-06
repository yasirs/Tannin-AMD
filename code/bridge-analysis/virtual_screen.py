import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
BULK_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/RPE_gene pvals.xlsx")
PERTURB_FILE = os.path.join(BASE_DIR, "data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad")
RPE_EXPRESSED_FILE = os.path.join(BASE_DIR, "results/baseline-expression/expressed_genes_RPE1.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/virtual-screen")

# Signature Definition
SIG_COL_FC = "H2O2PRG4_vs_H2O2.log2FoldChange"
SIG_COL_PVAL = "H2O2PRG4_vs_H2O2.pvalue"
FDR_THRESH = 0.05

# ==========================================
# FUNCTIONS
# ==========================================

def load_signature():
    print("Loading bulk signature...")
    df = pd.read_excel(BULK_FILE)
    
    # Calculate FDR
    mask = df[SIG_COL_PVAL].notna()
    pvals = df.loc[mask, SIG_COL_PVAL].values
    _, adj_pvals, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
    df.loc[mask, "padj"] = adj_pvals
    
    # Filter
    sig = df[df["padj"] < FDR_THRESH]
    print(f"  Signature size (FDR < {FDR_THRESH}): {len(sig)} genes")
    
    return sig.set_index("ensembl_gene_id")[SIG_COL_FC].to_dict()

def load_rpe_filter():
    print("Loading RPE1 expressed gene list...")
    df = pd.read_csv(RPE_EXPRESSED_FILE)
    # Return set of Ensembl IDs
    return set(df["ensembl_gene_id"].tolist())

def run_screen(sig_dict, rpe_filter):
    print(f"Loading K562 GWPS data from {PERTURB_FILE}...")
    ad = sc.read_h5ad(PERTURB_FILE)
    
    # 1. Intersect Genes (Signature vs Measured)
    common_genes = list(set(sig_dict.keys()) & set(ad.var_names))
    print(f"  Genes in signature and measured in K562: {len(common_genes)}")
    
    if len(common_genes) < 10:
        print("  Error: Too few common genes.")
        return None
    
    # 2. Prepare Data
    sig_vec = np.array([sig_dict[g] for g in common_genes])
    ad_sub = ad[:, common_genes]
    X = ad_sub.X
    if hasattr(X, "toarray"):
        X = X.toarray()
        
    # 3. Calculate Correlations
    print("  Calculating correlations for ~11,000 perturbations...")
    corrs = []
    
    # Optimization: Vectorized correlation?
    # Normalize inputs for fast cosine/correlation
    # (x - mean) / std
    
    # Doing manual loop for safety and simplicity with scipy
    for i in range(X.shape[0]):
        r, _ = spearmanr(sig_vec, X[i, :])
        corrs.append(r)
        
    res_df = pd.DataFrame({
        "perturbation": ad_sub.obs_names,
        "spearman_rho": corrs
    })
    
    # 4. Parse Metadata
    # Obs name format: [ID]_[Target_Gene]_[Site]_[EnsemblID]
    def parse(x):
        if "non-targeting" in x:
            return "non-targeting", "non-targeting"
        parts = x.split('_')
        if len(parts) >= 4:
            return parts[1], parts[-1]
        return x, x
        
    res_df[["gene_symbol", "ensembl_id"]] = res_df["perturbation"].apply(lambda x: pd.Series(parse(x)))
    
    # 5. Apply Filters
    # Filter out non-targeting
    res_df = res_df[res_df["gene_symbol"] != "non-targeting"]
    
    # Filter for RPE expression
    # Note: K562 GWPS has target genes (perturbations) that might NOT be in RPE.
    # We only want targets that exist in RPE, otherwise they aren't druggable targets in RPE.
    n_total = len(res_df)
    res_df = res_df[res_df["ensembl_id"].isin(rpe_filter)]
    print(f"  Filtered {n_total} -> {len(res_df)} perturbations based on RPE1 expression.")
    
    # 6. Rank
    res_df = res_df.sort_values("spearman_rho", ascending=False)
    
    return res_df

# ==========================================
# MAIN
# ==========================================

def main():
    sig = load_signature()
    rpe_filter = load_rpe_filter()
    
    results = run_screen(sig, rpe_filter)
    
    if results is not None:
        out_path = os.path.join(OUTPUT_DIR, "prg4_virtual_screen_results.csv")
        results.to_csv(out_path, index=False)
        print(f"\nSaved full results to {out_path}")
        
        print("\nTop 10 PRG4 Mimetics (Positive Correlation):")
        print(results.head(10)[["gene_symbol", "spearman_rho"]])
        
        print("\nTop 10 PRG4 Antagonists (Negative Correlation):")
        print(results.tail(10)[["gene_symbol", "spearman_rho"]])

if __name__ == "__main__":
    main()
