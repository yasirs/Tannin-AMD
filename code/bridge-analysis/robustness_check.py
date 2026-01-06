import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
BULK_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/RPE_gene pvals.xlsx")
PERTURB_FILE = os.path.join(BASE_DIR, "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/robustness-analysis")

CONTRASTS = {
    "H2O2_Stress": {
        "col_fc": "H2O2_vs_CTRL.log2FoldChange",
        "col_pval": "H2O2_vs_CTRL.pvalue"
    },
    "PRG4_Baseline": {
        "col_fc": "PRG4_vs_CTRL.log2FoldChange",
        "col_pval": "PRG4_vs_CTRL.pvalue"
    },
    "PRG4_Rescue": {
        "col_fc": "H2O2PRG4_vs_H2O2.log2FoldChange",
        "col_pval": "H2O2PRG4_vs_H2O2.pvalue"
    }
}

THRESHOLDS = [
    {"type": "nominal_p", "value": 0.05, "label": "p < 0.05"},
    {"type": "nominal_p", "value": 0.01, "label": "p < 0.01"},
    {"type": "fdr", "value": 0.1, "label": "FDR < 0.1"},
    {"type": "fdr", "value": 0.05, "label": "FDR < 0.05"},
    {"type": "top_n", "value": 1000, "label": "Top 1000"},
    {"type": "top_n", "value": 2000, "label": "Top 2000"},
    {"type": "top_n", "value": 3000, "label": "Top 3000"},
]

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================================
# FUNCTIONS
# ==========================================

def load_bulk_data():
    print(f"Loading bulk data from {BULK_FILE}...")
    df = pd.read_excel(BULK_FILE)
    
    # Calculate FDR for each contrast if not present (it is not present based on checks)
    for name, cols in CONTRASTS.items():
        pval_col = cols["col_pval"]
        # Drop NaNs for calculation
        mask = df[pval_col].notna()
        pvals = df.loc[mask, pval_col].values
        
        # BH Correction
        _, adj_pvals, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        
        # Assign back
        padj_col = pval_col.replace("pvalue", "padj")
        df.loc[mask, padj_col] = adj_pvals
        CONTRASTS[name]["col_padj"] = padj_col
        
    return df

def load_perturb_data():
    print(f"Loading Perturb-seq data from {PERTURB_FILE}...")
    return sc.read_h5ad(PERTURB_FILE)

def get_signature(df, contrast_name, threshold_cfg):
    cols = CONTRASTS[contrast_name]
    fc_col = cols["col_fc"]
    pval_col = cols["col_pval"]
    padj_col = cols["col_padj"]
    
    # Base filter: remove NaNs in FC or Pval
    sub = df.dropna(subset=[fc_col, pval_col])
    
    if threshold_cfg["type"] == "nominal_p":
        mask = sub[pval_col] < threshold_cfg["value"]
        sig = sub[mask]
    elif threshold_cfg["type"] == "fdr":
        mask = sub[padj_col] < threshold_cfg["value"]
        sig = sub[mask]
    elif threshold_cfg["type"] == "top_n":
        # Sort by absolute log2FC descending
        sub["abs_fc"] = sub[fc_col].abs()
        sig = sub.sort_values("abs_fc", ascending=False).head(threshold_cfg["value"])
    
    return sig.set_index("ensembl_gene_id")[fc_col].to_dict()

def compute_bridge_correlation(sig_dict, ad):
    # Intersect genes
    common_genes = list(set(sig_dict.keys()) & set(ad.var_names))
    if len(common_genes) < 10:
        return None, 0
    
    # Extract vectors
    sig_vec = np.array([sig_dict[g] for g in common_genes])
    
    # Subset AnnData
    ad_sub = ad[:, common_genes]
    X = ad_sub.X
    if hasattr(X, "toarray"):
        X = X.toarray()
        
    # Compute Spearman correlation for all perturbations
    # X shape: (n_perturbations, n_genes)
    # sig_vec shape: (n_genes,)
    
    # We can use apply_along_axis or a loop. Loop is safer for memory if massive, but X is 2k x 8k, it fits.
    # Scipy spearmanr can take matrix, but computes pairwise. We want 1-to-many.
    
    # Let's use a loop or vectorized approach if possible.
    # Manual rank correlation to speed up?
    # For now, simple loop is robust.
    
    corrs = []
    perturbations = ad_sub.obs_names.tolist()
    
    # Pre-rank signature for speed? No, spearmanr does it.
    
    # Optimization: Use pandas corrwith if we construct a DF?
    # Or just loop. 2000 perturbations is fast.
    
    for i in range(X.shape[0]):
        r, _ = spearmanr(sig_vec, X[i, :])
        corrs.append(r)
        
    res_df = pd.DataFrame({
        "perturbation": perturbations,
        "spearman_rho": corrs
    })
    
    # Parse target gene
    res_df["target_gene"] = res_df["perturbation"].apply(lambda x: x.split("_")[1] if "_" in x else x)
    
    return res_df.sort_values("spearman_rho", ascending=False), len(common_genes)

def jaccard_similarity(list1, list2):
    s1 = set(list1)
    s2 = set(list2)
    if len(s1) == 0 or len(s2) == 0:
        return 0.0
    return len(s1.intersection(s2)) / len(s1.union(s2))

# ==========================================
# MAIN EXECUTION
# ==========================================

def main():
    bulk_df = load_bulk_data()
    ad = load_perturb_data()
    
    summary_stats = []
    top_hits_dict = {} # Key: (contrast, threshold_label) -> list of top 50 KDs
    
    for contrast_name in CONTRASTS.keys():
        print(f"\nProcessing Contrast: {contrast_name}")
        
        for thresh in THRESHOLDS:
            label = thresh["label"]
            print(f"  Threshold: {label}")
            
            # 1. Get Signature
            sig_dict = get_signature(bulk_df, contrast_name, thresh)
            n_degs = len(sig_dict)
            
            # 2. Run Correlation
            res_df, n_overlap = compute_bridge_correlation(sig_dict, ad)
            
            top_corr = 0
            top_50_kds = []
            
            if res_df is not None:
                top_corr = res_df.iloc[0]["spearman_rho"]
                top_50_kds = res_df.head(50)["target_gene"].tolist()
                
                # Save full results for this run? Maybe just top 50 to save space/time
                # Or output strictly what is requested.
                # Let's save the top 100 for each combo to a CSV for inspection
                out_name = f"{contrast_name}_{label.replace(' ', '_').replace('<', 'lt')}_top100.csv"
                res_df.head(100).to_csv(os.path.join(OUTPUT_DIR, out_name), index=False)
            
            # 3. Store Stats
            summary_stats.append({
                "contrast": contrast_name,
                "threshold": label,
                "n_degs": n_degs,
                "n_overlap": n_overlap,
                "top_correlation": top_corr,
                "top_hit": top_50_kds[0] if top_50_kds else "N/A"
            })
            
            top_hits_dict[(contrast_name, label)] = top_50_kds

    # ==========================================
    # STABILITY ANALYSIS
    # ==========================================
    print("\nCalculating stability metrics...")
    stability_results = []
    
    for contrast_name in CONTRASTS.keys():
        # Compare all pairs of thresholds within this contrast
        t_labels = [t["label"] for t in THRESHOLDS]
        
        for i in range(len(t_labels)):
            for j in range(i+1, len(t_labels)):
                t1 = t_labels[i]
                t2 = t_labels[j]
                
                hits1 = top_hits_dict.get((contrast_name, t1), [])
                hits2 = top_hits_dict.get((contrast_name, t2), [])
                
                jaccard = jaccard_similarity(hits1, hits2)
                
                stability_results.append({
                    "contrast": contrast_name,
                    "threshold_1": t1,
                    "threshold_2": t2,
                    "jaccard_top50": jaccard
                })

    # ==========================================
    # OUTPUTS
    # ==========================================
    pd.DataFrame(summary_stats).to_csv(os.path.join(OUTPUT_DIR, "threshold_comparison.csv"), index=False)
    pd.DataFrame(stability_results).to_csv(os.path.join(OUTPUT_DIR, "stability_metrics.csv"), index=False)
    
    print("\nAnalysis Complete.")
    print(f"Results saved to {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
