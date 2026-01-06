import os
import pandas as pd
import numpy as np
import scipy.stats as stats

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
COUNTS_FILE = os.path.join(BASE_DIR, "data/external/GSE129964/GSE129964_serumCountsTable.tsv.gz")
OUR_SIG_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/RPE_gene pvals.xlsx")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/external-validation")

# ==========================================
# FUNCTIONS
# ==========================================

def load_and_process_counts():
    print("Loading GSE129964 counts...")
    df = pd.read_csv(COUNTS_FILE, sep='\t', index_col=0)
    
    # CPM Normalization
    print("Normalizing to CPM...")
    libsizes = df.sum(axis=0)
    cpm = df.div(libsizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)
    
    # Define Groups
    # C.Day0 (Control) vs S.Day9 (Late Stress)
    ctrl_cols = [c for c in df.columns if "C.Day0" in c]
    stress_cols = [c for c in df.columns if "S.Day9" in c]
    
    print(f"Control samples: {len(ctrl_cols)}, Stress samples: {len(stress_cols)}")
    
    # Differential Expression (Simple T-test for speed)
    # Or just LogFC of means
    ctrl_mean = log_cpm[ctrl_cols].mean(axis=1)
    stress_mean = log_cpm[stress_cols].mean(axis=1)
    logfc = stress_mean - ctrl_mean
    
    # T-test
    pvals = []
    for gene in df.index:
        c_vals = log_cpm.loc[gene, ctrl_cols]
        s_vals = log_cpm.loc[gene, stress_cols]
        # Check variance
        if np.var(c_vals) < 1e-6 and np.var(s_vals) < 1e-6:
            pvals.append(1.0)
        else:
            _, p = stats.ttest_ind(c_vals, s_vals)
            pvals.append(p)
            
    res_df = pd.DataFrame({
        "gene_symbol": df.index,
        "logfc_serum_stress": logfc,
        "pval_serum_stress": pvals
    }).set_index("gene_symbol")
    
    return res_df

def compare_signatures(ext_df, our_df):
    print("\nComparing signatures...")
    
    # Merge on symbol
    # Our data: H2O2_vs_CTRL.log2FoldChange
    merged = ext_df.join(our_df.set_index("hgnc_symbol"), how="inner")
    
    print(f"Common genes: {len(merged)}")
    
    # Filter for significant in BOTH (to compare valid signals)
    # Or just significant in one
    sig_ext = merged[merged["pval_serum_stress"] < 0.05]
    sig_our_h2o2 = merged[merged["H2O2_vs_CTRL.pvalue"] < 0.05]
    
    # 1. Correlation with H2O2 Stress
    # Do H2O2 and Serum Starvation induce similar changes?
    # Use genes significant in at least one
    union_sig = merged[(merged["pval_serum_stress"] < 0.01) | (merged["H2O2_vs_CTRL.pvalue"] < 0.01)]
    
    corr_h2o2 = union_sig["logfc_serum_stress"].corr(union_sig["H2O2_vs_CTRL.log2FoldChange"])
    print(f"Correlation (Serum Stress vs H2O2 Stress): {corr_h2o2:.3f} (n={len(union_sig)})")
    
    # 2. Correlation with PRG4 Rescue
    # Does PRG4 rescue the Serum Starvation signature?
    # We expect NEGATIVE correlation if PRG4 fixes the stress.
    corr_rescue = union_sig["logfc_serum_stress"].corr(union_sig["H2O2PRG4_vs_H2O2.log2FoldChange"])
    print(f"Correlation (Serum Stress vs PRG4 Rescue): {corr_rescue:.3f} (n={len(union_sig)})")
    
    # Save merged for plotting
    merged.to_csv(os.path.join(OUTPUT_DIR, "gse129964_comparison.csv"))

# ==========================================
# MAIN
# ==========================================

def main():
    ext_res = load_and_process_counts()
    
    print("Loading our bulk data...")
    our_res = pd.read_excel(OUR_SIG_FILE)
    
    compare_signatures(ext_res, our_res)

if __name__ == "__main__":
    main()
