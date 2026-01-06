import os
import pandas as pd
import numpy as np
import scanpy as sc

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
RAW_PERTURB_FILES = {
    "RPE1": os.path.join(BASE_DIR, "data/external/perturbseq/rpe1_raw_bulk_01.h5ad"),
    "K562_Ess": os.path.join(BASE_DIR, "data/external/perturbseq/K562_essential_raw_bulk_01.h5ad"),
    "K562_GWPS": os.path.join(BASE_DIR, "data/external/perturbseq/K562_gwps_raw_bulk_01.h5ad")
}
OUTPUT_DIR = os.path.join(BASE_DIR, "results/baseline-expression")

MARKERS = {
    "RPE_Specific": ["BEST1", "RPE65", "RLBP1", "TTR"],
    "K562_Specific": ["HBE1", "HBG1", "HBG2", "CD34"],
    "Housekeeping": ["ACTB", "GAPDH", "TBP"]
}

# ==========================================
# FUNCTIONS
# ==========================================

def get_baseline_expression(file_path, name):
    print(f"Loading {name} raw data from {file_path}...")
    ad = sc.read_h5ad(file_path)
    
    nt_mask = ad.obs_names.str.contains("non-targeting", case=False)
    n_nt = nt_mask.sum()
    
    if n_nt > 0:
        X_nt = ad[nt_mask].X
        if hasattr(X_nt, "toarray"):
            X_nt = X_nt.toarray()
        # Use nanmean to handle any sparse issues/nans
        baseline = np.nanmean(X_nt, axis=0)
    else:
        X = ad.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        baseline = np.nanmedian(X, axis=0)
        
    res_df = pd.DataFrame({
        "ensembl_gene_id": ad.var_names,
        "gene_name": ad.var["gene_name"] if "gene_name" in ad.var.columns else "N/A",
        f"{name}_baseline": baseline
    })
    
    return res_df

# ==========================================
# MAIN
# ==========================================

def main():
    dfs = []
    for name, path in RAW_PERTURB_FILES.items():
        dfs.append(get_baseline_expression(path, name))
        
    # Build a master mapping of ID to Name
    id_to_name = {}
    for df in dfs:
        for _, row in df.iterrows():
            if row["gene_name"] != "N/A":
                id_to_name[row["ensembl_gene_id"]] = row["gene_name"]

    # Merge
    merged = dfs[0][["ensembl_gene_id", "RPE1_baseline"]]
    for df in dfs[1:]:
        col_name = [c for c in df.columns if "_baseline" in c][0]
        merged = merged.merge(df[["ensembl_gene_id", col_name]], on="ensembl_gene_id", how="outer")
    
    merged["gene_name"] = merged["ensembl_gene_id"].map(id_to_name).fillna("N/A")
    
    # Threshold for "Expressed"
    # Given RPE1 ACTB ~50 and TBP ~0.1, let's use 0.05 as a baseline.
    THRESH = 0.05
    
    # Save tables
    merged.to_csv(os.path.join(OUTPUT_DIR, "all_baseline_expression.csv"), index=False)
    
    # Marker results
    marker_results = []
    print("\nMarker Check:")
    for cat, genes in MARKERS.items():
        for g in genes:
            match = merged[merged["gene_name"] == g]
            if not match.empty:
                rpe_val = match["RPE1_baseline"].values[0]
                k562_val = match["K562_GWPS_baseline"].values[0]
                print(f"  {g:8}: RPE1={rpe_val:.4f}, K562={k562_val:.4f}")
                marker_results.append({"category": cat, "gene": g, "RPE1": rpe_val, "K562": k562_val})
            else:
                # Try partial match or case-insensitive?
                pass
    
    pd.DataFrame(marker_results).to_csv(os.path.join(OUTPUT_DIR, "marker_baseline_check.csv"), index=False)
    
    # Summary
    summary = pd.DataFrame({
        "Dataset": ["RPE1", "K562_Ess", "K562_GWPS"],
        "Expressed_GT_0.05": [
            (merged["RPE1_baseline"] > THRESH).sum(),
            (merged["K562_Ess_baseline"] > THRESH).sum(),
            (merged["K562_GWPS_baseline"] > THRESH).sum()
        ]
    })
    summary.to_csv(os.path.join(OUTPUT_DIR, "baseline_expression_summary.csv"), index=False)
    print("\n", summary)

if __name__ == "__main__":
    main()