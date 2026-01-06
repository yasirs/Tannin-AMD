import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import pearsonr, spearmanr

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
FILE_RPE1 = os.path.join(BASE_DIR, "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad")
FILE_K562 = os.path.join(BASE_DIR, "data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/concordance-analysis")

# ==========================================
# FUNCTIONS
# ==========================================

def parse_obs_name(obs_name):
    # Format: [ID]_[Target_Gene]_[Site]_[EnsemblID]
    # We want Target_Gene and EnsemblID
    parts = obs_name.split('_')
    if len(parts) >= 4:
        return parts[1], parts[-1] # Gene Symbol, Ensembl ID
    elif "non-targeting" in obs_name:
        return "non-targeting", "non-targeting"
    return None, None

def load_and_prep(file_path, name):
    print(f"Loading {name}...")
    ad = sc.read_h5ad(file_path)
    
    # Map perturbations to common ID (Ensembl Gene ID of the target)
    obs_ids = []
    valid_indices = []
    
    for i, obs in enumerate(ad.obs_names):
        sym, ens = parse_obs_name(obs)
        if ens and ens != "non-targeting":
            obs_ids.append(ens)
            valid_indices.append(i)
            
    # Subset to valid gene perturbations
    ad_sub = ad[valid_indices].copy()
    ad_sub.obs["target_ensembl"] = obs_ids
    
    # Average multiple guides per gene?
    # The files are "pseudo-bulk", so usually 1 row per guide.
    # To compare cell types, it's best to aggregate to gene-level LFC.
    # Let's average Z-scores for the same target gene.
    
    print(f"  Aggregating {len(ad_sub)} perturbations to gene level...")
    # Convert to df for groupby
    df = pd.DataFrame(ad_sub.X, index=ad_sub.obs["target_ensembl"], columns=ad_sub.var_names)
    df_mean = df.groupby(df.index).mean()
    
    return df_mean

# ==========================================
# MAIN
# ==========================================

def main():
    # 1. Load Data
    df_rpe = load_and_prep(FILE_RPE1, "RPE1")
    df_k562 = load_and_prep(FILE_K562, "K562")
    
    # 2. Intersect Data
    # Common Target Genes (Rows)
    common_targets = df_rpe.index.intersection(df_k562.index)
    print(f"\nCommon Targets (KD Genes): {len(common_targets)}")
    
    # Common Measured Genes (Cols)
    common_measured = df_rpe.columns.intersection(df_k562.columns)
    print(f"Common Measured Genes: {len(common_measured)}")
    
    # Align Dataframes
    mat_rpe = df_rpe.loc[common_targets, common_measured]
    mat_k562 = df_k562.loc[common_targets, common_measured]
    
    # 3. Calculate Correlations
    print("\nCalculating concordance...")
    corrs = []
    pvals = []
    
    for gene in common_targets:
        vec_rpe = mat_rpe.loc[gene].values
        vec_k562 = mat_k562.loc[gene].values
        
        # Pearson correlation of the profile vectors
        r, p = pearsonr(vec_rpe, vec_k562)
        corrs.append(r)
        pvals.append(p)
        
    res_df = pd.DataFrame({
        "target_ensembl_id": common_targets,
        "concordance_r": corrs,
        "pval": pvals
    })
    
    # Add Symbols (retrieve from one of the files or index lookup if possible, 
    # but we dropped obs columns. We can infer from index if we had a map, 
    # or just trust the ID for now. Let's try to map back using the raw logic if needed,
    # but for now ID is fine. Actually, let's just create a map from the original load).
    
    # Clean up results
    res_df = res_df.sort_values("concordance_r", ascending=False)
    
    # Save
    res_df.to_csv(os.path.join(OUTPUT_DIR, "gene_concordance.csv"), index=False)
    
    # Summary Stats
    print("\nConcordance Summary:")
    print(f"  Mean Correlation: {np.mean(corrs):.4f}")
    print(f"  Median Correlation: {np.median(corrs):.4f}")
    print(f"  Genes with r > 0.5: {(np.array(corrs) > 0.5).sum()}")
    print(f"  Genes with r > 0.3: {(np.array(corrs) > 0.3).sum()}")
    
    # Save Top 20 Concordant Genes
    print("\nTop 10 Concordant KDs:")
    print(res_df.head(10))
    
    # Save Bottom 10 (Cell-Type Specific?)
    print("\nBottom 10 Concordant KDs:")
    print(res_df.tail(10))

if __name__ == "__main__":
    main()
