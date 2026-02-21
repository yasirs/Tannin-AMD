
import pandas as pd
import scanpy as sc
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns

# Constants
OUT_DIR = "results/three-axes/perturbseq"
GENE_SETS_DIR = "results/three-axes/gene-sets"
DATASETS = {
    "RPE1": "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad",
    "K562_GWPS": "data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad"
}

# Known AMD Genes (KD Targets)
GWAS_GENES = [
    "CFH", "C3", "ARMS2", "HTRA1", "APOE", "TIMP3", "RPE65", "BEST1", "RLBP1", "TTR",
    "CFI", "C9", "LIPC", "CETP", "COL8A1", "COL10A1", "TNFRSF10A", "IER3", "SLC16A8", "B3GALTL",
    "CDKN1A", "CDKN2A" # Add controls
]

def analyze_gwas(name, h5ad_path):
    print(f"Analyzing GWAS genes in {name}...")
    try:
        adata = sc.read_h5ad(h5ad_path, backed='r')
    except:
        return

    # Find Knockdowns for GWAS genes
    # Parse obs_names
    target_map = {} # Gene -> [ObsName, ...]
    for x in adata.obs_names:
        parts = x.split('_')
        if len(parts) >= 2:
            gene = parts[1]
            if gene in GWAS_GENES:
                if gene not in target_map: target_map[gene] = []
                target_map[gene].append(x)
                
    found_genes = list(target_map.keys())
    print(f"  Found {len(found_genes)} GWAS genes in KDs")
    
    if not found_genes: return

    # Load Axis Genes
    axes_genes = {}
    for ax in ["Oxidative", "Inflammatory", "Senescence"]:
        path = os.path.join(GENE_SETS_DIR, f"{ax.lower()}_pro_disease_mapped.csv")
        if os.path.exists(path):
            df = pd.read_csv(path)
            axes_genes[ax] = df['ensembl_gene_id'].dropna().unique().tolist()

    # Calculate Axis Scores for these KDs
    # Matrix: (KDs, Genes)
    results = [] # Rows for heatmap
    
    # We need to slice per KD or batch.
    # Collect all needed obs_names
    all_obs = []
    for g in found_genes:
        all_obs.extend(target_map[g])
        
    # Indices
    obs_indices = [i for i, x in enumerate(adata.obs_names) if x in all_obs]
    
    # We also only need columns matching axis genes
    all_axis_genes = set()
    for g in axes_genes.values(): all_axis_genes.update(g)
    
    # Filter valid vars
    var_names = list(adata.var_names)
    var_map = {name: i for i, name in enumerate(var_names)}
    
    valid_vars = [g for g in all_axis_genes if g in var_map]
    valid_var_indices = [var_map[g] for g in valid_vars]
    
    if not valid_var_indices: return
    
    # Load Subset
    # adata[obs, var] - split for backed mode
    subset_obs = adata[obs_indices].to_memory()
    subset = subset_obs[:, valid_vars]
    
    mat = subset.X
    if hasattr(mat, 'toarray'): mat = mat.toarray()
    
    # Remap matrix columns to gene names for scoring
    # subset.var_names corresponds to valid_vars
    sub_var = list(subset.var_names)
    
    # Score
    # Iterate over rows (KDs)
    for i in range(mat.shape[0]):
        row_vec = mat[i]
        obs_name = subset.obs_names[i]
        
        # Identify gene target
        # obs_name -> gene
        parts = obs_name.split('_')
        gene_target = parts[1] if len(parts) >= 2 else "Unknown"
        
        scores = {}
        for ax, agenes in axes_genes.items():
            # Indices in subset
            # intersection of agenes and sub_var
            # Get indices of agenes in sub_var
            # Precompute mask?
            # fast way:
            indices = [j for j, v in enumerate(sub_var) if v in agenes]
            if indices:
                scores[ax] = np.mean(row_vec[indices])
            else:
                scores[ax] = 0 # or NaN
        
        scores['Gene_Target'] = gene_target
        scores['Perturbation'] = obs_name
        results.append(scores)
        
    res_df = pd.DataFrame(results)
    
    # Aggregate by Gene (mean score across guides)
    agg_df = res_df.groupby('Gene_Target')[['Oxidative', 'Inflammatory', 'Senescence']].mean()
    
    agg_df.to_csv(os.path.join(OUT_DIR, f"gwas_gene_axis_effects_{name}.csv"))
    
    # Plot Heatmap
    if not agg_df.empty:
        plt.figure(figsize=(6, 8))
        sns.heatmap(agg_df, cmap="RdBu_r", center=0, annot=True, fmt=".2f")
        plt.title(f"GWAS KD Effects on Axes ({name})")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, f"Fig_GWAS_gene_axis_effects_{name}.pdf"))
        plt.close()

if __name__ == "__main__":
    for k, v in DATASETS.items():
        analyze_gwas(k, v)
