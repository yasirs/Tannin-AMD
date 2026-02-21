
import pandas as pd
import scanpy as sc
import numpy as np
import os
import glob

# Constants
GENE_SETS_DIR = "results/three-axes/gene-sets"
OUT_DIR = "results/three-axes/perturbseq"
os.makedirs(OUT_DIR, exist_ok=True)

# Datasets
DATASETS = {
    "RPE1": "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad",
    "K562_Essential": "data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad",
    "K562_GWPS": "data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad"
}

def load_axis_genes():
    axes = {}
    for ax in ["Oxidative", "Inflammatory", "Senescence"]:
        path = os.path.join(GENE_SETS_DIR, f"{ax.lower()}_pro_disease_mapped.csv")
        if os.path.exists(path):
            df = pd.read_csv(path)
            # Use Ensembl IDs for Perturb-seq
            axes[ax] = df['ensembl_gene_id'].dropna().unique().tolist()
        else:
            print(f"Warning: {path} not found")
    return axes

def analyze_dataset(name, h5ad_path, axes_genes):
    print(f"Analyzing {name}...")
    try:
        adata = sc.read_h5ad(h5ad_path, backed='r')
    except Exception as e:
        print(f"  Error loading {name}: {e}")
        return
    
    # Check if X is sparse or dense
    # Backed mode with sparse X needs care.
    # Scores = Mean of Z-scores of axis genes.
    # Data is likely Z-scores (mean ~0).
    
    results = {}
    
    # Parse gene targets from obs_names
    # Format: ID_GENE_...
    # We assume Gene Symbol is 2nd token.
    obs_names = adata.obs_names
    gene_targets = []
    for x in obs_names:
        parts = x.split('_')
        if len(parts) >= 2:
            gene_targets.append(parts[1])
        else:
            gene_targets.append(x) # Fallback
            
    # Calculate scores per axis
    var_names = set(adata.var_names)
    
    for ax, genes in axes_genes.items():
        # Filter genes present in dataset
        valid_genes = [g for g in genes if g in var_names]
        if not valid_genes:
            print(f"  No genes found for {ax} in {name}")
            continue
            
        # Get indices
        # Reading specific columns from backed h5ad is slow if random access?
        # Better to read entire matrix if fits in memory?
        # K562 GWPS is large (2.2K x 11K?). Wait, 11K KDs x Genome?
        # Shape: (Obs, Vars).
        # We need to average across Genes (Vars).
        # slice: adata[:, valid_genes]
        
        # We assume dataset fits in memory (L40 GPU machine, 48GB GPU, likely enough RAM).
        # Actually RPE1 is small. K562 GWPS is larger.
        # I'll convert to memory for speed if possible.
        # But safest is chunked or just slicing.
        
        # Slice
        subset = adata[:, valid_genes]
        # X might be sparse CSR.
        # Mean across axis 1 (genes)
        # We need to compute mean(subset.X, axis=1)
        
        if hasattr(subset.X, 'toarray'): # Sparse
            # Only materialize the subset
            data_mat = subset.X # This might be h5py dataset
            # If backed, subset.X is not in memory yet?
            # scanpy backed mode returns a view.
            # Convert to memory
            data_mat = subset.to_memory().X
        else:
            data_mat = subset.X # Dense or already loaded?
            
        # If still on disk, to_memory() brings it in.
        # If it's backed, subset is an AnnData view on disk. .to_memory() loads it.
        
        if hasattr(data_mat, "toarray"):
            data_mat = data_mat.toarray()
            
        scores = np.mean(data_mat, axis=1) # Mean score per perturbation
        
        if hasattr(scores, "A1"): # matrix
            scores = scores.A1
            
        # Create DF
        res_df = pd.DataFrame({
            "Perturbation": obs_names,
            "Gene_Target": gene_targets,
            "Axis_Score": scores
        })
        
        # Rank
        res_df = res_df.sort_values("Axis_Score", ascending=False)
        
        # Save rankings
        res_df.to_csv(os.path.join(OUT_DIR, f"axis_regulators_{name}_{ax}.csv"), index=False)
        results[ax] = res_df
        
        print(f"  Saved rankings for {ax} in {name}. Top regulator: {res_df.iloc[0]['Gene_Target']}")

    return results

def main():
    axes_genes = load_axis_genes()
    
    combine_list = []
    
    for name, path in DATASETS.items():
        res = analyze_dataset(name, path, axes_genes)
        # Store for concordance
        # (Implementing basic concordance analysis later if needed or here)
        
if __name__ == "__main__":
    main()
