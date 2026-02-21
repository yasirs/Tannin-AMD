
import pandas as pd
import scanpy as sc
import os
import glob

# Constants
GENE_SETS_DIR = "results/three-axes/gene-sets"
BULK_FILE = "data/RPE_cells/code/RPE_gene pvals.xlsx"
PERTURB_RPE1 = "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad"
PERTURB_K562_ESS = "data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad"
PERTURB_K562_GWPS = "data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad"

def load_features(h5ad_path):
    try:
        adata = sc.read_h5ad(h5ad_path, backed='r')
        return set(adata.var_names)
    except Exception as e:
        print(f"Error loading {h5ad_path}: {e}")
        return set()

def main():
    print("Loading Master Mapping (Bulk RNA-seq)...")
    bulk_df = pd.read_excel(BULK_FILE)
    # Ensure columns exist
    if 'hgnc_symbol' not in bulk_df.columns or 'ensembl_gene_id' not in bulk_df.columns:
        print("Error: Required columns not found in bulk file.")
        return

    # Create mapping dicts
    sym2ens = dict(zip(bulk_df['hgnc_symbol'], bulk_df['ensembl_gene_id']))
    ens2sym = dict(zip(bulk_df['ensembl_gene_id'], bulk_df['hgnc_symbol']))
    
    valid_ensembl_bulk = set(bulk_df['ensembl_gene_id'].dropna())
    valid_symbol_bulk = set(bulk_df['hgnc_symbol'].dropna())
    
    print(f"Bulk Features: {len(valid_ensembl_bulk)} Ensembl IDs")

    print("\nLoading Perturb-seq Features...")
    feat_rpe1 = load_features(PERTURB_RPE1)
    feat_k562_ess = load_features(PERTURB_K562_ESS)
    feat_k562_gwps = load_features(PERTURB_K562_GWPS)
    
    print(f"RPE1: {len(feat_rpe1)}")
    print(f"K562 Ess: {len(feat_k562_ess)}")
    print(f"K562 GWPS: {len(feat_k562_gwps)}")

    # Process Gene Sets
    summary_data = []
    
    gene_set_files = glob.glob(os.path.join(GENE_SETS_DIR, "*_disease.csv")) + \
                     glob.glob(os.path.join(GENE_SETS_DIR, "known_amd_genes.csv"))
    
    print("\nProcessing Gene Sets...")
    for file_path in sorted(gene_set_files):
        if "filtered" in file_path: continue
        
        filename = os.path.basename(file_path)
        gs_name = filename.replace(".csv", "")
        
        df = pd.read_csv(file_path)
        if 'gene_symbol' not in df.columns:
            print(f"Skipping {filename}: no gene_symbol column")
            continue
            
        genes = df['gene_symbol'].dropna().unique()
        total_genes = len(genes)
        
        # Map to Ensembl
        mapped_genes = []
        for g in genes:
            ens = sym2ens.get(g)
            if ens:
                mapped_genes.append({'gene_symbol': g, 'ensembl_gene_id': ens})
                
        mapped_df = pd.DataFrame(mapped_genes)
        
        if mapped_df.empty:
            print(f"Warning: No genes mapped for {filename}")
            in_bulk = 0
            in_rpe1 = 0
            in_k562_ess = 0
            in_k562_gwps = 0
        else:
            # Check coverage
            ensembl_ids = mapped_df['ensembl_gene_id'].tolist()
            in_bulk = len(set(ensembl_ids) & valid_ensembl_bulk)
            in_rpe1 = len(set(ensembl_ids) & feat_rpe1)
            in_k562_ess = len(set(ensembl_ids) & feat_k562_ess)
            in_k562_gwps = len(set(ensembl_ids) & feat_k562_gwps)
            
            # Save mapped file
            out_path = file_path.replace(".csv", "_mapped.csv")
            mapped_df.to_csv(out_path, index=False)
            
        summary_data.append({
            'GeneSet': gs_name,
            'Total_Input': total_genes,
            'Mapped_to_Ensembl': len(mapped_df),
            'In_Bulk_RNA': in_bulk,
            'In_Perturb_RPE1': in_rpe1,
            'In_Perturb_K562_Ess': in_k562_ess,
            'In_Perturb_K562_GWPS': in_k562_gwps
        })
        print(f"  {gs_name}: {in_bulk}/{total_genes} in bulk")

    # Save summary
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(GENE_SETS_DIR, "coverage_summary.csv"), index=False)
    print("\nCoverage Summary Saved.")
    print(summary_df)

if __name__ == "__main__":
    main()
