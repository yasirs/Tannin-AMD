"""
Perturb-seq Data Loader

Load and cache Perturb-seq h5ad datasets for analysis.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import sys
import scanpy as sc

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))
from config import PERTURBSEQ_DATA, OUTPUT_DIRS, ensure_dirs

#%%
def load_perturbseq_data(dataset_name: str, use_cache: bool = True) -> dict:
    """
    Load Perturb-seq dataset and prepare for analysis.
    
    Parameters
    ----------
    dataset_name : str
        Name of dataset: "RPE1" or "K562_GWPS"
    use_cache : bool
        Whether to use cached data if available
        
    Returns
    -------
    dict with keys:
        - expression: np.ndarray (n_perturbations x n_genes)
        - gene_names: list of gene symbols
        - perturbation_ids: list of perturbation IDs
        - target_genes: list of target gene symbols for each perturbation
        - metadata: pd.DataFrame of perturbation metadata
    """
    ensure_dirs()
    
    # Check cache
    cache_file = OUTPUT_DIRS["cache"] / f"{dataset_name}_processed.pkl"
    
    if use_cache and cache_file.exists():
        print(f"Loading cached {dataset_name} data from {cache_file}")
        with open(cache_file, "rb") as f:
            return pickle.load(f)
    
    # Load from h5ad
    h5ad_path = PERTURBSEQ_DATA.get(dataset_name)
    if not h5ad_path or not h5ad_path.exists():
        raise FileNotFoundError(f"Perturb-seq data not found: {h5ad_path}")
    
    print(f"Loading {dataset_name} from {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"  Shape: {adata.shape} (perturbations x genes)")
    
    # Extract data
    # Expression matrix
    if hasattr(adata.X, 'toarray'):
        expression = adata.X.toarray()
    else:
        expression = np.array(adata.X)
    
    # Gene names from var - prioritize gene_name column
    if 'gene_name' in adata.var.columns:
        gene_names = adata.var['gene_name'].tolist()
        print(f"  Using 'gene_name' column from var")
    elif 'gene_symbol' in adata.var.columns:
        gene_names = adata.var['gene_symbol'].tolist()
    elif 'gene' in adata.var.columns:
        gene_names = adata.var['gene'].tolist()
    else:
        gene_names = adata.var_names.tolist()
    
    # Clean gene names - remove any None/NaN
    gene_names = [str(g) if g is not None else "" for g in gene_names]
    print(f"  Sample gene names: {gene_names[:5]}")
    
    # Perturbation IDs and target genes from obs
    perturbation_ids = adata.obs_names.tolist()
    
    # Extract target gene from perturbation ID
    # Common formats: "ID_GENESYMBOL_P1_ENSGID" (e.g., "10005_ZBTB4_P1_ENSG00000174282")
    if 'target_gene' in adata.obs.columns:
        target_genes = adata.obs['target_gene'].tolist()
    elif 'gene' in adata.obs.columns:
        target_genes_raw = adata.obs['gene'].tolist()
        # Check if these are gene symbols or numeric IDs
        if len(target_genes_raw) > 0 and str(target_genes_raw[0]).isdigit():
            # Numeric - need to parse from perturbation_ids
            target_genes = []
            for pid in perturbation_ids:
                # Format: ID_GENESYMBOL_P1_ENSGID or ID_GENESYMBOL_ENSGID
                parts = str(pid).split("_")
                if len(parts) >= 2:
                    # Second part is usually the gene symbol
                    potential_gene = parts[1]
                    # Skip if it looks like an ENSG ID or P1/P2 marker
                    if potential_gene.startswith("ENSG") or potential_gene.startswith("ENST"):
                        # Try third part
                        if len(parts) >= 3:
                            potential_gene = parts[2]
                    target_genes.append(potential_gene)
                else:
                    target_genes.append(parts[0])
            print(f"  Parsed gene symbols from perturbation IDs")
            print(f"  Sample target genes: {target_genes[:5]}")
        else:
            target_genes = target_genes_raw
    else:
        # Try to parse from index - common format is "ID_GENE_P1_ENSG" 
        target_genes = []
        for pid in perturbation_ids:
            parts = str(pid).split("_")
            if len(parts) >= 2:
                potential_gene = parts[1]
                if potential_gene.startswith("ENSG") or potential_gene.startswith("ENST"):
                    if len(parts) >= 3:
                        potential_gene = parts[2]
                target_genes.append(potential_gene)
            else:
                target_genes.append(parts[0])
        print(f"  Parsed gene symbols from perturbation IDs (no obs column)")
        print(f"  Sample target genes: {target_genes[:5]}")
    
    # Build result dict
    result = {
        "expression": expression,
        "gene_names": gene_names,
        "perturbation_ids": perturbation_ids,
        "target_genes": target_genes,
        "metadata": adata.obs.copy(),
        "n_perturbations": len(perturbation_ids),
        "n_genes": len(gene_names),
    }
    
    # Cache the processed data
    print(f"Caching processed data to {cache_file}...")
    with open(cache_file, "wb") as f:
        pickle.dump(result, f)
    
    return result

#%%
def get_gene_indices(gene_names: list, query_genes: list) -> tuple:
    """
    Get indices of query genes in the gene_names list.
    
    Returns
    -------
    tuple: (indices, matched_genes, missing_genes)
    """
    gene_name_set = {g.upper(): i for i, g in enumerate(gene_names)}
    
    indices = []
    matched = []
    missing = []
    
    for gene in query_genes:
        gene_upper = str(gene).upper().strip()
        if gene_upper in gene_name_set:
            indices.append(gene_name_set[gene_upper])
            matched.append(gene)
        else:
            missing.append(gene)
    
    return indices, matched, missing

#%%
def get_perturbation_indices(target_genes: list, query_genes: list) -> tuple:
    """
    Get indices of perturbations targeting query genes.
    
    Returns
    -------
    tuple: (indices, matched_genes)
    """
    # Create case-insensitive lookup
    target_set = {}
    for i, g in enumerate(target_genes):
        g_upper = str(g).upper().strip()
        if g_upper not in target_set:
            target_set[g_upper] = []
        target_set[g_upper].append(i)
    
    indices = []
    matched = []
    
    for gene in query_genes:
        gene_upper = str(gene).upper().strip()
        if gene_upper in target_set:
            indices.extend(target_set[gene_upper])
            matched.append(gene)
    
    return indices, matched

#%%
def compute_gene_coverage(data: dict, gene_set: list, gene_set_name: str = "") -> dict:
    """
    Compute coverage of a gene set in the Perturb-seq data.
    
    Returns dict with coverage statistics.
    """
    # Gene coverage (how many genes from set are in expression matrix)
    gene_indices, matched_genes, missing_genes = get_gene_indices(
        data["gene_names"], gene_set
    )
    
    gene_coverage = len(matched_genes) / len(gene_set) if gene_set else 0
    
    # Perturbation coverage (how many genes from set have knockdown data)
    pert_indices, matched_targets = get_perturbation_indices(
        data["target_genes"], gene_set
    )
    
    pert_coverage = len(matched_targets) / len(gene_set) if gene_set else 0
    
    result = {
        "gene_set_name": gene_set_name,
        "gene_set_size": len(gene_set),
        "genes_in_expression": len(matched_genes),
        "gene_coverage_pct": gene_coverage * 100,
        "perturbations_available": len(matched_targets),
        "pert_coverage_pct": pert_coverage * 100,
    }
    
    return result

#%%
def validate_data(data: dict) -> bool:
    """Validate loaded Perturb-seq data."""
    print("\nData Validation:")
    
    all_valid = True
    
    # Check dimensions
    n_pert, n_genes = data["expression"].shape
    print(f"  Expression matrix: {n_pert} x {n_genes}")
    
    # Check gene names
    n_gene_names = len(data["gene_names"])
    if n_gene_names != n_genes:
        print(f"  ⚠ Gene name count mismatch: {n_gene_names} vs {n_genes}")
        all_valid = False
    else:
        print(f"  ✓ Gene names: {n_gene_names}")
    
    # Check perturbation IDs
    n_pert_ids = len(data["perturbation_ids"])
    if n_pert_ids != n_pert:
        print(f"  ⚠ Perturbation ID count mismatch: {n_pert_ids} vs {n_pert}")
        all_valid = False
    else:
        print(f"  ✓ Perturbation IDs: {n_pert_ids}")
    
    # Check target genes
    n_targets = len(data["target_genes"])
    unique_targets = len(set(data["target_genes"]))
    print(f"  ✓ Target genes: {n_targets} total, {unique_targets} unique")
    
    # Check for NaN values
    nan_count = np.isnan(data["expression"]).sum()
    if nan_count > 0:
        print(f"  ⚠ NaN values in expression: {nan_count}")
    else:
        print(f"  ✓ No NaN values")
    
    return all_valid

#%%
def main():
    """Load and validate both Perturb-seq datasets."""
    print("="*60)
    print("PERTURB-SEQ DATA LOADING")
    print("="*60)
    
    for dataset_name in ["RPE1", "K562_GWPS"]:
        print(f"\n{'='*40}")
        print(f"Loading {dataset_name}...")
        print("="*40)
        
        try:
            data = load_perturbseq_data(dataset_name, use_cache=False)
            validate_data(data)
        except Exception as e:
            print(f"ERROR loading {dataset_name}: {e}")
            continue
    
    print("\n" + "="*60)
    print("Data loading complete!")

#%%
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate", action="store_true", help="Validate loaded data")
    parser.add_argument("--dataset", type=str, default=None, help="Specific dataset to load")
    args = parser.parse_args()
    
    if args.dataset:
        data = load_perturbseq_data(args.dataset, use_cache=True)
        validate_data(data)
    else:
        main()
