"""
Export Perturb-seq data for fgsea (R)

Exports expression matrix and metadata in formats that fgsea can read.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import argparse

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))

from config import OUTPUT_DIRS, ensure_dirs
from perturbseq_loader import load_perturbseq_data

#%%
def export_for_fgsea(dataset_name: str):
    """
    Export Perturb-seq data for fgsea.
    
    Creates:
    - {dataset}_expression.csv.gz : perturbation_id, gene1, gene2, ... (expression values)
    - {dataset}_metadata.csv : perturbation_id, target_gene
    """
    print(f"Exporting {dataset_name} for fgsea...")
    
    ensure_dirs()
    
    # Load data
    data = load_perturbseq_data(dataset_name, use_cache=True)
    
    expression = data["expression"]
    gene_names = data["gene_names"]
    perturbation_ids = data["perturbation_ids"]
    target_genes = data["target_genes"]
    
    print(f"  Expression: {expression.shape}")
    print(f"  Genes: {len(gene_names)}")
    print(f"  Perturbations: {len(perturbation_ids)}")
    
    # Create expression DataFrame
    expr_df = pd.DataFrame(
        expression,
        index=perturbation_ids,
        columns=gene_names
    )
    expr_df.index.name = "perturbation_id"
    expr_df = expr_df.reset_index()
    
    # Create metadata DataFrame
    meta_df = pd.DataFrame({
        "perturbation_id": perturbation_ids,
        "target_gene": target_genes,
    })
    
    # Save
    cache_dir = OUTPUT_DIRS["cache"]
    
    expr_file = cache_dir / f"{dataset_name}_expression.csv.gz"
    meta_file = cache_dir / f"{dataset_name}_metadata.csv"
    
    print(f"  Saving expression to {expr_file}...")
    expr_df.to_csv(expr_file, index=False, compression="gzip")
    
    print(f"  Saving metadata to {meta_file}...")
    meta_df.to_csv(meta_file, index=False)
    
    print(f"  Export complete!")
    print(f"    Expression: {expr_file.stat().st_size / 1e6:.1f} MB")
    print(f"    Metadata: {meta_file.stat().st_size / 1e3:.1f} KB")

#%%
def main():
    parser = argparse.ArgumentParser(description="Export Perturb-seq data for fgsea")
    parser.add_argument("--dataset", type=str, default=None, 
                       help="Dataset to export (RPE1 or K562_GWPS). If not specified, exports both.")
    args = parser.parse_args()
    
    print("="*60)
    print("EXPORTING PERTURB-SEQ DATA FOR fgsea")
    print("="*60)
    
    if args.dataset:
        export_for_fgsea(args.dataset)
    else:
        for dataset in ["RPE1", "K562_GWPS"]:
            print()
            export_for_fgsea(dataset)
    
    print("\n" + "="*60)
    print("Export complete!")
    print("Now run: Rscript run_fgsea_perturbseq.R --dataset [DATASET]")
    print("="*60)

#%%
if __name__ == "__main__":
    main()
