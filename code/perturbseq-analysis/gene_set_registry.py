"""
Gene Set Registry

Load, extract, and manage all 13 gene sets for the Perturb-seq analysis.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))
from config import (
    GENE_SET_SOURCES, OUTPUT_DIRS, ANALYSIS_PARAMS, 
    GENE_SET_CATEGORIES, ensure_dirs
)

#%%
def load_age_gene_sets() -> dict:
    """Load Age-UP and Age-DOWN gene sets from GSE29801 DE results."""
    age_de = pd.read_csv(GENE_SET_SOURCES["age_de"])
    
    # Extract gene symbols and directions
    age_up = age_de[age_de["direction"] == "UP"]["gene"].dropna().unique().tolist()
    age_down = age_de[age_de["direction"] == "DOWN"]["gene"].dropna().unique().tolist()
    
    return {
        "Age-UP": age_up,
        "Age-DOWN": age_down,
    }

#%%
def load_are_gene_set() -> list:
    """Load ARE/NRF2 gene set."""
    are_df = pd.read_csv(GENE_SET_SOURCES["are"])
    return are_df["gene"].dropna().unique().tolist()

#%%
def extract_amd_gene_sets() -> dict:
    """Extract AMD gene sets from GSE135092 RPE DE results."""
    p_threshold = ANALYSIS_PARAMS["amd_pvalue_threshold"]
    gene_sets = {}
    
    for tissue, source_key in [("Macula", "amd_macula"), ("nonMacula", "amd_nonmacula")]:
        # Load DE results
        de_file = GENE_SET_SOURCES[source_key]
        if not de_file.exists():
            print(f"WARNING: {de_file} not found")
            continue
            
        de_df = pd.read_csv(de_file)
        
        # Identify the correct column names
        # P-value column
        if "P.Value" in de_df.columns:
            pval_col = "P.Value"
        elif "pvalue" in de_df.columns:
            pval_col = "pvalue"
        else:
            pval_col = [c for c in de_df.columns if "pval" in c.lower()][0]
        
        # LogFC column
        if "logFC" in de_df.columns:
            lfc_col = "logFC"
        elif "logfc" in de_df.columns:
            lfc_col = "logfc"
        else:
            lfc_col = [c for c in de_df.columns if "log" in c.lower()][0]
        
        # Gene symbol column - prefer gene_symbol, then hgnc_symbol, then gene
        if "gene_symbol" in de_df.columns:
            gene_col = "gene_symbol"
        elif "hgnc_symbol" in de_df.columns:
            gene_col = "hgnc_symbol"
        else:
            gene_col = "gene"
        
        print(f"    {tissue}: Using columns: pval={pval_col}, lfc={lfc_col}, gene={gene_col}")
        
        # Filter by p-value and split by direction
        sig_genes = de_df[de_df[pval_col] < p_threshold].copy()
        
        up_genes = sig_genes[sig_genes[lfc_col] > 0][gene_col].dropna()
        down_genes = sig_genes[sig_genes[lfc_col] < 0][gene_col].dropna()
        
        # Clean gene names - keep only valid gene symbols
        up_genes = [str(g).strip() for g in up_genes 
                   if pd.notna(g) and str(g).strip() != "" 
                   and not str(g).startswith("ENSG")]
        down_genes = [str(g).strip() for g in down_genes 
                     if pd.notna(g) and str(g).strip() != "" 
                     and not str(g).startswith("ENSG")]
        
        gene_sets[f"AMD-{tissue}-UP"] = list(set(up_genes))
        gene_sets[f"AMD-{tissue}-DOWN"] = list(set(down_genes))
        
        print(f"    {tissue}: Found {len(gene_sets[f'AMD-{tissue}-UP'])} UP, {len(gene_sets[f'AMD-{tissue}-DOWN'])} DOWN genes at p < {p_threshold}")
    
    return gene_sets

#%%
def load_axis_gene_sets() -> dict:
    """Load pre-curated axis gene sets (Senescence, Redox, Inflammation)."""
    gene_sets = {}
    
    axis_mapping = {
        "Senescence-PRO": "senescence_pro",
        "Senescence-ANTI": "senescence_anti",
        "Redox-PRO": "redox_pro",
        "Redox-ANTI": "redox_anti",
        "Inflammation-PRO": "inflammation_pro",
        "Inflammation-ANTI": "inflammation_anti",
    }
    
    for name, source_key in axis_mapping.items():
        source_file = GENE_SET_SOURCES[source_key]
        if source_file.exists():
            df = pd.read_csv(source_file)
            # Get the gene column (first column usually)
            gene_col = df.columns[0]
            genes = df[gene_col].dropna().unique().tolist()
            gene_sets[name] = genes
        else:
            print(f"WARNING: {source_file} not found")
            gene_sets[name] = []
    
    return gene_sets

#%%
def load_all_gene_sets() -> dict:
    """Load all 13 gene sets from their respective sources."""
    print("Loading gene sets...")
    
    all_sets = {}
    
    # Age gene sets
    print("  Loading Age gene sets...")
    age_sets = load_age_gene_sets()
    all_sets.update(age_sets)
    
    # ARE gene set
    print("  Loading ARE gene set...")
    all_sets["ARE"] = load_are_gene_set()
    
    # AMD gene sets
    print("  Extracting AMD gene sets (p < 0.01)...")
    amd_sets = extract_amd_gene_sets()
    all_sets.update(amd_sets)
    
    # Axis gene sets
    print("  Loading Axis gene sets...")
    axis_sets = load_axis_gene_sets()
    all_sets.update(axis_sets)
    
    return all_sets

#%%
def create_registry(gene_sets: dict) -> pd.DataFrame:
    """Create a registry DataFrame with metadata for all gene sets."""
    records = []
    
    for name, genes in gene_sets.items():
        meta = GENE_SET_CATEGORIES.get(name, {})
        records.append({
            "name": name,
            "category": meta.get("category", "Unknown"),
            "direction": meta.get("direction", "Unknown"),
            "description": meta.get("description", ""),
            "size": len(genes),
        })
    
    return pd.DataFrame(records)

#%%
def save_gene_sets(gene_sets: dict, output_dir: Path):
    """Save all gene sets as individual CSV files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for name, genes in gene_sets.items():
        df = pd.DataFrame({"gene": sorted(genes)})
        out_file = output_dir / f"{name}.csv"
        df.to_csv(out_file, index=False)
        print(f"  Saved {name}: {len(genes)} genes -> {out_file.name}")

#%%
def validate_gene_sets(gene_sets: dict) -> bool:
    """Validate that gene sets meet expected criteria."""
    print("\n" + "="*60)
    print("GENE SET VALIDATION")
    print("="*60)
    
    all_valid = True
    
    # Expected approximate sizes
    expected = {
        "Age-UP": (300, 350),
        "Age-DOWN": (140, 160),
        "ARE": (230, 260),
        "AMD-Macula-UP": (200, 1000),
        "AMD-Macula-DOWN": (200, 1000),
        "AMD-nonMacula-UP": (400, 1500),
        "AMD-nonMacula-DOWN": (400, 1500),
        "Senescence-PRO": (240, 260),
        "Senescence-ANTI": (320, 340),
        "Redox-PRO": (240, 260),
        "Redox-ANTI": (10, 30),
        "Inflammation-PRO": (560, 590),
        "Inflammation-ANTI": (100, 130),
    }
    
    for name, genes in gene_sets.items():
        size = len(genes)
        exp_range = expected.get(name, (0, 10000))
        status = "✓" if exp_range[0] <= size <= exp_range[1] else "⚠"
        
        if status == "⚠":
            all_valid = False
            print(f"  {status} {name}: {size} genes (expected {exp_range[0]}-{exp_range[1]})")
        else:
            print(f"  {status} {name}: {size} genes")
    
    print("="*60)
    
    # Check for empty gene sets
    empty_sets = [name for name, genes in gene_sets.items() if len(genes) == 0]
    if empty_sets:
        print(f"WARNING: Empty gene sets: {empty_sets}")
        all_valid = False
    
    return all_valid

#%%
def main():
    """Main function to load, validate, and save all gene sets."""
    ensure_dirs()
    
    # Load all gene sets
    gene_sets = load_all_gene_sets()
    
    # Validate
    is_valid = validate_gene_sets(gene_sets)
    
    # Create and save registry
    registry = create_registry(gene_sets)
    registry_file = OUTPUT_DIRS["gene_sets"] / "registry.csv"
    registry.to_csv(registry_file, index=False)
    print(f"\nRegistry saved to: {registry_file}")
    
    # Display registry
    print("\n" + "="*60)
    print("GENE SET REGISTRY")
    print("="*60)
    print(registry.to_string(index=False))
    
    # Save individual gene sets
    print("\nSaving individual gene set files...")
    save_gene_sets(gene_sets, OUTPUT_DIRS["gene_sets"])
    
    print(f"\n{'='*60}")
    print("Gene set preparation complete!")
    print(f"Total gene sets: {len(gene_sets)}")
    print(f"Output directory: {OUTPUT_DIRS['gene_sets']}")
    
    return gene_sets, registry

#%%
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate", action="store_true", help="Validate gene sets only")
    args = parser.parse_args()
    
    if args.validate:
        gene_sets = load_all_gene_sets()
        validate_gene_sets(gene_sets)
    else:
        main()
