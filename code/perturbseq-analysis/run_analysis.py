"""
Main Orchestration Script

Run complete Perturb-seq analysis on all 13 gene sets across 2 datasets.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import time
import argparse

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    OUTPUT_DIRS, GENE_SET_CATEGORIES, ANALYSIS_PARAMS,
    ensure_dirs, get_dataset_output_dir
)
from gene_set_registry import load_all_gene_sets, create_registry, validate_gene_sets
from perturbseq_loader import load_perturbseq_data, compute_gene_coverage, validate_data
from enrichment_analysis import (
    run_forward_analysis, run_backward_analysis,
    validate_backward_results, identify_top_mimetics,
    compute_summary_statistics
)
from convergence_analysis import (
    compute_rpe1_k562_concordance, generate_concordance_table,
    compute_all_pairwise_overlaps, identify_convergent_genes
)
from visualization import (
    plot_mimetic_distribution, plot_convergence_scatter,
    plot_enrichment_distribution, plot_top_genes_bar, save_figure
)

#%%
def get_paired_gene_sets(gene_sets: dict) -> dict:
    """
    Organize gene sets into UP/DOWN pairs for mimetic analysis.
    
    Returns dict: {pair_name: {"up": up_genes, "down": down_genes}}
    """
    pairs = {}
    
    # Age: Age-UP vs Age-DOWN
    if "Age-UP" in gene_sets and "Age-DOWN" in gene_sets:
        pairs["Age"] = {
            "up": gene_sets["Age-UP"],
            "down": gene_sets["Age-DOWN"],
            "category": "Age",
        }
    
    # AMD Macula
    if "AMD-Macula-UP" in gene_sets and "AMD-Macula-DOWN" in gene_sets:
        pairs["AMD-Macula"] = {
            "up": gene_sets["AMD-Macula-UP"],
            "down": gene_sets["AMD-Macula-DOWN"],
            "category": "AMD",
        }
    
    # AMD nonMacula
    if "AMD-nonMacula-UP" in gene_sets and "AMD-nonMacula-DOWN" in gene_sets:
        pairs["AMD-nonMacula"] = {
            "up": gene_sets["AMD-nonMacula-UP"],
            "down": gene_sets["AMD-nonMacula-DOWN"],
            "category": "AMD",
        }
    
    # Senescence
    if "Senescence-PRO" in gene_sets and "Senescence-ANTI" in gene_sets:
        pairs["Senescence"] = {
            "up": gene_sets["Senescence-PRO"],
            "down": gene_sets["Senescence-ANTI"],
            "category": "Senescence",
        }
    
    # Redox
    if "Redox-PRO" in gene_sets and "Redox-ANTI" in gene_sets:
        pairs["Redox"] = {
            "up": gene_sets["Redox-PRO"],
            "down": gene_sets["Redox-ANTI"],
            "category": "Redox",
        }
    
    # Inflammation
    if "Inflammation-PRO" in gene_sets and "Inflammation-ANTI" in gene_sets:
        pairs["Inflammation"] = {
            "up": gene_sets["Inflammation-PRO"],
            "down": gene_sets["Inflammation-ANTI"],
            "category": "Inflammation",
        }
    
    return pairs

#%%
def run_single_analysis(
    perturbseq_data: dict,
    dataset_name: str,
    gene_set_pair: dict,
    pair_name: str,
    are_genes: list = None
) -> dict:
    """
    Run complete analysis for a single gene set pair on one dataset.
    
    Returns dict with all results.
    """
    print(f"\n{'='*60}")
    print(f"Analyzing: {pair_name} on {dataset_name}")
    print("="*60)
    
    up_genes = gene_set_pair["up"]
    down_genes = gene_set_pair["down"]
    category = gene_set_pair["category"]
    
    # Get output directory
    out_dir = get_dataset_output_dir(pair_name, dataset_name, "backward")
    
    # --- Backward Analysis ---
    print("\n[1] Running backward analysis (all knockdowns)...")
    start_time = time.time()
    
    backward_df = run_backward_analysis(
        perturbseq_data,
        up_genes,
        down_genes,
        gene_set_name=pair_name
    )
    
    elapsed = time.time() - start_time
    print(f"    Completed in {elapsed:.1f} seconds")
    
    # Validate results
    validate_backward_results(backward_df, pair_name)
    
    # Save results
    backward_file = out_dir / f"backward_{pair_name}_mimetics.csv"
    backward_df.to_csv(backward_file, index=False)
    print(f"    Saved: {backward_file}")
    
    # Get top mimetics and antagonists
    top_n = ANALYSIS_PARAMS["top_n_mimetics"]
    top_mimetics, top_antagonists = identify_top_mimetics(backward_df, n_top=top_n)
    
    # Save top genes
    top_mimetics.to_csv(out_dir / f"top_{top_n}_mimetics.csv", index=False)
    top_antagonists.to_csv(out_dir / f"top_{top_n}_antagonists.csv", index=False)
    
    # --- Generate Figures ---
    print("\n[2] Generating figures...")
    
    # Mimetic score distribution
    fig, ax = plot_mimetic_distribution(
        backward_df,
        title=f"{pair_name} Mimetic Score Distribution ({dataset_name})",
        output_path=out_dir / f"Fig_{pair_name}_distribution"
    )
    
    # Top genes bar plot
    fig, ax = plot_top_genes_bar(
        backward_df,
        n_top=15,
        title=f"{pair_name} Top Mimetics/Antagonists ({dataset_name})",
        output_path=out_dir / f"Fig_{pair_name}_top_genes"
    )
    
    # --- Forward Analysis (if ARE genes provided) ---
    forward_df = None
    if are_genes and len(up_genes) > 0:
        print("\n[3] Running forward analysis...")
        
        forward_out_dir = get_dataset_output_dir(pair_name, dataset_name, "forward")
        
        # Check if query genes have knockdown data
        from perturbseq_loader import get_perturbation_indices
        pert_indices, matched = get_perturbation_indices(
            perturbseq_data["target_genes"],
            up_genes + down_genes
        )
        
        if len(pert_indices) > 5:
            forward_df = run_forward_analysis(
                perturbseq_data,
                query_gene_set=up_genes + down_genes,
                target_gene_set=are_genes,
                query_name=pair_name,
                target_name="ARE"
            )
            
            forward_file = forward_out_dir / f"forward_{pair_name}_ARE_enrichment.csv"
            forward_df.to_csv(forward_file, index=False)
            print(f"    Saved: {forward_file}")
        else:
            print(f"    Skipping forward analysis - insufficient knockdown data")
    
    # Compute summary statistics
    stats = compute_summary_statistics(backward_df, "mimetic_score")
    
    return {
        "pair_name": pair_name,
        "dataset": dataset_name,
        "backward_df": backward_df,
        "forward_df": forward_df,
        "top_mimetics": top_mimetics,
        "top_antagonists": top_antagonists,
        "stats": stats,
        "output_dir": out_dir,
    }

#%%
def run_all_analyses(use_cache: bool = True) -> dict:
    """
    Run complete analysis on all gene sets across both datasets.
    """
    print("="*70)
    print("COMPREHENSIVE PERTURB-SEQ GENE SET ANALYSIS")
    print("="*70)
    
    total_start = time.time()
    
    # Ensure output directories exist
    ensure_dirs()
    
    # --- Load Gene Sets ---
    print("\n[PHASE 1] Loading gene sets...")
    gene_sets = load_all_gene_sets()
    validate_gene_sets(gene_sets)
    
    # Get ARE genes for forward analysis
    are_genes = gene_sets.get("ARE", [])
    
    # Organize into pairs
    pairs = get_paired_gene_sets(gene_sets)
    print(f"\nPrepared {len(pairs)} gene set pairs for analysis")
    
    # --- Load Perturb-seq Data ---
    print("\n[PHASE 2] Loading Perturb-seq data...")
    
    rpe1_data = load_perturbseq_data("RPE1", use_cache=use_cache)
    k562_data = load_perturbseq_data("K562_GWPS", use_cache=use_cache)
    
    # Compute coverage for all gene sets
    print("\n[PHASE 2.1] Computing gene set coverage...")
    
    coverage_records = []
    for name, genes in gene_sets.items():
        for dataset_name, data in [("RPE1", rpe1_data), ("K562_GWPS", k562_data)]:
            cov = compute_gene_coverage(data, genes, name)
            cov["dataset"] = dataset_name
            coverage_records.append(cov)
    
    coverage_df = pd.DataFrame(coverage_records)
    coverage_file = OUTPUT_DIRS["gene_sets"] / "coverage_in_perturbseq.csv"
    coverage_df.to_csv(coverage_file, index=False)
    print(f"Saved coverage summary: {coverage_file}")
    
    # --- Run Analyses ---
    print(f"\n[PHASE 3] Running analyses ({len(pairs)} pairs × 2 datasets = {len(pairs)*2} analyses)...")
    
    all_results = {
        "RPE1": {},
        "K562_GWPS": {},
    }
    
    for pair_name, pair_info in pairs.items():
        # Run on RPE1
        result_rpe1 = run_single_analysis(
            rpe1_data, "RPE1", pair_info, pair_name, are_genes
        )
        all_results["RPE1"][pair_name] = result_rpe1
        
        # Run on K562
        result_k562 = run_single_analysis(
            k562_data, "K562_GWPS", pair_info, pair_name, are_genes
        )
        all_results["K562_GWPS"][pair_name] = result_k562
    
    # --- Convergence Analysis ---
    print(f"\n[PHASE 4] Running convergence analysis...")
    
    conc_dir = OUTPUT_DIRS["concordance"]
    conc_dir.mkdir(parents=True, exist_ok=True)
    
    # RPE1 vs K562 concordance for each gene set
    rpe1_backward = {name: r["backward_df"] for name, r in all_results["RPE1"].items()}
    k562_backward = {name: r["backward_df"] for name, r in all_results["K562_GWPS"].items()}
    
    concordance_df = generate_concordance_table(rpe1_backward, k562_backward)
    concordance_df.to_csv(conc_dir / "rpe1_vs_k562_concordance.csv", index=False)
    print(f"Saved: rpe1_vs_k562_concordance.csv")
    
    # Cross-gene-set overlaps (within each dataset)
    for dataset_name, results in all_results.items():
        backward_dfs = {name: r["backward_df"] for name, r in results.items()}
        overlap_df = compute_all_pairwise_overlaps(backward_dfs)
        overlap_df.to_csv(conc_dir / f"cross_geneset_overlaps_{dataset_name}.csv", index=False)
        print(f"Saved: cross_geneset_overlaps_{dataset_name}.csv")
    
    # Identify convergent genes across all analyses
    all_backward_dfs = []
    all_names = []
    
    for dataset_name, results in all_results.items():
        for pair_name, r in results.items():
            all_backward_dfs.append(r["backward_df"])
            all_names.append(f"{pair_name}_{dataset_name}")
    
    convergent_genes = identify_convergent_genes(all_backward_dfs, all_names)
    convergent_genes.to_csv(conc_dir / "convergent_genes_all.csv", index=False)
    print(f"Saved: convergent_genes_all.csv")
    print(f"Found {len(convergent_genes)} genes appearing in multiple analyses")
    
    # --- Summary ---
    total_elapsed = time.time() - total_start
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"Total time: {total_elapsed/60:.1f} minutes")
    print(f"Gene set pairs analyzed: {len(pairs)}")
    print(f"Datasets: RPE1, K562_GWPS")
    print(f"Output directory: {OUTPUT_DIRS['Age'].parent}")
    
    return all_results

#%%
def main():
    parser = argparse.ArgumentParser(description="Run Perturb-seq gene set analysis")
    parser.add_argument("--no-cache", action="store_true", help="Don't use cached data")
    parser.add_argument("--pair", type=str, help="Run only specific pair (e.g., 'Age', 'Senescence')")
    args = parser.parse_args()
    
    run_all_analyses(use_cache=not args.no_cache)

#%%
if __name__ == "__main__":
    main()
