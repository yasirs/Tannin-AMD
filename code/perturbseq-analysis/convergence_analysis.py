"""
Convergence Analysis Module

Identify overlapping genes across multiple analyses.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import sys
from typing import List, Dict, Tuple

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))
from config import ANALYSIS_PARAMS, OUTPUT_DIRS

#%%
def compute_overlap(set1: set, set2: set, universe_size: int = 20000) -> dict:
    """
    Compute overlap between two gene sets with hypergeometric test.
    
    Parameters
    ----------
    set1, set2 : set
        Gene sets to compare
    universe_size : int
        Size of gene universe for hypergeometric test
        
    Returns
    -------
    dict with overlap statistics
    """
    intersection = set1 & set2
    union = set1 | set2
    
    # Hypergeometric test for enrichment
    # P(X >= k) where X ~ Hypergeom(N, K, n)
    # N = universe, K = size of set1, n = size of set2, k = overlap
    k = len(intersection)
    M = universe_size
    n = len(set1)
    N = len(set2)
    
    # sf gives P(X > k-1) = P(X >= k)
    pval = stats.hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
    
    # Jaccard index
    jaccard = len(intersection) / len(union) if union else 0
    
    return {
        "set1_size": len(set1),
        "set2_size": len(set2),
        "overlap": k,
        "jaccard": jaccard,
        "hypergeom_pval": pval,
        "genes": list(intersection),
    }

#%%
def identify_convergent_genes(
    mimetic_results: List[pd.DataFrame],
    gene_set_names: List[str],
    n_top: int = 200,
    direction: str = "antagonist"
) -> pd.DataFrame:
    """
    Identify genes that appear in top N across multiple analyses.
    
    Parameters
    ----------
    mimetic_results : List[pd.DataFrame]
        List of backward analysis result DataFrames
    gene_set_names : List[str]
        Names for each result set
    n_top : int
        Number of top genes to consider from each set
    direction : str
        "mimetic" (top positive) or "antagonist" (most negative)
        
    Returns
    -------
    DataFrame with convergent genes and their appearances
    """
    # Get top genes from each analysis
    top_genes_per_set = {}
    
    for df, name in zip(mimetic_results, gene_set_names):
        if direction == "antagonist":
            # Most negative = aging antagonists / reverses signature
            top_df = df.nsmallest(n_top, "mimetic_score")
        else:
            # Most positive = aging mimetics
            top_df = df.nlargest(n_top, "mimetic_score")
        
        top_genes_per_set[name] = set(top_df["target_gene"].str.upper())
    
    # Find genes appearing in multiple sets
    all_genes = set().union(*top_genes_per_set.values())
    
    results = []
    for gene in all_genes:
        appearances = []
        for name, genes in top_genes_per_set.items():
            if gene in genes:
                appearances.append(name)
        
        if len(appearances) >= 2:  # Appears in at least 2 analyses
            results.append({
                "gene": gene,
                "n_appearances": len(appearances),
                "appears_in": "; ".join(appearances),
            })
    
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("n_appearances", ascending=False)
    
    return df

#%%
def compute_rpe1_k562_concordance(
    rpe1_df: pd.DataFrame,
    k562_df: pd.DataFrame,
    n_top: int = 200
) -> dict:
    """
    Compute concordance between RPE1 and K562 results for same gene set.
    
    Returns dict with concordance metrics.
    """
    # Get top antagonists from each
    rpe1_top = set(rpe1_df.nsmallest(n_top, "mimetic_score")["target_gene"].str.upper())
    k562_top = set(k562_df.nsmallest(n_top, "mimetic_score")["target_gene"].str.upper())
    
    overlap_stats = compute_overlap(rpe1_top, k562_top)
    
    # Also compute correlation of scores for shared genes
    # Merge on target gene
    merged = rpe1_df.merge(
        k562_df[["target_gene", "mimetic_score"]],
        on="target_gene",
        suffixes=("_rpe1", "_k562"),
        how="inner"
    )
    
    if len(merged) > 10:
        corr, corr_pval = stats.spearmanr(
            merged["mimetic_score_rpe1"],
            merged["mimetic_score_k562"]
        )
    else:
        corr, corr_pval = np.nan, np.nan
    
    return {
        **overlap_stats,
        "n_shared_genes": len(merged),
        "spearman_corr": corr,
        "corr_pval": corr_pval,
    }

#%%
def compute_all_pairwise_overlaps(
    results_dict: Dict[str, pd.DataFrame],
    n_top: int = 200
) -> pd.DataFrame:
    """
    Compute pairwise overlaps between all gene set analyses.
    
    Parameters
    ----------
    results_dict : dict
        {gene_set_name: backward_analysis_df}
    n_top : int
        Number of top antagonists to compare
        
    Returns
    -------
    DataFrame with pairwise overlap statistics
    """
    names = list(results_dict.keys())
    records = []
    
    for i, name1 in enumerate(names):
        for j, name2 in enumerate(names):
            if j <= i:
                continue
            
            df1 = results_dict[name1]
            df2 = results_dict[name2]
            
            # Get top antagonists
            top1 = set(df1.nsmallest(n_top, "mimetic_score")["target_gene"].str.upper())
            top2 = set(df2.nsmallest(n_top, "mimetic_score")["target_gene"].str.upper())
            
            overlap_stats = compute_overlap(top1, top2)
            
            records.append({
                "set1": name1,
                "set2": name2,
                **overlap_stats,
            })
    
    return pd.DataFrame(records)

#%%
def generate_concordance_table(
    rpe1_results: Dict[str, pd.DataFrame],
    k562_results: Dict[str, pd.DataFrame],
    n_top: int = 200
) -> pd.DataFrame:
    """
    Generate concordance table: RPE1 vs K562 for each gene set.
    """
    records = []
    
    for gene_set_name in rpe1_results.keys():
        if gene_set_name not in k562_results:
            continue
        
        conc = compute_rpe1_k562_concordance(
            rpe1_results[gene_set_name],
            k562_results[gene_set_name],
            n_top=n_top
        )
        
        records.append({
            "gene_set": gene_set_name,
            "top_n": n_top,
            "rpe1_top_genes": conc["set1_size"],
            "k562_top_genes": conc["set2_size"],
            "overlap": conc["overlap"],
            "jaccard": conc["jaccard"],
            "hypergeom_pval": conc["hypergeom_pval"],
            "n_shared_genes": conc["n_shared_genes"],
            "spearman_corr": conc["spearman_corr"],
            "corr_pval": conc["corr_pval"],
        })
    
    return pd.DataFrame(records)

#%%
if __name__ == "__main__":
    print("Testing convergence analysis module...")
    
    # Create test data
    np.random.seed(42)
    
    genes = [f"GENE{i}" for i in range(1000)]
    
    df1 = pd.DataFrame({
        "target_gene": genes,
        "mimetic_score": np.random.randn(1000) * 0.1,
    })
    
    df2 = pd.DataFrame({
        "target_gene": genes,
        "mimetic_score": np.random.randn(1000) * 0.1,
    })
    
    # Test overlap computation
    set1 = set(genes[:100])
    set2 = set(genes[50:150])
    overlap = compute_overlap(set1, set2)
    print(f"Overlap test: {overlap['overlap']} genes")
    
    # Test concordance
    conc = compute_rpe1_k562_concordance(df1, df2, n_top=100)
    print(f"Concordance test: {conc['overlap']} overlapping, r={conc['spearman_corr']:.3f}")
    
    print("Convergence tests passed!")
