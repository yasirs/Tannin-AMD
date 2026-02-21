"""
Enrichment Analysis Module

Core GSEA-based enrichment scoring for Perturb-seq analysis.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import sys
from typing import List, Tuple, Dict
import warnings

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))
from config import ANALYSIS_PARAMS, VIZ_SETTINGS

#%%
def compute_gsea_enrichment(
    ranked_genes: np.ndarray,
    gene_indices: List[int],
    weighted: bool = True
) -> float:
    """
    Compute GSEA-style enrichment score for a gene set.
    
    Parameters
    ----------
    ranked_genes : np.ndarray
        Indices of genes sorted by their rank (from highest to lowest expression)
    gene_indices : List[int]
        Indices of genes in the query set
    weighted : bool
        Whether to use weighted scoring (uses rank position)
        
    Returns
    -------
    float : Enrichment score between -1 and 1
    """
    n_genes = len(ranked_genes)
    n_hit = len(gene_indices)
    
    if n_hit == 0 or n_genes == 0:
        return 0.0
    
    # Create hit indicator
    hit_indicator = np.isin(ranked_genes, gene_indices)
    
    # Running sum calculation
    n_miss = n_genes - n_hit
    
    if n_miss == 0:
        return 0.0
    
    # Step sizes
    hit_step = 1.0 / n_hit
    miss_step = 1.0 / n_miss
    
    # Calculate running sum
    running_sum = np.zeros(n_genes)
    current_sum = 0.0
    
    for i in range(n_genes):
        if hit_indicator[i]:
            current_sum += hit_step
        else:
            current_sum -= miss_step
        running_sum[i] = current_sum
    
    # Return the maximum deviation from zero
    max_pos = np.max(running_sum)
    max_neg = np.min(running_sum)
    
    # Return whichever has larger absolute value
    if abs(max_pos) > abs(max_neg):
        return max_pos
    else:
        return max_neg

#%%
def rank_genes_by_expression(expression_vector: np.ndarray) -> np.ndarray:
    """
    Rank genes by expression level (highest to lowest).
    
    Returns indices that would sort the genes from highest to lowest expression.
    """
    return np.argsort(-expression_vector)

#%%
def compute_enrichment_for_perturbation(
    expression_vector: np.ndarray,
    gene_set_indices: List[int]
) -> float:
    """
    Compute enrichment score for a single perturbation against a gene set.
    """
    ranked_genes = rank_genes_by_expression(expression_vector)
    return compute_gsea_enrichment(ranked_genes, gene_set_indices)

#%%
def compute_enrichment_with_pvalue(
    expression_vector: np.ndarray,
    gene_set_indices: List[int],
    n_permutations: int = 1000
) -> Tuple[float, float]:
    """
    Compute enrichment score with permutation-based p-value.
    
    Parameters
    ----------
    expression_vector : np.ndarray
        Gene expression for this perturbation
    gene_set_indices : List[int]
        Indices of genes in the query set
    n_permutations : int
        Number of permutations for p-value calculation
        
    Returns
    -------
    tuple: (enrichment_score, p_value)
    """
    n_genes = len(expression_vector)
    n_set = len(gene_set_indices)
    
    if n_set == 0:
        return 0.0, 1.0
    
    # Real enrichment score
    ranked_genes = rank_genes_by_expression(expression_vector)
    real_es = compute_gsea_enrichment(ranked_genes, gene_set_indices)
    
    # Null distribution by permuting gene set membership
    null_scores = np.zeros(n_permutations)
    
    for i in range(n_permutations):
        # Random gene set of same size
        random_indices = np.random.choice(n_genes, n_set, replace=False).tolist()
        null_scores[i] = compute_gsea_enrichment(ranked_genes, random_indices)
    
    # Two-tailed p-value
    if real_es >= 0:
        p_value = np.mean(null_scores >= real_es)
    else:
        p_value = np.mean(null_scores <= real_es)
    
    # Ensure minimum p-value based on number of permutations
    p_value = max(p_value, 1.0 / (n_permutations + 1))
    
    return real_es, p_value

#%%
def compute_mimetic_score_with_pvalue(
    expression_vector: np.ndarray,
    up_gene_indices: List[int],
    down_gene_indices: List[int],
    n_permutations: int = 100
) -> Tuple[float, float, float, float, float, float]:
    """
    Compute mimetic score with p-values for UP and DOWN enrichments.
    
    Returns
    -------
    tuple: (mimetic_score, up_es, down_es, up_pval, down_pval, mimetic_pval)
    """
    # Compute enrichments with p-values
    up_es, up_pval = compute_enrichment_with_pvalue(
        expression_vector, up_gene_indices, n_permutations
    ) if up_gene_indices else (0.0, 1.0)
    
    down_es, down_pval = compute_enrichment_with_pvalue(
        expression_vector, down_gene_indices, n_permutations
    ) if down_gene_indices else (0.0, 1.0)
    
    mimetic_score = up_es - down_es
    
    # Combined p-value using Fisher's method
    if up_pval > 0 and down_pval > 0:
        chi2_stat = -2 * (np.log(up_pval) + np.log(down_pval))
        mimetic_pval = 1 - stats.chi2.cdf(chi2_stat, df=4)
    else:
        mimetic_pval = min(up_pval, down_pval)
    
    return mimetic_score, up_es, down_es, up_pval, down_pval, mimetic_pval


#%%
def compute_mimetic_score(
    expression_vector: np.ndarray,
    up_gene_indices: List[int],
    down_gene_indices: List[int]
) -> Tuple[float, float, float]:
    """
    Compute mimetic score for a perturbation.
    
    Mimetic score = Enrichment(UP genes) - Enrichment(DOWN genes)
    
    Positive: perturbation upregulates UP genes and downregulates DOWN genes
    Negative: perturbation has opposite effect (antagonist)
    
    Returns
    -------
    tuple: (mimetic_score, up_enrichment, down_enrichment)
    """
    ranked_genes = rank_genes_by_expression(expression_vector)
    
    up_enrichment = compute_gsea_enrichment(ranked_genes, up_gene_indices) if up_gene_indices else 0.0
    down_enrichment = compute_gsea_enrichment(ranked_genes, down_gene_indices) if down_gene_indices else 0.0
    
    mimetic_score = up_enrichment - down_enrichment
    
    return mimetic_score, up_enrichment, down_enrichment

#%%
def run_forward_analysis(
    perturbseq_data: dict,
    query_gene_set: List[str],
    target_gene_set: List[str],
    query_name: str = "Query",
    target_name: str = "Target"
) -> pd.DataFrame:
    """
    Forward Analysis: For genes in query set that have knockdown data,
    compute their effect on the target gene set.
    
    Parameters
    ----------
    perturbseq_data : dict
        Loaded Perturb-seq data from perturbseq_loader
    query_gene_set : List[str]
        Genes to knock down (query set)
    target_gene_set : List[str]
        Genes to measure enrichment of (target set)
    query_name : str
        Name of query set for labeling
    target_name : str
        Name of target set for labeling
        
    Returns
    -------
    pd.DataFrame with columns:
        - perturbation_id
        - target_gene (knocked down gene)
        - enrichment_score
        - in_query_set (which query genes have knockdown data)
    """
    from perturbseq_loader import get_gene_indices, get_perturbation_indices
    
    expression = perturbseq_data["expression"]
    gene_names = perturbseq_data["gene_names"]
    target_genes = perturbseq_data["target_genes"]
    perturbation_ids = perturbseq_data["perturbation_ids"]
    
    # Get indices of target gene set in expression matrix
    target_indices, matched_targets, _ = get_gene_indices(gene_names, target_gene_set)
    
    if len(target_indices) < ANALYSIS_PARAMS["min_gene_set_size"]:
        warnings.warn(f"Target gene set too small after mapping: {len(target_indices)} genes")
    
    # Get perturbations that knock down query genes
    pert_indices, matched_query = get_perturbation_indices(target_genes, query_gene_set)
    
    print(f"  Forward analysis: {len(matched_query)} query genes have knockdown data")
    print(f"  Target set: {len(target_indices)} genes mapped")
    
    # Compute enrichment for each perturbation
    results = []
    for idx in pert_indices:
        expr = expression[idx, :]
        es = compute_enrichment_for_perturbation(expr, target_indices)
        
        results.append({
            "perturbation_id": perturbation_ids[idx],
            "target_gene": target_genes[idx],
            "enrichment_score": es,
            "in_query_set": True,
        })
    
    df = pd.DataFrame(results)
    return df

#%%
def run_backward_analysis(
    perturbseq_data: dict,
    up_genes: List[str],
    down_genes: List[str],
    gene_set_name: str = "Query"
) -> pd.DataFrame:
    """
    Backward Analysis: For ALL knockdowns, compute mimetic scores.
    
    Mimetic score = how much the knockdown phenocopies the signature
    (upregulates UP genes and downregulates DOWN genes)
    
    Parameters
    ----------
    perturbseq_data : dict
        Loaded Perturb-seq data
    up_genes : List[str]
        Genes that are UP in the signature
    down_genes : List[str]
        Genes that are DOWN in the signature
    gene_set_name : str
        Name for labeling
        
    Returns
    -------
    pd.DataFrame with mimetic scores for all knockdowns
    """
    from perturbseq_loader import get_gene_indices
    
    expression = perturbseq_data["expression"]
    gene_names = perturbseq_data["gene_names"]
    target_genes = perturbseq_data["target_genes"]
    perturbation_ids = perturbseq_data["perturbation_ids"]
    n_perts = len(perturbation_ids)
    
    # Get gene set indices
    up_indices, matched_up, _ = get_gene_indices(gene_names, up_genes)
    down_indices, matched_down, _ = get_gene_indices(gene_names, down_genes)
    
    print(f"  Backward analysis: {gene_set_name}")
    print(f"    UP genes: {len(matched_up)}/{len(up_genes)} mapped")
    print(f"    DOWN genes: {len(matched_down)}/{len(down_genes)} mapped")
    
    # Compute mimetic scores for all perturbations
    results = []
    n_skipped = 0
    
    for i in range(n_perts):
        if i % 1000 == 0:
            print(f"    Processing perturbation {i}/{n_perts}...", end="\r")
        
        # Skip non-targeting controls
        target = str(target_genes[i]).lower()
        if "non-targeting" in target or "nontargeting" in target or target == "control":
            n_skipped += 1
            continue
        
        expr = expression[i, :]
        mimetic, up_es, down_es = compute_mimetic_score(expr, up_indices, down_indices)
        
        results.append({
            "perturbation_id": perturbation_ids[i],
            "target_gene": target_genes[i],
            "mimetic_score": mimetic,
            "up_enrichment": up_es,
            "down_enrichment": down_es,
        })
    
    print(f"    Processed {n_perts - n_skipped} perturbations (skipped {n_skipped} controls).")
    
    df = pd.DataFrame(results)
    
    # Sort by mimetic score
    df = df.sort_values("mimetic_score", ascending=False)
    
    return df

#%%
def identify_top_mimetics(df: pd.DataFrame, n_top: int = 200) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Identify top mimetics and antagonists from backward analysis results.
    
    Returns
    -------
    tuple: (top_mimetics_df, top_antagonists_df)
    """
    df_sorted = df.sort_values("mimetic_score", ascending=False)
    
    top_mimetics = df_sorted.head(n_top).copy()
    top_antagonists = df_sorted.tail(n_top).copy()
    
    return top_mimetics, top_antagonists

#%%
def compute_summary_statistics(df: pd.DataFrame, score_col: str = "mimetic_score") -> dict:
    """Compute summary statistics for quality checks."""
    scores = df[score_col].values
    
    return {
        "n": len(scores),
        "mean": float(np.mean(scores)),
        "std": float(np.std(scores)),
        "min": float(np.min(scores)),
        "max": float(np.max(scores)),
        "median": float(np.median(scores)),
        "skewness": float(stats.skew(scores)),
        "kurtosis": float(stats.kurtosis(scores)),
    }

#%%
def validate_backward_results(df: pd.DataFrame, gene_set_name: str = "") -> bool:
    """
    Validate backward analysis results with sanity checks.
    """
    stats_dict = compute_summary_statistics(df, "mimetic_score")
    
    print(f"\n  Sanity checks for {gene_set_name}:")
    
    all_valid = True
    
    # Mean should be near zero
    if abs(stats_dict["mean"]) > 0.1:
        print(f"    ⚠ Mean mimetic score not near zero: {stats_dict['mean']:.4f}")
        all_valid = False
    else:
        print(f"    ✓ Mean: {stats_dict['mean']:.4f}")
    
    # Std should be reasonable (0.05 - 0.3)
    if not (0.02 < stats_dict["std"] < 0.5):
        print(f"    ⚠ Unusual std: {stats_dict['std']:.4f}")
    else:
        print(f"    ✓ Std: {stats_dict['std']:.4f}")
    
    # Distribution should be roughly symmetric
    if abs(stats_dict["skewness"]) > 1.0:
        print(f"    ⚠ High skewness: {stats_dict['skewness']:.4f}")
    else:
        print(f"    ✓ Skewness: {stats_dict['skewness']:.4f}")
    
    print(f"    Range: [{stats_dict['min']:.4f}, {stats_dict['max']:.4f}]")
    
    return all_valid

#%%
if __name__ == "__main__":
    # Quick test
    print("Testing enrichment functions...")
    
    # Create test data
    n_genes = 1000
    ranked = np.arange(n_genes)
    
    # Gene set at top = positive enrichment
    top_genes = list(range(50))
    es_top = compute_gsea_enrichment(ranked, top_genes)
    print(f"Top genes enrichment: {es_top:.3f} (expected positive)")
    
    # Gene set at bottom = negative enrichment
    bottom_genes = list(range(950, 1000))
    es_bottom = compute_gsea_enrichment(ranked, bottom_genes)
    print(f"Bottom genes enrichment: {es_bottom:.3f} (expected negative)")
    
    # Random genes = near zero
    np.random.seed(42)
    random_genes = np.random.choice(n_genes, 50, replace=False).tolist()
    es_random = compute_gsea_enrichment(ranked, random_genes)
    print(f"Random genes enrichment: {es_random:.3f} (expected ~0)")
    
    print("\nTests passed!")
