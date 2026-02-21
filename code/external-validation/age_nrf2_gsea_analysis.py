#!/usr/bin/env python3
"""
Age-NRF2-ARE Analysis using GSEA-based methods
Implements Forward, Backward (Age), Backward (PRG4), and Convergence analyses
"""

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
import gseapy as gp
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = Path("/home/ysuhail/work/Tannin-AMD")
RESULTS_DIR = BASE_DIR / "results/cohort-GSE29801/age_nrf2_analysis"
PERTURBSEQ_FILE = BASE_DIR / "data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad"

print("=" * 60)
print("AGE-NRF2-ARE GSEA-BASED ANALYSIS")
print("=" * 60)

#%% ===========================================================================
# STEP 1: LOAD DATA
# ===========================================================================

print("\n[1] LOADING DATA...")

# Load Perturb-seq
print("  Loading Perturb-seq data...")
adata = ad.read_h5ad(PERTURBSEQ_FILE)
print(f"    Perturbations: {adata.n_obs}")
print(f"    Genes: {adata.n_vars}")

# Extract gene names from var
gene_names = adata.var['gene_name'].values
perturbation_ids = adata.obs.index.values

# Create gene name to index mapping
gene_to_idx = {g: i for i, g in enumerate(gene_names)}

# Extract perturbation target genes from obs index
# Format: "0_A1BG_P1_ENSG00000121410" -> "A1BG"
perturbation_targets = []
for pid in perturbation_ids:
    parts = pid.split('_')
    if len(parts) >= 2:
        perturbation_targets.append(parts[1])
    else:
        perturbation_targets.append(pid)

perturbation_df = pd.DataFrame({
    'perturbation_id': perturbation_ids,
    'target_gene': perturbation_targets
})
print(f"    Unique target genes: {perturbation_df['target_gene'].nunique()}")

# Load age-DE genes
print("  Loading age-DE genes...")
age_de = pd.read_csv(RESULTS_DIR / "age_de_extramacular.csv")
age_UP = set(age_de[age_de['direction'] == 'UP']['gene'].dropna().unique())
age_DOWN = set(age_de[age_de['direction'] == 'DOWN']['gene'].dropna().unique())
print(f"    Age-UP genes: {len(age_UP)}")
print(f"    Age-DOWN genes: {len(age_DOWN)}")

# Load ARE genes
print("  Loading ARE gene set...")
are_genes = set(pd.read_csv(RESULTS_DIR / "ARE_gene_set.csv")['gene'].dropna().unique())
print(f"    ARE genes: {len(are_genes)}")

# Load PRG4 signature
print("  Loading PRG4 signature...")
prg4_de = pd.read_csv(BASE_DIR / "data/RPE_DE_results.csv")
prg4_de = prg4_de[prg4_de['hgnc_symbol'].notna() & (prg4_de['hgnc_symbol'] != '')]
# Use rescue comparison (H2O2+PRG4 vs H2O2)
prg4_sig = prg4_de[['hgnc_symbol', 'H2O2PRG4_vs_H2O2.log2FoldChange', 'H2O2PRG4_vs_H2O2.pvalue']].copy()
prg4_sig.columns = ['gene', 'lfc', 'pval']
prg4_sig = prg4_sig.dropna(subset=['lfc'])
# Define PRG4-UP and PRG4-DOWN based on direction (using lenient threshold for exploratory)
prg4_UP = set(prg4_sig[prg4_sig['lfc'] > 0.5]['gene'].unique())
prg4_DOWN = set(prg4_sig[prg4_sig['lfc'] < -0.5]['gene'].unique())
print(f"    PRG4-UP genes (lfc>0.5): {len(prg4_UP)}")
print(f"    PRG4-DOWN genes (lfc<-0.5): {len(prg4_DOWN)}")

# Check overlap with perturbseq genes
ps_genes = set(gene_names)
print("\n  Gene overlaps with Perturb-seq:")
print(f"    Age-UP in PS: {len(age_UP & ps_genes)}/{len(age_UP)}")
print(f"    Age-DOWN in PS: {len(age_DOWN & ps_genes)}/{len(age_DOWN)}")
print(f"    ARE in PS: {len(are_genes & ps_genes)}/{len(are_genes)}")

#%% ===========================================================================
# STEP 2: HELPER FUNCTIONS
# ===========================================================================

def get_knockdown_profile(adata, perturbation_idx, gene_names):
    """Get the expression profile for a knockdown as a ranked gene list."""
    # Expression values for this perturbation
    expr = adata.X[perturbation_idx, :].toarray().flatten() if hasattr(adata.X, 'toarray') else adata.X[perturbation_idx, :]
    
    # Create ranked gene list (higher expression = higher rank)
    # For GSEA, we typically want a score, not raw expression
    # We'll use z-scores or the expression values directly
    ranked = pd.Series(expr, index=gene_names).sort_values(ascending=False)
    return ranked

def run_gsea_prerank(ranked_genes, gene_set, gene_set_name="query"):
    """Run GSEA prerank on a ranked gene list with a gene set."""
    # Filter gene set to only include genes in the ranked list
    gene_set_filtered = [g for g in gene_set if g in ranked_genes.index]
    
    if len(gene_set_filtered) < 5:
        return {'NES': np.nan, 'pval': np.nan, 'fdr': np.nan, 'size': len(gene_set_filtered), 'leading_edge': []}
    
    try:
        # Use gseapy's prerank
        res = gp.prerank(
            rnk=ranked_genes,
            gene_sets={gene_set_name: gene_set_filtered},
            min_size=5,
            max_size=500,
            permutation_num=100,  # Fast for exploratory
            outdir=None,
            verbose=False
        )
        
        if gene_set_name in res.res2d.index:
            row = res.res2d.loc[gene_set_name]
            return {
                'NES': row['NES'],
                'pval': row['NOM p-val'],
                'fdr': row['FDR q-val'],
                'size': int(row['Tag %'].split('/')[1]) if isinstance(row['Tag %'], str) else len(gene_set_filtered),
                'leading_edge': row['Lead_genes'].split(';') if pd.notna(row['Lead_genes']) else []
            }
    except Exception as e:
        pass
    
    return {'NES': np.nan, 'pval': np.nan, 'fdr': np.nan, 'size': len(gene_set_filtered), 'leading_edge': []}

def quick_enrichment_score(ranked_genes, gene_set):
    """Quick enrichment approximation using rank-based method."""
    gene_set_in_list = [g for g in gene_set if g in ranked_genes.index]
    if len(gene_set_in_list) < 3:
        return np.nan, len(gene_set_in_list)
    
    # Get ranks of gene set members (lower rank = higher in list)
    n_genes = len(ranked_genes)
    gene_ranks = ranked_genes.rank(ascending=False)  # 1 = highest
    
    set_ranks = gene_ranks[gene_set_in_list].values
    
    # Compare mean rank of set to expected (middle of list)
    expected_mean = (n_genes + 1) / 2
    observed_mean = np.mean(set_ranks)
    
    # Normalize to -1 to +1 scale (positive = enriched at top)
    enrichment = (expected_mean - observed_mean) / expected_mean
    
    return enrichment, len(gene_set_in_list)

#%% ===========================================================================
# ANALYSIS 1: FORWARD APPROACH
# Age-DE gene knockdowns -> ARE enrichment
# ===========================================================================

print("\n" + "=" * 60)
print("[ANALYSIS 1] FORWARD: Age-DE knockdowns -> ARE effect")
print("=" * 60)

# Find age-DE genes that have knockdowns
age_de_genes = age_UP | age_DOWN
age_de_with_kd = []

for idx, row in perturbation_df.iterrows():
    if row['target_gene'] in age_de_genes:
        direction = 'UP' if row['target_gene'] in age_UP else 'DOWN'
        age_de_with_kd.append({
            'idx': idx,
            'perturbation_id': row['perturbation_id'],
            'target_gene': row['target_gene'],
            'age_direction': direction
        })

age_de_kd_df = pd.DataFrame(age_de_with_kd)
print(f"Age-DE genes with knockdowns: {len(age_de_kd_df)} perturbations ({age_de_kd_df['target_gene'].nunique()} unique genes)")
print(f"  Age-UP: {(age_de_kd_df['age_direction'] == 'UP').sum()}")
print(f"  Age-DOWN: {(age_de_kd_df['age_direction'] == 'DOWN').sum()}")

# Run GSEA for each knockdown (using quick approximation for speed)
print("\nRunning enrichment analysis for each knockdown...")
forward_results = []

for i, row in age_de_kd_df.iterrows():
    # Get knockdown profile
    profile = get_knockdown_profile(adata, row['idx'], gene_names)
    
    # Quick enrichment score
    es, n_are = quick_enrichment_score(profile, are_genes)
    
    forward_results.append({
        'target_gene': row['target_gene'],
        'perturbation_id': row['perturbation_id'],
        'age_direction': row['age_direction'],
        'ARE_enrichment': es,
        'n_ARE_genes': n_are
    })

forward_df = pd.DataFrame(forward_results)
forward_df = forward_df.sort_values('ARE_enrichment', ascending=False)

print("\n--- SANITY CHECK: Known biology ---")
# Check KEAP1 (should activate ARE when knocked down)
keap1_kd = forward_df[forward_df['target_gene'] == 'KEAP1']
if len(keap1_kd) > 0:
    print(f"  KEAP1 knockdown ARE enrichment: {keap1_kd['ARE_enrichment'].values[0]:.3f}")
    print(f"    Expected: POSITIVE (removing KEAP1 activates NRF2/ARE)")
else:
    print("  KEAP1 not in age-DE (expected - checking general knockdowns later)")

# Check NFE2L2 
nfe2l2_kd = forward_df[forward_df['target_gene'] == 'NFE2L2']
if len(nfe2l2_kd) > 0:
    print(f"  NFE2L2 knockdown ARE enrichment: {nfe2l2_kd['ARE_enrichment'].values[0]:.3f}")
    print(f"    Expected: NEGATIVE (removing master regulator suppresses ARE)")
else:
    print("  NFE2L2 not in age-DE knockdowns")

print("\n--- Results by Age Direction ---")
for direction in ['UP', 'DOWN']:
    subset = forward_df[forward_df['age_direction'] == direction]
    print(f"\n  {direction} (n={len(subset)}):")
    print(f"    Mean ARE enrichment: {subset['ARE_enrichment'].mean():.4f}")
    print(f"    Top 5 (strongest ARE activation when knocked down):")
    for _, r in subset.head(5).iterrows():
        print(f"      {r['target_gene']}: {r['ARE_enrichment']:.3f}")
    print(f"    Bottom 5 (strongest ARE suppression when knocked down):")
    for _, r in subset.tail(5).iterrows():
        print(f"      {r['target_gene']}: {r['ARE_enrichment']:.3f}")

# Save
forward_df.to_csv(RESULTS_DIR / "forward_age_de_are_gsea.csv", index=False)
print(f"\nSaved: forward_age_de_are_gsea.csv")

#%% ===========================================================================
# ANALYSIS 2: BACKWARD - AGING MIMETICS
# Which knockdowns phenocopy aging?
# ===========================================================================

print("\n" + "=" * 60)
print("[ANALYSIS 2] BACKWARD: Which knockdowns mimic aging?")
print("=" * 60)

print("Running enrichment analysis for ALL knockdowns...")

# For speed, we'll do this for all ~11k knockdowns using the quick method
aging_mimetic_results = []

for idx in range(adata.n_obs):
    if idx % 1000 == 0:
        print(f"  Processing perturbation {idx}/{adata.n_obs}...")
    
    profile = get_knockdown_profile(adata, idx, gene_names)
    target_gene = perturbation_targets[idx]
    
    # Enrichment of age-UP genes at top of knockdown profile
    es_up, n_up = quick_enrichment_score(profile, age_UP)
    
    # Enrichment of age-DOWN genes (we want them at bottom, so negative enrichment means mimetic)
    es_down, n_down = quick_enrichment_score(profile, age_DOWN)
    
    # Combined mimetic score: age-UP enriched at top AND age-DOWN depleted at top
    # Mimetic = high es_up AND low (negative) es_down
    # So: mimetic_score = es_up - es_down
    mimetic_score = es_up - es_down if not (np.isnan(es_up) or np.isnan(es_down)) else np.nan
    
    aging_mimetic_results.append({
        'perturbation_id': perturbation_ids[idx],
        'target_gene': target_gene,
        'age_UP_enrichment': es_up,
        'age_DOWN_enrichment': es_down,
        'mimetic_score': mimetic_score,
        'is_ARE_gene': target_gene in are_genes,
        'is_age_UP': target_gene in age_UP,
        'is_age_DOWN': target_gene in age_DOWN
    })

aging_mim_df = pd.DataFrame(aging_mimetic_results)
aging_mim_df = aging_mim_df.sort_values('mimetic_score', ascending=False)

print("\n--- SANITY CHECKS ---")
print(f"Mimetic score distribution:")
print(f"  Mean: {aging_mim_df['mimetic_score'].mean():.4f}")
print(f"  Std: {aging_mim_df['mimetic_score'].std():.4f}")
print(f"  Min: {aging_mim_df['mimetic_score'].min():.4f}")
print(f"  Max: {aging_mim_df['mimetic_score'].max():.4f}")

print(f"\nAge-UP genes as knockdowns should NOT be strong mimetics (knocking down removes the age effect):")
age_up_kds = aging_mim_df[aging_mim_df['is_age_UP'] == True]
print(f"  Mean mimetic score of age-UP knockdowns: {age_up_kds['mimetic_score'].mean():.4f}")

print("\n--- Top Aging MIMETICS (knockdown phenocopies aging) ---")
top_mimetics = aging_mim_df.head(20)
for _, r in top_mimetics.iterrows():
    are_flag = " [ARE]" if r['is_ARE_gene'] else ""
    print(f"  {r['target_gene']}: {r['mimetic_score']:.3f}{are_flag}")

print("\n--- Top Aging ANTAGONISTS (knockdown reverses aging) ---")
top_antagonists = aging_mim_df.tail(20).iloc[::-1]
for _, r in top_antagonists.iterrows():
    are_flag = " [ARE]" if r['is_ARE_gene'] else ""
    print(f"  {r['target_gene']}: {r['mimetic_score']:.3f}{are_flag}")

# ARE enrichment among mimetics
print("\n--- ARE Enrichment Among Mimetics/Antagonists ---")
top_200 = set(aging_mim_df.head(200)['target_gene'])
bot_200 = set(aging_mim_df.tail(200)['target_gene'])

are_in_top = len(top_200 & are_genes)
are_in_bot = len(bot_200 & are_genes)
print(f"  ARE genes in top 200 mimetics: {are_in_top}/200 = {are_in_top/200*100:.1f}%")
print(f"  ARE genes in top 200 antagonists: {are_in_bot}/200 = {are_in_bot/200*100:.1f}%")
print(f"  Background rate: {len(are_genes & set(perturbation_targets))}/{len(set(perturbation_targets))} = {len(are_genes & set(perturbation_targets))/len(set(perturbation_targets))*100:.1f}%")

# Save
aging_mim_df.to_csv(RESULTS_DIR / "backward_aging_mimetics.csv", index=False)
print(f"\nSaved: backward_aging_mimetics.csv")

#%% ===========================================================================
# ANALYSIS 3: BACKWARD - PRG4 MIMETICS  
# Which knockdowns do what PRG4 does?
# ===========================================================================

print("\n" + "=" * 60)
print("[ANALYSIS 3] BACKWARD: Which knockdowns mimic PRG4?")
print("=" * 60)

print("Running enrichment analysis for ALL knockdowns...")

prg4_mimetic_results = []

for idx in range(adata.n_obs):
    if idx % 1000 == 0:
        print(f"  Processing perturbation {idx}/{adata.n_obs}...")
    
    profile = get_knockdown_profile(adata, idx, gene_names)
    target_gene = perturbation_targets[idx]
    
    # Enrichment of PRG4-UP genes at top
    es_up, n_up = quick_enrichment_score(profile, prg4_UP)
    
    # Enrichment of PRG4-DOWN genes
    es_down, n_down = quick_enrichment_score(profile, prg4_DOWN)
    
    # PRG4 mimetic = PRG4-UP enriched at top AND PRG4-DOWN at bottom
    prg4_mim_score = es_up - es_down if not (np.isnan(es_up) or np.isnan(es_down)) else np.nan
    
    prg4_mimetic_results.append({
        'perturbation_id': perturbation_ids[idx],
        'target_gene': target_gene,
        'PRG4_UP_enrichment': es_up,
        'PRG4_DOWN_enrichment': es_down,
        'prg4_mimetic_score': prg4_mim_score,
        'is_ARE_gene': target_gene in are_genes,
        'is_age_UP': target_gene in age_UP,
        'is_age_DOWN': target_gene in age_DOWN
    })

prg4_mim_df = pd.DataFrame(prg4_mimetic_results)
prg4_mim_df = prg4_mim_df.sort_values('prg4_mimetic_score', ascending=False)

print("\n--- Top PRG4 MIMETICS ---")
for _, r in prg4_mim_df.head(20).iterrows():
    flags = []
    if r['is_ARE_gene']: flags.append("ARE")
    if r['is_age_UP']: flags.append("age-UP")
    if r['is_age_DOWN']: flags.append("age-DOWN")
    flag_str = f" [{', '.join(flags)}]" if flags else ""
    print(f"  {r['target_gene']}: {r['prg4_mimetic_score']:.3f}{flag_str}")

# Save
prg4_mim_df.to_csv(RESULTS_DIR / "backward_prg4_mimetics.csv", index=False)
print(f"\nSaved: backward_prg4_mimetics.csv")

#%% ===========================================================================
# ANALYSIS 4: CONVERGENCE
# Where do aging antagonists, PRG4 mimetics, and ARE overlap?
# ===========================================================================

print("\n" + "=" * 60)
print("[ANALYSIS 4] CONVERGENCE ANALYSIS")
print("=" * 60)

# Define sets
aging_antagonists = set(aging_mim_df.tail(200)['target_gene'])
aging_mimetics = set(aging_mim_df.head(200)['target_gene'])
prg4_mimetics = set(prg4_mim_df.head(200)['target_gene'])
are_set = are_genes & set(perturbation_targets)  # ARE genes that are knockdowns

print("Set sizes:")
print(f"  Aging antagonists (top 200): {len(aging_antagonists)}")
print(f"  PRG4 mimetics (top 200): {len(prg4_mimetics)}")
print(f"  ARE genes (with knockdowns): {len(are_set)}")

print("\n--- Pairwise Overlaps ---")
print(f"  Aging antagonists ∩ PRG4 mimetics: {len(aging_antagonists & prg4_mimetics)}")
print(f"  Aging antagonists ∩ ARE: {len(aging_antagonists & are_set)}")
print(f"  PRG4 mimetics ∩ ARE: {len(prg4_mimetics & are_set)}")

triple = aging_antagonists & prg4_mimetics & are_set
print(f"\n--- Triple Intersection ---")
print(f"  Aging antagonists ∩ PRG4 mimetics ∩ ARE: {len(triple)}")
if len(triple) > 0:
    print(f"  Genes: {sorted(triple)}")

# Correlation between aging antagonist score and PRG4 mimetic score
merged = aging_mim_df.merge(prg4_mim_df[['target_gene', 'prg4_mimetic_score']], on='target_gene')
merged = merged.dropna(subset=['mimetic_score', 'prg4_mimetic_score'])
if len(merged) > 100:
    # Aging antagonist = negative mimetic score, PRG4 mimetic = positive score
    # If PRG4 reverses aging, expect negative correlation
    r, p = stats.spearmanr(-merged['mimetic_score'], merged['prg4_mimetic_score'])
    print(f"\n--- Correlation: Aging Antagonist vs PRG4 Mimetic ---")
    print(f"  Spearman r = {r:.3f}, p = {p:.2e}")
    if r > 0 and p < 0.05:
        print("  INTERPRETATION: Knockdowns that reverse aging also mimic PRG4 ✓")
    elif r < 0 and p < 0.05:
        print("  INTERPRETATION: Knockdowns that reverse aging do NOT mimic PRG4")
    else:
        print("  INTERPRETATION: No significant relationship")

# Save convergence
convergence = {
    'aging_antagonists_in_prg4_mimetics': len(aging_antagonists & prg4_mimetics),
    'aging_antagonists_in_ARE': len(aging_antagonists & are_set),
    'prg4_mimetics_in_ARE': len(prg4_mimetics & are_set),
    'triple_intersection': len(triple),
    'triple_genes': ';'.join(sorted(triple)) if triple else '',
    'antagonist_prg4_correlation_r': r if len(merged) > 100 else np.nan,
    'antagonist_prg4_correlation_p': p if len(merged) > 100 else np.nan
}
pd.DataFrame([convergence]).to_csv(RESULTS_DIR / "convergence_overlaps.csv", index=False)
print(f"\nSaved: convergence_overlaps.csv")

#%% ===========================================================================
# FINAL SUMMARY
# ===========================================================================

print("\n" + "=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)

print("\nOutput files:")
for f in RESULTS_DIR.glob("*.csv"):
    if 'forward' in f.name or 'backward' in f.name or 'convergence' in f.name:
        print(f"  {f.name}")

print("\n--- SUMMARY OF KEY FINDINGS ---")
print(f"1. Forward: {len(forward_df)} age-DE knockdowns analyzed for ARE effect")
print(f"2. Backward-Age: {len(aging_mim_df)} knockdowns scored for aging mimetic/antagonist")
print(f"3. Backward-PRG4: {len(prg4_mim_df)} knockdowns scored for PRG4 mimetic")
print(f"4. Convergence: {len(triple)} genes are aging antagonists, PRG4 mimetics, AND ARE genes")
