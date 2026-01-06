import os
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.stats as stats

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
LEADING_EDGE_FILE = os.path.join(BASE_DIR, "results/gsea-analysis/leading_edge_PRG4_Rescue_UP.txt")
PERTURB_FILE = os.path.join(BASE_DIR, "data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad")
GMT_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/c2.cp.kegg.v7.0.symbols.gmt")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/gsea-analysis")

# Top Mimetics from Virtual Screen to check
MIMETICS = ["KEAP1", "DICER1", "UBR5", "SHOC2"]

# ==========================================
# FUNCTIONS
# ==========================================

def load_gmt(gmt_path):
    pathways = {}
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            genes = set(parts[2:])
            pathways[name] = genes
    return pathways

def run_enrichment(gene_list, pathways, background_size=20000):
    results = []
    gene_set = set(gene_list)
    k = len(gene_set)
    N = background_size
    
    for name, pathway_genes in pathways.items():
        m = len(pathway_genes)
        overlap = gene_set.intersection(pathway_genes)
        x = len(overlap)
        
        if x == 0: continue
        
        table = [[x, k - x], [m - x, N - k - m + x]]
        odds, pval = stats.fisher_exact(table, alternative='greater')
        
        results.append({
            "pathway": name,
            "pvalue": pval,
            "overlap_count": x,
            "overlap_genes": ",".join(list(overlap))
        })
        
    return pd.DataFrame(results).sort_values("pvalue")

def get_perturbation_effect(perturb_file, target_genes, leading_edge):
    print(f"Loading K562 GWPS to check {target_genes}...")
    ad = sc.read_h5ad(perturb_file)
    
    # Find observations for targets
    # Obs names are [ID]_[Target]...
    
    effects = {}
    
    for target in target_genes:
        # Regex or simple string match
        # obs_names contains target gene symbol usually
        mask = ad.obs_names.str.contains(f"_{target}_")
        if mask.sum() == 0:
            print(f"  Target {target} not found in Perturb-seq.")
            continue
            
        print(f"  Found {mask.sum()} profiles for {target}.")
        
        # Get mean profile
        X_sub = ad[mask].X
        if hasattr(X_sub, "toarray"):
            X_sub = X_sub.toarray()
        profile = np.mean(X_sub, axis=0)
        
        # Create Series
        prof_series = pd.Series(profile, index=ad.var_names)
        
        # Filter for Leading Edge genes (assuming var_names are IDs, need mapping?)
        # ad.var usually has gene_name column.
        
        if "gene_name" in ad.var.columns:
            prof_series.index = ad.var["gene_name"]
        
        # Intersection
        common_le = list(set(leading_edge) & set(prof_series.index))
        
        if not common_le:
            continue
            
        # Scores of LE genes in this profile
        le_scores = prof_series.loc[common_le]
        
        # Do they tend to be UP? (Since Leading Edge is Rescue UP)
        mean_le = le_scores.mean()
        pct_up = (le_scores > 0).mean()
        
        # Enrichment of UP regulation?
        # Compare LE scores vs Background scores
        bg_scores = prof_series
        stat, pval = stats.mannwhitneyu(le_scores, bg_scores, alternative='greater')
        
        effects[target] = {
            "Mean_LFC_in_LE": mean_le,
            "Pct_UP_in_LE": pct_up,
            "MW_Pval": pval
        }
        
    return pd.DataFrame(effects).T

# ==========================================
# MAIN
# ==========================================

def main():
    # 1. Load Leading Edge
    with open(LEADING_EDGE_FILE, 'r') as f:
        le_genes = [line.strip() for line in f]
    print(f"Loaded {len(le_genes)} Leading Edge genes.")
    
    # 2. Enrichment
    print("Running Pathway Enrichment on Leading Edge...")
    pathways = load_gmt(GMT_FILE)
    enr = run_enrichment(le_genes, pathways)
    if not enr.empty:
        print(enr.head(10)[["pathway", "pvalue", "overlap_genes"]])
        enr.to_csv(os.path.join(OUTPUT_DIR, "leading_edge_enrichment.csv"), index=False)
        
    # 3. Check Mimetics
    eff = get_perturbation_effect(PERTURB_FILE, MIMETICS, le_genes)
    if not eff.empty:
        print("\nDo Mimetics induce the Leading Edge genes?")
        print(eff)
        eff.to_csv(os.path.join(OUTPUT_DIR, "leading_edge_mimetic_validation.csv"))

if __name__ == "__main__":
    main()
