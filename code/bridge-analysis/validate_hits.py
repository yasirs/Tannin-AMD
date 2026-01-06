import os
import pandas as pd
import numpy as np
import scipy.stats as stats

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
SCREEN_FILE = os.path.join(BASE_DIR, "results/virtual-screen/prg4_virtual_screen_results.csv")
GMT_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/c2.cp.kegg.v7.0.symbols.gmt")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/virtual-screen")

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
        
        # Hypergeometric test
        table = [[x, k - x], [m - x, N - k - m + x]]
        odds, pval = stats.fisher_exact(table, alternative='greater')
        
        results.append({
            "pathway": name,
            "pvalue": pval,
            "overlap_count": x,
            "overlap_genes": ",".join(list(overlap))
        })
        
    res_df = pd.DataFrame(results)
    if not res_df.empty:
        res_df = res_df.sort_values("pvalue")
    return res_df

# ==========================================
# MAIN
# ==========================================

def main():
    print("Loading screen results...")
    df = pd.read_csv(SCREEN_FILE)
    pathways = load_gmt(GMT_FILE)
    
    # Top 100 Mimetics
    top_genes = df.head(100)["gene_symbol"].tolist()
    print(f"Enrichment for Top 100 Mimetics (e.g., {top_genes[:3]})...")
    res_top = run_enrichment(top_genes, pathways)
    if not res_top.empty:
        print(res_top.head(5)[["pathway", "pvalue", "overlap_genes"]])
        res_top.to_csv(os.path.join(OUTPUT_DIR, "enrichment_top_mimetics.csv"), index=False)
        
    # Top 100 Antagonists (Bottom of list)
    bot_genes = df.tail(100)["gene_symbol"].tolist()
    print(f"\nEnrichment for Top 100 Antagonists (e.g., {bot_genes[-3:]})...")
    res_bot = run_enrichment(bot_genes, pathways)
    if not res_bot.empty:
        print(res_bot.head(5)[["pathway", "pvalue", "overlap_genes"]])
        res_bot.to_csv(os.path.join(OUTPUT_DIR, "enrichment_top_antagonists.csv"), index=False)

if __name__ == "__main__":
    main()
