import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
BULK_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/RPE_gene pvals.xlsx")
GMT_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/c2.cp.kegg.v7.0.symbols.gmt")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/robustness-analysis")

# ==========================================
# FUNCTIONS
# ==========================================

def load_gmt(gmt_path):
    pathways = {}
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            genes = set(parts[2:]) # Skip URL
            pathways[name] = genes
    return pathways

def load_bulk_genes():
    # We just need the gene symbols and p-values/logfc to recreate the lists
    # Note: The bulk file has 'hgnc_symbol'
    df = pd.read_excel(BULK_FILE)
    return df

def get_gene_list(df, contrast_col_pval, contrast_col_fc, threshold_p=0.05):
    # Filter by nominal p-value
    mask = df[contrast_col_pval] < threshold_p
    # Return symbols
    return set(df.loc[mask, "hgnc_symbol"].dropna().tolist())

def run_enrichment(gene_set, pathways, background_set):
    results = []
    N = len(background_set)
    k = len(gene_set)
    
    for name, pathway_genes in pathways.items():
        # Genes in pathway and in background
        pathway_genes_bg = pathway_genes.intersection(background_set)
        m = len(pathway_genes_bg)
        
        if m < 5: continue # Skip small pathways
        
        # Overlap
        overlap = gene_set.intersection(pathway_genes_bg)
        x = len(overlap)
        
        if x == 0: continue
        
        # Contingency table
        #      InPath  NotInPath
        # InSet   x       k-x
        # NotIn   m-x     (N-k)-(m-x) -> N-k-m+x
        
        table = [[x, k - x], [m - x, N - k - m + x]]
        odds, pval = stats.fisher_exact(table, alternative='greater')
        
        results.append({
            "pathway": name,
            "pvalue": pval,
            "overlap_count": x,
            "pathway_size": m,
            "overlap_genes": ",".join(list(overlap))
        })
        
    res_df = pd.DataFrame(results)
    if not res_df.empty:
        res_df = res_df.sort_values("pvalue")
    return res_df

def plot_top_pathways(res_df, title, filename):
    if res_df.empty: return
    
    top = res_df.head(10).iloc[::-1] # Top 10, reversed for plotting
    
    plt.figure(figsize=(10, 6))
    plt.barh(top["pathway"], -np.log10(top["pvalue"] ), color='skyblue')
    plt.xlabel("-Log10 P-value")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# ==========================================
# MAIN
# ==========================================

def main():
    print("Loading data...")
    pathways = load_gmt(GMT_FILE)
    df = load_bulk_genes()
    
    # Define Background: All genes in Bulk that have a symbol
    background = set(df["hgnc_symbol"].dropna().tolist())
    print(f"Background size: {len(background)}")
    
    # Analysis Config
    configs = [
        {"name": "H2O2_p05", "col_p": "H2O2_vs_CTRL.pvalue", "col_fc": "H2O2_vs_CTRL.log2FoldChange", "thresh": 0.05},
        {"name": "H2O2_p01", "col_p": "H2O2_vs_CTRL.pvalue", "col_fc": "H2O2_vs_CTRL.log2FoldChange", "thresh": 0.01},
        {"name": "Rescue_p05", "col_p": "H2O2PRG4_vs_H2O2.pvalue", "col_fc": "H2O2PRG4_vs_H2O2.log2FoldChange", "thresh": 0.05},
        {"name": "Rescue_p01", "col_p": "H2O2PRG4_vs_H2O2.pvalue", "col_fc": "H2O2PRG4_vs_H2O2.log2FoldChange", "thresh": 0.01},
    ]
    
    all_res = []
    
    for cfg in configs:
        print(f"Running enrichment for {cfg['name']}...")
        gene_list = get_gene_list(df, cfg["col_p"], cfg["col_fc"], cfg["thresh"])
        # Intersect with background (should be subset anyway but safe to do)
        gene_list = gene_list.intersection(background)
        
        res = run_enrichment(gene_list, pathways, background)
        
        if not res.empty:
            res["signature"] = cfg["name"]
            all_res.append(res.head(20)) # Keep top 20 for summary CSV
            
            # Plot
            plot_path = os.path.join(OUTPUT_DIR, f"enrichment_{cfg['name']}.png")
            plot_top_pathways(res, f"Top Pathways: {cfg['name']}", plot_path)
    
    # Save Summary
    if all_res:
        final_df = pd.concat(all_res)
        final_df.to_csv(os.path.join(OUTPUT_DIR, "signature_enrichment_summary.csv"), index=False)
        print(f"Saved summary to {os.path.join(OUTPUT_DIR, 'signature_enrichment_summary.csv')}")

if __name__ == "__main__":
    main()
