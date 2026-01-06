import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
GSE_FILE = os.path.join(BASE_DIR, "results/external-validation/gse129964_comparison.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/gsea-analysis")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# ==========================================
# GSEA ALGORITHM
# ==========================================

def calc_es(ranked_genes, gene_set, weight=1.0):
    """
    Calculate Enrichment Score (ES) and Leading Edge.
    ranked_genes: list of (gene, score) sorted by score descending.
    gene_set: set of genes.
    """
    N = len(ranked_genes)
    Nh = len(gene_set)
    if Nh == 0: return 0, []
    
    Nm = N - Nh
    
    # Extract scores for genes in set
    # Create a vector of hits (1 if in set, 0 else)
    # And a vector of scores
    
    hits = np.array([1 if g[0] in gene_set else 0 for g in ranked_genes])
    scores = np.array([abs(g[1])**weight for g in ranked_genes])
    
    # P_hit: sum of scores for hits
    hit_score_sum = np.sum(scores[hits == 1])
    if hit_score_sum == 0: return 0, []
    
    P_hit = np.cumsum(hits * scores) / hit_score_sum
    P_miss = np.cumsum(1 - hits) / Nm
    
    RES = P_hit - P_miss
    
    max_idx = np.argmax(np.abs(RES))
    ES = RES[max_idx]
    
    # Leading Edge
    # If ES > 0, leading edge is 0 to max_idx
    # If ES < 0, leading edge is max_idx to end (conceptually)
    
    if ES >= 0:
        le_idx = max_idx
        leading_edge_genes = [ranked_genes[i][0] for i in range(le_idx + 1) if ranked_genes[i][0] in gene_set]
    else:
        # For negative ES, the peak is the minimum (most negative deviation)
        # Leading edge is usually considered the genes *after* the peak if we look at the tail
        # But standard GSEA defines it relative to the peak.
        # Let's simplify: Return genes contributing to the peak.
        le_idx = max_idx
        leading_edge_genes = [ranked_genes[i][0] for i in range(le_idx, N) if ranked_genes[i][0] in gene_set]

    return ES, leading_edge_genes, RES

def permutation_test(ranked_genes, gene_set, observed_es, n_perm=1000):
    null_es = []
    genes = [g[0] for g in ranked_genes]
    scores = [g[1] for g in ranked_genes]
    N = len(ranked_genes)
    
    for _ in range(n_perm):
        # Shuffle genes
        shuffled_genes = np.random.permutation(genes)
        # Zip with original scores
        shuffled_ranked = list(zip(shuffled_genes, scores))
        es, _, _ = calc_es(shuffled_ranked, gene_set)
        null_es.append(es)
        
    null_es = np.array(null_es)
    if observed_es > 0:
        pval = np.sum(null_es >= observed_es) / n_perm
    else:
        pval = np.sum(null_es <= observed_es) / n_perm
        
    return pval, np.mean(null_es), np.std(null_es)

# ==========================================
# MAIN
# ==========================================

def main():
    print("Loading data...")
    df = pd.read_csv(GSE_FILE)
    
    # Prepare Ranked List (Serum Starvation LogFC)
    # Sort descending
    df_rank = df.sort_values("logfc_serum_stress", ascending=False)
    ranked_list = list(zip(df_rank["gene_symbol"], df_rank["logfc_serum_stress"]))
    
    # Prepare Gene Sets (PRG4 Rescue UP/DOWN)
    # Use significant DEGs
    # Rescue UP: logFC > 0 and p < 0.05
    rescue_up = set(df[(df["H2O2PRG4_vs_H2O2.log2FoldChange"] > 0) & (df["H2O2PRG4_vs_H2O2.pvalue"] < 0.05)]["gene_symbol"])
    rescue_down = set(df[(df["H2O2PRG4_vs_H2O2.log2FoldChange"] < 0) & (df["H2O2PRG4_vs_H2O2.pvalue"] < 0.05)]["gene_symbol"])
    
    print(f"Ranked List Size: {len(ranked_list)}")
    print(f"Gene Set Rescue_UP Size: {len(rescue_up)}")
    print(f"Gene Set Rescue_DOWN Size: {len(rescue_down)}")
    
    sets = {
        "PRG4_Rescue_UP": rescue_up,
        "PRG4_Rescue_DOWN": rescue_down
    }
    
    results = []
    
    for name, gset in sets.items():
        print(f"\nRunning GSEA for {name}...")
        common = gset.intersection(set(df_rank["gene_symbol"]))
        print(f"  Overlap with ranked list: {len(common)}")
        
        if len(common) < 10:
            print("  Skipping (too few genes)")
            continue
            
        es, le_genes, res_curve = calc_es(ranked_list, common)
        print(f"  ES: {es:.4f}")
        
        pval, _, _ = permutation_test(ranked_list, common, es, n_perm=1000)
        print(f"  P-value: {pval:.4f}")
        
        results.append({
            "GeneSet": name,
            "ES": es,
            "Pval": pval,
            "Leading_Edge_Count": len(le_genes),
            "Leading_Edge_Genes": ",".join(le_genes)
        })
        
        # Plot
        plt.figure(figsize=(10, 4))
        plt.plot(res_curve, label=f"ES={es:.3f}, p={pval:.3f}", color="green" if es > 0 else "red")
        plt.axhline(0, color='black', linewidth=0.5)
        plt.title(f"GSEA: {name} vs Serum Starvation")
        plt.xlabel("Rank in Serum Starvation (Left=UP, Right=DOWN)")
        plt.ylabel("Enrichment Score")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f"gsea_plot_{name}.png"))
        
        # Save Leading Edge
        with open(os.path.join(OUTPUT_DIR, f"leading_edge_{name}.txt"), "w") as f:
            for g in le_genes:
                f.write(g + "\n")

    # Save Summary
    pd.DataFrame(results).to_csv(os.path.join(OUTPUT_DIR, "gsea_summary.csv"), index=False)

if __name__ == "__main__":
    main()
