import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
BULK_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/RPE_gene pvals.xlsx")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/sc-analysis")

# Literature-derived markers for Human RPE/Choroid (Voigt et al., etc.)
CELL_TYPE_MARKERS = {
    "RPE": ["RPE65", "BEST1", "RLBP1", "TYR", "MITF", "OTX2", "TTR", "PMEL", "SERPINF1", "CRALBP"],
    "Endothelial": ["PECAM1", "VWF", "CD34", "PLVAP", "ENG", "CLDN5", "ICAM2", "FLT1", "KDR"],
    "Macrophage/Microglia": ["CD14", "AIF1", "CD68", "C1QA", "C1QB", "C1QC", "CSF1R", "CD163", "TYROBP"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "PDGFRB", "LUM", "FBLN1", "COL3A1"],
    "T-Cell": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "B-Cell": ["CD79A", "MS4A1", "CD19"],
    "Melanocyte": ["MLANA", "DCT", "TYRP1"] # Choroidal melanocytes
}

# ==========================================
# FUNCTIONS
# ==========================================

def load_data():
    print("Loading bulk data...")
    df = pd.read_excel(BULK_FILE)
    return df

def run_enrichment(gene_set, markers, background_size=20000):
    results = []
    k = len(gene_set)
    N = background_size
    
    for cell_type, marker_genes in markers.items():
        m = len(marker_genes)
        overlap = gene_set.intersection(set(marker_genes))
        x = len(overlap)
        
        # Fisher's exact
        table = [[x, k - x], [m - x, N - k - m + x]]
        odds, pval = stats.fisher_exact(table, alternative='greater')
        
        results.append({
            "Cell_Type": cell_type,
            "P_Value": pval,
            "Overlap_Count": x,
            "Marker_Set_Size": m,
            "Overlap_Genes": ",".join(list(overlap))
        })
        
    return pd.DataFrame(results).sort_values("P_Value")

# ==========================================
# MAIN
# ==========================================

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    df = load_data()
    
    # Define Signatures
    # H2O2 Stress: Genes UP in H2O2 vs CTRL (p < 0.05)
    stress_up = set(df[(df["H2O2_vs_CTRL.log2FoldChange"] > 0) & (df["H2O2_vs_CTRL.pvalue"] < 0.05)]["hgnc_symbol"].dropna())
    stress_down = set(df[(df["H2O2_vs_CTRL.log2FoldChange"] < 0) & (df["H2O2_vs_CTRL.pvalue"] < 0.05)]["hgnc_symbol"].dropna())
    
    # PRG4 Rescue: Genes UP in Rescue vs H2O2 (p < 0.05) -> Restored genes
    rescue_up = set(df[(df["H2O2PRG4_vs_H2O2.log2FoldChange"] > 0) & (df["H2O2PRG4_vs_H2O2.pvalue"] < 0.05)]["hgnc_symbol"].dropna())
    # Rescue Down: Genes suppressed by PRG4 (Counter-acted stress genes)
    rescue_down = set(df[(df["H2O2PRG4_vs_H2O2.log2FoldChange"] < 0) & (df["H2O2PRG4_vs_H2O2.pvalue"] < 0.05)]["hgnc_symbol"].dropna())
    
    signatures = {
        "Stress_UP": stress_up,
        "Stress_DOWN": stress_down,
        "Rescue_UP": rescue_up,
        "Rescue_DOWN": rescue_down
    }
    
    all_res = []
    
    print("\nRunning Cell Type Enrichment...")
    for sig_name, genes in signatures.items():
        print(f"  Signature: {sig_name} (n={len(genes)})")
        res = run_enrichment(genes, CELL_TYPE_MARKERS)
        res["Signature"] = sig_name
        all_res.append(res)
        print(res.head(3)[["Cell_Type", "P_Value", "Overlap_Genes"]])
        
    final_df = pd.concat(all_res)
    final_df.to_csv(os.path.join(OUTPUT_DIR, "cell_type_enrichment.csv"), index=False)
    
    # Check PRG4 and CD44 expression
    print("\nTarget Gene Expression (Bulk):")
    targets = ["PRG4", "CD44", "TLR2", "TLR4"] # PRG4 receptors
    sub = df[df["hgnc_symbol"].isin(targets)][["hgnc_symbol", "H2O2_vs_CTRL.log2FoldChange", "H2O2_vs_CTRL.pvalue", "H2O2PRG4_vs_H2O2.log2FoldChange", "H2O2PRG4_vs_H2O2.pvalue"]]
    print(sub)
    sub.to_csv(os.path.join(OUTPUT_DIR, "prg4_receptor_bulk.csv"), index=False)
    
    # Visualization
    # Heatmap of -log10 P-values
    pvals = final_df.pivot(index="Cell_Type", columns="Signature", values="P_Value")
    pvals = -np.log10(pvals + 1e-10) # Log transform
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(pvals, cmap="Reds", annot=True, fmt=".1f")
    plt.title("Cell Type Marker Enrichment in Bulk Signatures")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "cell_type_enrichment_heatmap.png"))

if __name__ == "__main__":
    main()
