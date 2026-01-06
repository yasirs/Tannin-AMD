import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
BULK_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/RPE_gene pvals.xlsx")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/external-validation")

# Known AMD Risk / Pathogenic Genes (GWAS & Literature)
AMD_GENES = {
    "Complement": ["CFH", "C3", "CFI", "C2", "CFB", "CFD"],
    "Oxidative/Mito": ["ARMS2", "HTRA1", "SOD2", "MT-ND2"],
    "Lipid/ECM": ["APOE", "ABCA4", "TIMP3", "LIPC"],
    "RPE_Function": ["BEST1", "RPE65", "RLBP1"]
}

# ==========================================
# MAIN
# ==========================================

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    print("Loading bulk results...")
    df = pd.read_excel(BULK_FILE)
    
    # Check for presence of AMD genes
    flat_genes = [g for sublist in AMD_GENES.values() for g in sublist]
    
    # Filter
    sub = df[df["hgnc_symbol"].isin(flat_genes)].copy()
    
    if sub.empty:
        print("No AMD genes found in dataset.")
        return

    # Extract LogFCs
    # H2O2 vs CTRL (Disease Model)
    # PRG4 vs CTRL (Drug Baseline)
    # Rescue vs H2O2 (Drug Effect)
    
    print("\nAMD Gene Behavior:")
    print(f"{ 'Gene':<10} | { 'H2O2_FC':<10} | { 'Rescue_FC':<10} | {'Category'}")
    print("-" * 50)
    
    res_rows = []
    
    for cat, genes in AMD_GENES.items():
        for g in genes:
            row = sub[sub["hgnc_symbol"] == g]
            if not row.empty:
                fc_dis = row.iloc[0]["H2O2_vs_CTRL.log2FoldChange"]
                pval_dis = row.iloc[0]["H2O2_vs_CTRL.pvalue"]
                fc_res = row.iloc[0]["H2O2PRG4_vs_H2O2.log2FoldChange"]
                pval_res = row.iloc[0]["H2O2PRG4_vs_H2O2.pvalue"]
                
                print(f"{g:<10} | {fc_dis:.3f}      | {fc_res:.3f}      | {cat}")
                
                res_rows.append({
                    "Gene": g,
                    "Category": cat,
                    "Disease_LogFC": fc_dis,
                    "Disease_P": pval_dis,
                    "Rescue_LogFC": fc_res,
                    "Rescue_P": pval_res,
                    "Reversed": (np.sign(fc_dis) != np.sign(fc_res)) and (abs(fc_dis) > 0.1)
                })
    
    # Save Results
    res_df = pd.DataFrame(res_rows)
    res_df.to_csv(os.path.join(OUTPUT_DIR, "amd_risk_gene_validation.csv"), index=False)
    
    # Plotting
    # Scatter: Disease vs Rescue
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=res_df, x="Disease_LogFC", y="Rescue_LogFC", hue="Category", s=100)
    plt.axhline(0, color='gray', linestyle='--')
    plt.axvline(0, color='gray', linestyle='--')
    # Add labels
    for i, row in res_df.iterrows():
        plt.text(row["Disease_LogFC"], row["Rescue_LogFC"], row["Gene"], fontsize=9)
        
    plt.title("PRG4 Rescue Effect on AMD Risk Genes")
    plt.xlabel("H2O2 Stress Log2FC (Disease Model)")
    plt.ylabel("PRG4 Rescue Log2FC (Treatment)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "amd_rescue_scatter.png"))
    
    # Correlation of Rescue
    # If PRG4 works, we expect negative correlation (Rescue opposes Disease)
    if len(res_df) > 2:
        corr = res_df["Disease_LogFC"].corr(res_df["Rescue_LogFC"])
        print(f"\nCorrelation between Disease and Rescue LFC for Risk Genes: {corr:.3f}")

if __name__ == "__main__":
    main()
