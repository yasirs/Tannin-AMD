import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# Paths
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
BULK_DIR = os.path.join(BASE_DIR, "data/RPE_cells/code")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/bridge-results")

os.makedirs(OUTPUT_DIR, exist_ok=True)
sns.set_theme(style="whitegrid", context="talk")

def plot_volcano(file_path, lfc_col, pval_col, label, top_n=10):
    df = pd.read_excel(file_path)
    
    # Filter for valid values
    df = df.dropna(subset=[lfc_col, pval_col])
    df['-log10p'] = -np.log10(df[pval_col] + 1e-300)
    
    plt.figure(figsize=(10, 8))
    
    # Define significance
    sig = (df[pval_col] < 0.05) & (df[lfc_col].abs() > 1)
    
    sns.scatterplot(data=df[~sig], x=lfc_col, y="-log10p", alpha=0.3, color="gray", s=10)
    sns.scatterplot(data=df[sig], x=lfc_col, y="-log10p", alpha=0.6, color="red", s=15)
    
    # Annotate top genes
    top_genes = df.sort_values(pval_col).head(top_n)
    for i, row in top_genes.iterrows():
        plt.text(row[lfc_col], row['-log10p'], row['hgnc_symbol'], fontsize=10)
        
    plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)
    plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    
    plt.title(f"Volcano Plot: {label}")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 P-value")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"Volcano_{label.replace(' ', '_')}.png"), dpi=300)
    plt.close()

def plot_pca(tpm_file):
    df = pd.read_excel(tpm_file)
    # TPM file has genes in rows, samples in columns?
    # Let's assume columns are TS1...TS16
    sample_cols = [c for c in df.columns if c.startswith("TS")]
    data = df[sample_cols].values.T # samples in rows for PCA
    
    pca = PCA(n_components=2)
    components = pca.fit_transform(np.log2(data + 1))
    
    # Sample mapping
    groups = {
        "CTRL": ["TS1", "TS2", "TS3"],
        "H2O2": ["TS5", "TS6", "TS7"],
        "PRG4": ["TS9", "TS11", "TS12"],
        "H2O2PRG4": ["TS13", "TS15", "TS16"]
    }
    sample_to_group = {}
    for g, ss in groups.items():
        for s in ss:
            sample_to_group[s] = g
            
    plot_df = pd.DataFrame({
        "PC1": components[:, 0],
        "PC2": components[:, 1],
        "Sample": sample_cols,
        "Condition": [sample_to_group.get(s, "Other") for s in sample_cols]
    })
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=plot_df, x="PC1", y="PC2", hue="Condition", s=200, style="Condition")
    
    # Add labels to points
    for i, row in plot_df.iterrows():
        plt.text(row['PC1']+0.1, row['PC2']+0.1, row['Sample'], fontsize=12)
        
    plt.title("PCA of RPE Bulk RNA-seq Samples")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "PCA_RPE_Bulk.png"), dpi=300)
    plt.close()

if __name__ == "__main__":
    pvals_file = os.path.join(BULK_DIR, "RPE_gene pvals.xlsx")
    tpm_file = os.path.join(BULK_DIR, "RPE_TPMS.xlsx")
    
    print("Generating PCA plot...")
    plot_pca(tpm_file)
    
    print("Generating Volcano plot for H2O2 vs CTRL...")
    plot_volcano(pvals_file, "H2O2_vs_CTRL.log2FoldChange", "H2O2_vs_CTRL.pvalue", "H2O2 vs CTRL")
    
    print("Generating Volcano plot for PRG4 Rescue...")
    plot_volcano(pvals_file, "H2O2PRG4_vs_H2O2.log2FoldChange", "H2O2PRG4_vs_H2O2.pvalue", "PRG4 Rescue")
    
    print(f"Bulk visualizations saved to {OUTPUT_DIR}")
