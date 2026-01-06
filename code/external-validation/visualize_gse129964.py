import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
DATA_FILE = os.path.join(BASE_DIR, "data/external/GSE129964/GSE129964_serumCountsTable.tsv.gz")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/external-validation")

# Set visualization style
sns.set_theme(style="whitegrid", context="paper")
plt.rcParams['font.family'] = 'serif'

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    print("Loading GSE129964 data...")
    df = pd.read_csv(DATA_FILE, sep='\t', index_col=0)
    
    # CPM Normalization
    libsizes = df.sum(axis=0)
    cpm = df.div(libsizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)
    
    # Metadata Parsing
    # Columns like C.Day0.Repeat1, S.Day1.Repeat1
    # C = Control, S = Starved
    samples = log_cpm.columns
    metadata = []
    for s in samples:
        parts = s.split('.')
        metadata.append({
            "Sample": s,
            "Condition": "Control" if parts[0] == "C" else "Starved",
            "Day": int(parts[1].replace("Day", "")),
            "Repeat": parts[2]
        })
    meta_df = pd.DataFrame(metadata)
    
    # 1. PCA
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    print("Computing PCA...")
    vars = log_cpm.var(axis=1).sort_values(ascending=False)
    top_genes = vars.head(2000).index
    
    pca = PCA(n_components=2)
    scaled_data = StandardScaler().fit_transform(log_cpm.loc[top_genes].T)
    pca_coords = pca.fit_transform(scaled_data)
    
    pca_plot_df = pd.DataFrame(pca_coords, columns=["PC1", "PC2"], index=samples)
    pca_plot_df = pca_plot_df.merge(meta_df, left_index=True, right_on="Sample")
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_plot_df, x="PC1", y="PC2", hue="Day", style="Condition", s=100, palette="viridis")
    plt.title("PCA of GSE129964 (Serum Starvation Timecourse)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "gse129964_pca.pdf"))
    plt.close()
    
    # 2. Key Genes Timecourse
    print("Plotting Timecourse...")
    # Map symbols? The file uses Symbols already (A1BG, etc.)
    risk_genes = ["CFH", "RPE65", "VEGFA", "LDHA"]
    
    tidy_data = []
    for g in risk_genes:
        if g in log_cpm.index:
            expr = log_cpm.loc[g]
            for s, val in expr.items():
                m = meta_df[meta_df["Sample"] == s].iloc[0]
                tidy_data.append({
                    "Gene": g,
                    "Expression": val,
                    "Day": m["Day"],
                    "Condition": m["Condition"]
                })
    
    if tidy_data:
        tidy_df = pd.DataFrame(tidy_data)
        plt.figure(figsize=(10, 6))
        sns.lineplot(data=tidy_df, x="Day", y="Expression", hue="Gene", style="Condition", markers=True, dashes=False)
        plt.title("Expression Timecourse under Serum Starvation")
        plt.savefig(os.path.join(OUTPUT_DIR, "gse129964_timecourse.pdf"))
        plt.close()

if __name__ == "__main__":
    main()
