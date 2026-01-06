import os
import pandas as pd
import numpy as np
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
DATA_DIR = os.path.join(BASE_DIR, "data/external/geo")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/cohort-GSE135092") # Saving here for now or separate? 
# Let's create a new folder for GSE29801 or keep in cohort folder but named distinctively.
# User asked for "another public dataset". GSE29801 is it.
# Maybe "results/cohort-GSE29801"?
OUTPUT_DIR = os.path.join(BASE_DIR, "results/cohort-GSE29801")

GPL_PATH = os.path.join(DATA_DIR, "GSE29801_extracted/GPL4133_old_annotations.txt.gz")
EXPR_PATH = os.path.join(DATA_DIR, "GSE29801_expression_matrix.csv.gz")
SIG_PATH = os.path.join(DATA_DIR, "GSE29801_amd_signature.csv")

# Set visualization style
sns.set_theme(style="whitegrid", context="paper")
plt.rcParams['font.family'] = 'serif'

# ==========================================
# FUNCTIONS
# ==========================================

def parse_gpl(gpl_path):
    print("Parsing GPL annotation...")
    id_to_sym = {}
    with gzip.open(gpl_path, 'rt', encoding='latin-1') as f:
        header_found = False
        headers = []
        for line in f:
            line = line.strip()
            if not header_found:
                if line.startswith("ID") and "GENE_SYMBOL" in line:
                    header_found = True
                    headers = line.split('\t')
                continue
            
            parts = line.split('\t')
            if len(parts) < len(headers): continue
            
            # Find indices
            try:
                id_idx = headers.index("ID")
                sym_idx = headers.index("GENE_SYMBOL")
                
                probe_id = parts[id_idx]
                sym = parts[sym_idx]
                
                if sym and sym != "nan":
                    id_to_sym[probe_id] = sym
            except:
                continue
                
    print(f"Mapped {len(id_to_sym)} probes to symbols.")
    return id_to_sym

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    # 1. Load Data
    print("Loading expression matrix...")
    expr = pd.read_csv(EXPR_PATH, compression='gzip', index_col=0)
    
    # 2. Load GPL
    id_map = parse_gpl(GPL_PATH)
    
    # 3. Map IDs to Symbols
    # Expression matrix columns are IDs (except 'Group')
    gene_cols = [c for c in expr.columns if c != "Group"]
    
    # Create a mapping for columns
    col_map = {c: id_map.get(c, c) for c in gene_cols}
    
    # Rename
    expr_mapped = expr.rename(columns=col_map)
    
    # Average duplicate symbols
    # Separate numeric data and metadata
    groups = expr_mapped["Group"]
    data = expr_mapped.drop(columns=["Group"])
    
    # Transpose to (Genes x Samples) for groupby
    data_t = data.T
    data_t.index.name = "Symbol"
    
    # Group by Symbol and mean
    print("Collapsing duplicates...")
    data_collapsed = data_t.groupby("Symbol").mean()
    
    # Transpose back
    expr_final = data_collapsed.T
    expr_final["Group"] = groups
    
    # 4. PCA
    print("Computing PCA...")
    # Filter numeric columns
    numeric_data = expr_final.drop(columns=["Group"])
    
    # Top variable genes
    vars = numeric_data.var().sort_values(ascending=False)
    top_genes = vars.head(2000).index
    
    pca = PCA(n_components=2)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(numeric_data[top_genes])
    pca_coords = pca.fit_transform(scaled_data)
    
    pca_df = pd.DataFrame(pca_coords, columns=["PC1", "PC2"], index=numeric_data.index)
    pca_df["Group"] = groups
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Group", alpha=0.7, palette={"AMD": "firebrick", "Control": "navy"})
    plt.title(f"PCA of GSE29801 (Top 2000 Genes)\nExplained Var: {pca.explained_variance_ratio_[0]:.1%} / {pca.explained_variance_ratio_[1]:.1%}")
    plt.savefig(os.path.join(OUTPUT_DIR, "gse29801_pca.pdf"))
    plt.close()
    
    # 5. Load Signature for Volcano
    print("Plotting Volcano...")
    sig = pd.read_csv(SIG_PATH, index_col=0)
    # Map IDs in signature
    sig["Symbol"] = sig.index.astype(str).map(id_map)
    
    # Calculate -log10 p
    sig["neg_log_p"] = -np.log10(sig["pvalue"])
    sig["Significant"] = (sig["pvalue"] < 0.05) & (sig["logfc_amd_vs_ctrl"].abs() > 0.5)
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=sig, x="logfc_amd_vs_ctrl", y="neg_log_p", hue="Significant", 
                    palette={True: "red", False: "grey"}, alpha=0.5, legend=False, s=10)
    
    # Label Risk Genes
    risk_genes = ["CFH", "ARMS2", "HTRA1", "C3", "RPE65", "TIMP3"]
    
    for sym in risk_genes:
        # Find row with this symbol
        row = sig[sig["Symbol"] == sym]
        if not row.empty:
            # Pick best p-value if multiple
            best_row = row.loc[row["pvalue"].idxmin()]
            plt.text(best_row["logfc_amd_vs_ctrl"], best_row["neg_log_p"], sym, fontsize=9, fontweight='bold')
            plt.scatter(best_row["logfc_amd_vs_ctrl"], best_row["neg_log_p"], color='black', s=20)
            
    plt.title("Volcano Plot: GSE29801 (AMD vs Control)")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 P-value")
    plt.savefig(os.path.join(OUTPUT_DIR, "gse29801_volcano.pdf"))
    plt.close()
    
    # 6. Risk Gene Boxplots
    print("Plotting Risk Genes...")
    tidy_data = []
    
    for sym in ["CFH", "RPE65", "HTRA1", "BEST1"]:
        if sym in expr_final.columns:
            vals = expr_final[sym]
            for i, val in enumerate(vals):
                tidy_data.append({
                    "Gene": sym,
                    "Expression": val,
                    "Group": groups.iloc[i]
                })
                
    if tidy_data:
        tidy_df = pd.DataFrame(tidy_data)
        plt.figure(figsize=(10, 6))
        sns.violinplot(data=tidy_df, x="Gene", y="Expression", hue="Group", split=True, inner="quart", palette={"AMD": "firebrick", "Control": "skyblue"})
        plt.title("Expression of Key AMD Risk Genes (GSE29801)")
        plt.savefig(os.path.join(OUTPUT_DIR, "gse29801_risk_violins.pdf"))
        plt.close()

if __name__ == "__main__":
    main()
