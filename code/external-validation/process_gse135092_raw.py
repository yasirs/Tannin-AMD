import os
import pandas as pd
import numpy as np
import gzip
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
RAW_DIR = os.path.join(BASE_DIR, "data/external/geo/GSE135092_extracted")
META_FILE = os.path.join(BASE_DIR, "data/external/geo/GSE135092_series_matrix.txt.gz")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/cohort-GSE135092")

# Set visualization style
sns.set_theme(style="whitegrid", context="paper")
plt.rcParams['font.family'] = 'serif'

# ==========================================
# FUNCTIONS
# ==========================================

def parse_metadata(meta_path):
    print("Parsing metadata from Series Matrix...")
    metadata = {}
    with gzip.open(meta_path, 'rt', encoding='latin-1') as f:
        for line in f:
            if line.startswith("!Sample_title"):
                metadata["title"] = line.strip().replace('"', '').split('\t')[1:]
            elif line.startswith("!Sample_geo_accession"):
                metadata["gsm"] = line.strip().replace('"', '').split('\t')[1:]
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.strip().replace('"', '').split('\t')[1:]
                
                # Identify Group
                if any("amd_status" in x or "AMD" in x for x in parts) and any("Control" in x for x in parts):
                    if "group" not in metadata:
                        cleaned = []
                        for p in parts:
                            p_lower = p.lower()
                            if "control" in p_lower or "normal" in p_lower: cleaned.append("Control")
                            elif "amd" in p_lower: cleaned.append("AMD")
                            else: cleaned.append("Other")
                        metadata["group"] = cleaned
                
                # Identify Tissue
                if any("tissue:" in x for x in parts):
                    if "tissue" not in metadata:
                        metadata["tissue"] = [p.replace("tissue: ", "").replace('"', '').strip() for p in parts]
                        
    if "gsm" in metadata and "group" in metadata:
        res = {}
        for i, gsm in enumerate(metadata["gsm"]):
            res[gsm] = {
                "group": metadata["group"][i] if i < len(metadata["group"]) else "Unknown",
                "tissue": metadata["tissue"][i] if "tissue" in metadata and i < len(metadata["tissue"]) else "Unknown"
            }
        return res
    return None

def build_counts_matrix(raw_dir, gsm_map):
    print("Building counts matrix from raw files...")
    files = [f for f in os.listdir(raw_dir) if f.endswith(".tsv.gz")]
    counts_dict = {}
    
    processed_count = 0
    for f in files:
        gsm = f.split('_')[0]
        if gsm not in gsm_map or gsm_map[gsm]["group"] not in ["AMD", "Control"]:
            continue
            
        path = os.path.join(raw_dir, f)
        try:
            df = pd.read_csv(path, sep='\t', comment='#', index_col="ID_REF", usecols=["ID_REF", "count"])
            counts_dict[gsm] = df["count"]
            processed_count += 1
            if processed_count % 100 == 0:
                print(f"  Processed {processed_count} files...")
        except Exception as e:
            print(f"Error reading {f}: {e}")
            
    counts_df = pd.DataFrame(counts_dict)
    counts_df.index.name = "ensembl_gene_id"
    return counts_df

def run_de(counts_df, gsm_map):
    print("Running Differential Expression...")
    libsizes = counts_df.sum(axis=0)
    cpm = counts_df.div(libsizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)
    
    amd_cols = [c for c in counts_df.columns if gsm_map[c]["group"] == "AMD"]
    ctrl_cols = [c for c in counts_df.columns if gsm_map[c]["group"] == "Control"]
    
    print(f"  AMD: {len(amd_cols)}, Control: {len(ctrl_cols)}")
    
    logfc = log_cpm[amd_cols].mean(axis=1) - log_cpm[ctrl_cols].mean(axis=1)
    t_stat, p_val = stats.ttest_ind(log_cpm[amd_cols], log_cpm[ctrl_cols], axis=1)
    
    res = pd.DataFrame({
        "logfc": logfc,
        "pvalue": p_val,
        "mean_expr": log_cpm.mean(axis=1)
    })
    return res, log_cpm

def visualize_results(res, log_cpm, gsm_map, ens_to_sym):
    print("Generating visualizations...")
    
    # 1. PCA by Group and Tissue
    print("  Computing PCA...")
    vars = log_cpm.var(axis=1).sort_values(ascending=False)
    top_genes = vars.head(2000).index
    
    pca = PCA(n_components=2)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(log_cpm.loc[top_genes].T)
    pca_coords = pca.fit_transform(scaled_data)
    
    pca_df = pd.DataFrame(pca_coords, columns=["PC1", "PC2"], index=log_cpm.columns)
    pca_df["Group"] = [gsm_map[c]["group"] for c in pca_df.index]
    pca_df["Tissue"] = [gsm_map[c]["tissue"] for c in pca_df.index]
    
    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Tissue", style="Group", alpha=0.7)
    plt.title(f"PCA of GSE135092 (N={len(pca_df)})\nColor by Tissue, Shape by AMD Status")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "gse135092_pca_detailed.pdf"))
    plt.close()
    
    # 2. Volcano Plot
    print("  Plotting Volcano...")
    res["neg_log_p"] = -np.log10(res["pvalue"])
    res["Significant"] = (res["pvalue"] < 0.05) & (res["logfc"].abs() > 0.5)
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=res, x="logfc", y="neg_log_p", hue="Significant", 
                    palette={True: "red", False: "grey"}, alpha=0.5, legend=False, s=10)
    
    risk_genes = ["CFH", "ARMS2", "HTRA1", "C3", "RPE65"]
    for gene_sym in risk_genes:
        ids = [k for k, v in ens_to_sym.items() if v == gene_sym]
        for gene_id in ids:
            if gene_id in res.index:
                row = res.loc[gene_id]
                plt.text(row["logfc"], row["neg_log_p"], gene_sym, fontsize=9, fontweight='bold')
                plt.scatter(row["logfc"], row["neg_log_p"], color='black', s=20)

    plt.title("Volcano Plot: AMD vs Control")
    plt.savefig(os.path.join(OUTPUT_DIR, "gse135092_volcano.pdf"))
    plt.close()
    
    # 3. Violin Plots
    print("  Plotting Risk Gene Densities...")
    risk_ids = {}
    for sym in ["CFH", "RPE65", "HTRA1", "BEST1"]:
        ids = [k for k, v in ens_to_sym.items() if v == sym]
        if ids and ids[0] in log_cpm.index:
            risk_ids[sym] = ids[0]
            
    if risk_ids:
        tidy_data = []
        for sym, gene_id in risk_ids.items():
            expr = log_cpm.loc[gene_id]
            for sample, val in expr.items():
                tidy_data.append({
                    "Gene": sym,
                    "Expression": val,
                    "Group": gsm_map[sample]["group"],
                    "Tissue": gsm_map[sample]["tissue"]
                })
        
        tidy_df = pd.DataFrame(tidy_data)
        plt.figure(figsize=(12, 6))
        sns.violinplot(data=tidy_df, x="Gene", y="Expression", hue="Group", split=True, inner="quart", palette={"AMD": "firebrick", "Control": "skyblue"})
        plt.title("Expression of Key AMD Risk Genes by Disease Status")
        plt.savefig(os.path.join(OUTPUT_DIR, "gse135092_risk_violins.pdf"))
        plt.close()
        
        # New plot: By Tissue
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=tidy_df, x="Gene", y="Expression", hue="Tissue")
        plt.title("Expression of Key AMD Risk Genes by Tissue Type")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, "gse135092_risk_by_tissue.pdf"))
        plt.close()

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    gsm_map = parse_metadata(META_FILE)
    if not gsm_map: return
        
    counts = build_counts_matrix(RAW_DIR, gsm_map)
    if counts.empty: return
        
    res, log_cpm = run_de(counts, gsm_map)
    
    # Map symbols
    map_df = pd.read_csv(os.path.join(BASE_DIR, "results/baseline-expression/all_baseline_expression.csv"))
    ens_to_sym = dict(zip(map_df["ensembl_gene_id"], map_df["gene_name"]))
    res["gene_symbol"] = res.index.map(ens_to_sym).fillna("N/A")
    
    # Save
    res.to_csv(os.path.join(OUTPUT_DIR, "GSE135092_DE_results.csv"))
    
    # Visualize
    visualize_results(res, log_cpm, gsm_map, ens_to_sym)

if __name__ == "__main__":
    main()