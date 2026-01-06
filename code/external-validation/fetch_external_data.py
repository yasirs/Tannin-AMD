import os
import pandas as pd
import numpy as np
import requests
import gzip
import io
import sys

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
OUTPUT_DIR = os.path.join(BASE_DIR, "data/external/geo")

DATASETS = [
    {"id": "GSE135092", "url_group": "GSE135nnn"},
    {"id": "GSE29801",  "url_group": "GSE29nnn"}
]

def download_file(url, local_path):
    print(f"Downloading {url}...")
    try:
        r = requests.get(url, stream=True)
        r.raise_for_status()
        with open(local_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete.")
        return True
    except Exception as e:
        print(f"Download failed: {e}")
        return False

def parse_geo_matrix(local_path, gse_id):
    print(f"Parsing Series Matrix for {gse_id}...")
    
    metadata = {}
    row_skip = 0
    
    # 1. Parse Metadata first
    with gzip.open(local_path, 'rt', encoding='latin-1') as f:
        for i, line in enumerate(f):
            if line.startswith("!Sample_title"):
                parts = line.strip().replace('"', '').split('\t') 
                metadata["title"] = parts[1:]
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.strip().replace('"', '').split('\t')
                if "characteristics" not in metadata:
                    metadata["characteristics"] = []
                metadata["characteristics"].append(parts[1:])
            elif line.startswith("!series_matrix_table_begin"):
                row_skip = i + 1
                break
    
    # 2. Parse Data Table
    try:
        df = pd.read_csv(local_path, compression='gzip', sep='\t', skiprows=row_skip, index_col=0, comment='!', low_memory=False)
    except Exception as e:
        print(f"Error reading dataframe: {e}")
        return None

    # Clean columns
    df.columns = [str(c).replace('"', '') for c in df.columns]
    
    # 3. Determine Groups
    n_samples = len(metadata["title"])
    groups = ["Unknown"] * n_samples
    
    # Heuristic for group assignment
    if "characteristics" in metadata:
        for row in metadata["characteristics"]:
            # Check if this row is relevant (contains AMD/Control info)
            row_str = " ".join(row).lower()
            if "amd" in row_str or "control" in row_str or "normal" in row_str:
                # This row likely contains the status. Process it.
                # We need to process this row for ALL samples.
                
                # Check strict values if possible, otherwise heuristic
                for i, val in enumerate(row):
                    val_lower = str(val).lower()
                    
                    # Prioritize Control/Normal to avoid "amd_status: control" matching "amd"
                    if "normal" in val_lower or "control" in val_lower:
                        groups[i] = "Control"
                    elif "amd" in val_lower or "age-related macular degeneration" in val_lower:
                         # Ensure it's not just the key "amd status:"
                         # But since we checked control first, "amd status: amd" should fall through here.
                         # "amd status: control" is caught above.
                        groups[i] = "AMD"
                
                # If we found meaningful assignments, stop looking at other rows
                # (unless we want to combine info, but usually one row is sufficient)
                if groups.count("AMD") > 0 and groups.count("Control") > 0:
                    break


    valid_indices = [i for i, g in enumerate(groups) if g in ["AMD", "Control"]]
    valid_groups = [groups[i] for i in valid_indices]
    
    print(f"  Found {groups.count('AMD')} AMD and {groups.count('Control')} Control samples.")
    
    if len(valid_indices) < 2:
        print("  Not enough samples found.")
        return None
        
    df_sub = df.iloc[:, valid_indices]
    
    # Numeric conversion
    df_sub = df_sub.apply(pd.to_numeric, errors='coerce')
    
    # Log transform check
    if df_sub.max().max() > 50:
        print("  Max value > 50, applying log2...")
        df_sub = np.log2(df_sub + 1)
        
    # Calculate LogFC
    amd_cols = [df_sub.columns[i] for i, g in enumerate(valid_groups) if g == "AMD"]
    ctrl_cols = [df_sub.columns[i] for i, g in enumerate(valid_groups) if g == "Control"]
    
    amd_mean = df_sub[amd_cols].mean(axis=1)
    ctrl_mean = df_sub[ctrl_cols].mean(axis=1)
    logfc = amd_mean - ctrl_mean
    
    res_df = pd.DataFrame({
        "gene_symbol": df_sub.index, 
        "logfc_amd_vs_ctrl": logfc,
        "pvalue": [np.nan] * len(logfc) # Placeholder
    })
    
    # T-test for p-values (simple)
    from scipy import stats
    t_stat, p_val = stats.ttest_ind(df_sub[amd_cols], df_sub[ctrl_cols], axis=1, nan_policy='omit')
    res_df["pvalue"] = p_val
    
    return res_df, df_sub, valid_groups

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    for ds in DATASETS:
        gse_id = ds['id']
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{ds['url_group']}/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
        local_path = os.path.join(OUTPUT_DIR, f"{gse_id}_series_matrix.txt.gz")
        
        if not os.path.exists(local_path):
            success = download_file(url, local_path)
            if not success: continue
            
        result = parse_geo_matrix(local_path, gse_id)
        if result is not None:
            res_df, expr_df, groups = result
            
            # Save Signature
            out_path = os.path.join(OUTPUT_DIR, f"{gse_id}_amd_signature.csv")
            res_df.to_csv(out_path)
            
            # Save Expression Matrix (Transposed: Samples x Genes) for easier loading
            # Add Group info
            expr_df = expr_df.T
            expr_df["Group"] = groups
            expr_path = os.path.join(OUTPUT_DIR, f"{gse_id}_expression_matrix.csv.gz")
            expr_df.to_csv(expr_path, compression='gzip')
            
            print(f"Saved {gse_id} signature to {out_path}")
            print(f"Saved {gse_id} expression matrix to {expr_path}")

if __name__ == "__main__":
    main()
