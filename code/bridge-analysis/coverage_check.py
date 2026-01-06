import os
import pandas as pd
import numpy as np
import scanpy as sc
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats

# ==========================================
# CONFIGURATION
# ==========================================
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
BULK_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/RPE_gene pvals.xlsx")
GMT_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/c2.cp.kegg.v7.0.symbols.gmt")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/coverage-analysis")

PERTURB_FILES = {
    "RPE1_Essential": os.path.join(BASE_DIR, "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad"),
    "K562_Essential": os.path.join(BASE_DIR, "data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad"),
    "K562_GWPS": os.path.join(BASE_DIR, "data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad")
}

CONTRASTS = {
    "H2O2_Stress": {
        "col_pval": "H2O2_vs_CTRL.pvalue",
        "col_padj": "H2O2_vs_CTRL.padj", # Will be created
        "col_symbol": "hgnc_symbol",
        "col_id": "ensembl_gene_id"
    },
    "PRG4_Baseline": {
        "col_pval": "PRG4_vs_CTRL.pvalue",
        "col_padj": "PRG4_vs_CTRL.padj",
        "col_symbol": "hgnc_symbol",
        "col_id": "ensembl_gene_id"
    },
    "PRG4_Rescue": {
        "col_pval": "H2O2PRG4_vs_H2O2.pvalue",
        "col_padj": "H2O2PRG4_vs_H2O2.padj",
        "col_symbol": "hgnc_symbol",
        "col_id": "ensembl_gene_id"
    }
}

THRESHOLD_FDR = 0.05

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

def load_bulk_data():
    print(f"Loading bulk data from {BULK_FILE}...")
    df = pd.read_excel(BULK_FILE)
    
    # Calculate FDR for each contrast
    for name, cols in CONTRASTS.items():
        pval_col = cols["col_pval"]
        mask = df[pval_col].notna()
        pvals = df.loc[mask, pval_col].values
        
        _, adj_pvals, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        
        padj_col = cols["col_padj"]
        df[padj_col] = np.nan
        df.loc[mask, padj_col] = adj_pvals
        
    return df

def get_perturb_info(file_path):
    print(f"Loading metadata from {file_path}...")
    ad = sc.read_h5ad(file_path, backed='r') # Load in backed mode to save memory initially?
    # Actually need to read obs and var.
    # Reading fully might be faster if memory allows. The files are small (<400MB).
    ad = sc.read_h5ad(file_path)
    
    # Genes Measured (var)
    measured_genes_id = set(ad.var_names) # Index is Ensembl ID
    
    # Genes Perturbed (obs)
    # Obs index format: [ID]_[Target_Gene]_[Target_Site]_[Gene_ID] usually, or just check format
    # Let's inspect typical obs name: '10005_ZBTB4_P1_ENSG00000174282'
    # We need the ENSG part at the end for robust matching.
    
    perturbed_ids = set()
    for obs_name in ad.obs_names:
        if "non-targeting" in obs_name: continue
        parts = obs_name.split('_')
        # Assuming format ends with ENSG...
        # Let's verify if the last part looks like an Ensembl ID
        if parts[-1].startswith("ENSG"):
            perturbed_ids.add(parts[-1])
        else:
            # Fallback or different format?
            # K562 GWPS might be different.
            # Let's try to extract gene symbol and map? No, stick to IDs if possible.
            pass
            
    # Also get perturbed symbols for Venn diagrams (easier for display) if needed, 
    # but analysis should be on IDs.
    
    return measured_genes_id, perturbed_ids

def run_enrichment(gene_symbols, pathways, background_symbols):
    results = []
    N = len(background_symbols)
    k = len(gene_symbols)
    
    for name, pathway_genes in pathways.items():
        pathway_genes_bg = pathway_genes.intersection(background_symbols)
        m = len(pathway_genes_bg)
        if m < 5: continue
        
        overlap = gene_symbols.intersection(pathway_genes_bg)
        x = len(overlap)
        if x == 0: continue
        
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
    bulk_df = load_bulk_data()
    pathways = load_gmt(GMT_FILE)
    
    # Pre-load Perturb-seq info
    perturb_data = {}
    for name, path in PERTURB_FILES.items():
        measured, perturbed = get_perturb_info(path)
        perturb_data[name] = {
            "measured": measured,
            "perturbed": perturbed
        }
        print(f"{name}: {len(measured)} measured, {len(perturbed)} perturbed")

    summary_rows = []
    
    all_bulk_symbols = set(bulk_df["hgnc_symbol"].dropna())
    
    for contrast_name, cols in CONTRASTS.items():
        print(f"\nAnalyzing contrast: {contrast_name}")
        
        # Get DEGs
        padj_col = cols["col_padj"]
        id_col = cols["col_id"]
        sym_col = cols["col_symbol"]
        
        # Filter DEGs
        sig_df = bulk_df[bulk_df[padj_col] < THRESHOLD_FDR]
        deg_ids = set(sig_df[id_col].dropna())
        deg_symbols = set(sig_df[sym_col].dropna())
        
        n_degs = len(deg_ids)
        print(f"  N DEGs (FDR < {THRESHOLD_FDR}): {n_degs}")
        
        row = {"Contrast": contrast_name, "N_DEGs": n_degs}
        
        # For Venn Diagram export
        venn_export = {
            "DEG_IDs": list(deg_ids)
        }
        
        # Check against each Perturb-seq dataset
        for p_name, p_info in perturb_data.items():
            measured = p_info["measured"]
            perturbed = p_info["perturbed"]
            
            n_meas = len(deg_ids.intersection(measured))
            n_pert = len(deg_ids.intersection(perturbed))
            
            row[f"Measured_{p_name}"] = n_meas
            row[f"KD_{p_name}"] = n_pert
            row[f"Pct_KD_{p_name}"] = (n_pert / n_degs * 100) if n_degs > 0 else 0
            
            venn_export[f"KD_{p_name}_IDs"] = list(perturbed)
            
            # Missing Genes Analysis
            # 1. Missing from RPE1 (Essential)
            if p_name == "RPE1_Essential":
                missing_ids = deg_ids - perturbed
                missing_symbols = set(bulk_df[bulk_df[id_col].isin(missing_ids)][sym_col].dropna())
                
                # Save missing list
                missing_df = pd.DataFrame({"ensembl_gene_id": list(missing_ids)})
                missing_df = missing_df.merge(bulk_df[[id_col, sym_col]], on=id_col, how="left").drop_duplicates()
                missing_df.to_csv(os.path.join(OUTPUT_DIR, f"missing_genes_{contrast_name}_{p_name}.csv"), index=False)
                
                # Enrichment on Missing
                if missing_symbols:
                    enr_res = run_enrichment(missing_symbols, pathways, all_bulk_symbols)
                    if not enr_res.empty:
                        enr_res["contrast"] = contrast_name
                        enr_res["missing_from"] = p_name
                        enr_res.head(20).to_csv(os.path.join(OUTPUT_DIR, f"enrichment_missing_{contrast_name}_{p_name}.csv"), index=False)

            # 2. Missing from K562 GWPS
            if p_name == "K562_GWPS":
                missing_ids = deg_ids - perturbed
                missing_symbols = set(bulk_df[bulk_df[id_col].isin(missing_ids)][sym_col].dropna())
                
                # Save missing list
                missing_df = pd.DataFrame({"ensembl_gene_id": list(missing_ids)})
                missing_df = missing_df.merge(bulk_df[[id_col, sym_col]], on=id_col, how="left").drop_duplicates()
                missing_df.to_csv(os.path.join(OUTPUT_DIR, f"missing_genes_{contrast_name}_{p_name}.csv"), index=False)
                
                # Enrichment
                if missing_symbols:
                    enr_res = run_enrichment(missing_symbols, pathways, all_bulk_symbols)
                    if not enr_res.empty:
                        enr_res["contrast"] = contrast_name
                        enr_res["missing_from"] = p_name
                        enr_res.head(20).to_csv(os.path.join(OUTPUT_DIR, f"enrichment_missing_{contrast_name}_{p_name}.csv"), index=False)

        summary_rows.append(row)

    # Save Summary
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(OUTPUT_DIR, "coverage_summary.csv"), index=False)
    print("\nCoverage Summary:")
    print(summary_df.to_string())

if __name__ == "__main__":
    main()
