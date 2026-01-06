import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Paths
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
ATAC_FILE = os.path.join(BASE_DIR, "data/external/geo/GSE99287/GSE99287_RPE_ATACSeq_peak_counts.txt.gz")
BULK_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/H2O2_DEGs.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/atac-integration")

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # 1. Load ATAC Data
    print("Loading ATAC-seq data...")
    atac_df = pd.read_csv(ATAC_FILE, sep='\t')
    
    # Pre-process genes: split comma-separated genes and explode
    atac_df['Gene'] = atac_df['Gene'].fillna('')
    atac_df['Gene_list'] = atac_df['Gene'].apply(lambda x: [g.strip() for g in x.split(',') if g.strip()])
    atac_exploded = atac_df.explode('Gene_list')
    atac_exploded = atac_exploded.rename(columns={'Gene_list': 'hgnc_symbol'})

    # 2. Load Bulk RNA-seq Data
    print("Loading Bulk RNA-seq data...")
    bulk_df = pd.read_csv(BULK_FILE)
    # Drop rows without hgnc_symbol
    bulk_df = bulk_df.dropna(subset=['hgnc_symbol'])

    # 3. Merge
    print("Merging datasets...")
    merged = pd.merge(atac_exploded, bulk_df, on='hgnc_symbol')
    
    # 4. Filter for Differentially Accessible Regions (DARs)
    # The file has a 'Differentially_accessible_region' column
    # Note: All DARs in this file are LOSSES in AMD (Coefficient < 0)
    
    # Aggregate ATAC data per gene to avoid duplicates
    # Use min coefficient to represent the strongest loss for a gene
    atac_gene = atac_exploded.groupby('hgnc_symbol').agg({
        'Coefficient_of_stage': 'min',
        'Differentially_accessible_region': 'max',
        'Fdr_of_stage_coefficient': 'min'
    }).reset_index()

    # Merge with Bulk
    merged = pd.merge(atac_gene, bulk_df, on='hgnc_symbol')
    
    dars = merged[merged['Differentially_accessible_region'] == 1].copy()
    non_dars = merged[merged['Differentially_accessible_region'] == 0].copy()
    
    print(f"Total overlapping genes: {len(merged)}")
    print(f"Total overlapping DAR-linked genes: {len(dars)}")
    print(f"Total overlapping non-DAR-linked genes: {len(non_dars)}")

    # 5. Analysis: ATAC AMD Change vs H2O2 Response (Scatter)
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=merged,
        x='Coefficient_of_stage',
        y='H2O2_vs_CTRL.log2FoldChange',
        hue='Differentially_accessible_region',
        palette={0: 'lightgrey', 1: 'red'},
        alpha=0.5
    )
    plt.axhline(0, color='black', linestyle='--')
    plt.axvline(0, color='black', linestyle='--')
    plt.title("Chromatin Accessibility (AMD) vs. Expression (H2O2)")
    plt.xlabel("ATAC-seq AMD Stage Coefficient (Negative = Loss in AMD)")
    plt.ylabel("Bulk RNA-seq H2O2 vs CTRL Log2FC")
    plt.savefig(os.path.join(OUTPUT_DIR, "atac_vs_h2o2_scatter.pdf"))
    plt.close()

    # 6. Analysis: PRG4 Rescue Effect on DAR vs Non-DAR
    plt.figure(figsize=(8, 6))
    data_to_plot = [
        dars['H2O2PRG4_vs_H2O2.log2FoldChange'].dropna(),
        non_dars['H2O2PRG4_vs_H2O2.log2FoldChange'].dropna()
    ]
    plt.boxplot(data_to_plot, tick_labels=['AMD DAR Genes', 'Other Genes'])
    plt.ylabel("PRG4 Rescue Effect (Log2FC H2O2+PRG4 vs H2O2)")
    plt.title("PRG4 Rescue Effect on AMD DAR-linked Genes")
    
    # Statistical test
    t_stat, p_val = stats.ttest_ind(data_to_plot[0], data_to_plot[1])
    plt.text(1.5, plt.ylim()[1]*0.9, f"t-test p = {p_val:.2e}", ha='center')
    
    plt.savefig(os.path.join(OUTPUT_DIR, "prg4_rescue_on_dars_boxplot.pdf"))
    plt.close()

    # 7. Identify top candidate targets
    # Genes that lose accessibility in AMD, are down in H2O2, and rescued by PRG4
    candidates = dars[
        (dars['H2O2_vs_CTRL.log2FoldChange'] < -0.5) &
        (dars['H2O2PRG4_vs_H2O2.log2FoldChange'] > 0.5)
    ]
    
    candidates = candidates.sort_values('H2O2PRG4_vs_H2O2.log2FoldChange', ascending=False)
    candidates[['hgnc_symbol', 'Coefficient_of_stage', 'H2O2_vs_CTRL.log2FoldChange', 'H2O2PRG4_vs_H2O2.log2FoldChange']].to_csv(
        os.path.join(OUTPUT_DIR, "top_prg4_rescue_dar_candidates.csv"), index=False
    )
    
    print(f"Found {len(candidates)} candidate genes rescued by PRG4 in AMD DAR regions.")

    # 8. Functional Enrichment of Candidates
    # Check if candidates are enriched for certain pathways (optional but good)
    # For now, let's just print top genes
    print("Top candidates:")
    print(candidates[['hgnc_symbol', 'H2O2PRG4_vs_H2O2.log2FoldChange']].head(10))

    # Save summary report
    with open(os.path.join(OUTPUT_DIR, "atac_integration_summary.md"), "w") as f:
        f.write("# ATAC-seq Integration Analysis Summary\n\n")
        f.write(f"- Total overlapping genes between GSE99287 and internal bulk: {len(merged)}\n")
        f.write(f"- Genes linked to AMD Differentially Accessible Regions (DARs): {len(dars)}\n")
        f.write(f"- Note: All DARs in this dataset represent **decreased** accessibility in AMD.\n")
        f.write(f"- Number of candidates (DAR-linked, Down in H2O2, Rescued by PRG4): {len(candidates)}\n")
        f.write(f"\n### Statistical Comparison\n")
        f.write(f"Does PRG4 rescue DAR-linked genes more than others?\n")
        f.write(f"- Mean Rescue (DAR): {data_to_plot[0].mean():.4f}\n")
        f.write(f"- Mean Rescue (Non-DAR): {data_to_plot[1].mean():.4f}\n")
        f.write(f"- t-test p-value: {p_val:.2e}\n")
        f.write(f"\nTop Candidates:\n\n")
        f.write(candidates[['hgnc_symbol', 'H2O2PRG4_vs_H2O2.log2FoldChange']].head(20).to_markdown(index=False))


if __name__ == "__main__":
    main()
