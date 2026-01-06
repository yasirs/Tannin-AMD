import os
import pandas as pd

# Paths
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
GWAS_FILE = os.path.join(BASE_DIR, "data/external/gwas/amd_gwas_loci.csv")
BULK_FILE = os.path.join(BASE_DIR, "data/RPE_cells/code/PRG4_Rescue_DEGs.csv") # Using the one with hgnc_symbol
SCREEN_FILE = os.path.join(BASE_DIR, "results/virtual-screen/prg4_virtual_screen_results.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/gwas-integration")

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # 1. Load Data
    gwas_df = pd.read_csv(GWAS_FILE)
    bulk_df = pd.read_csv(BULK_FILE)
    screen_df = pd.read_csv(SCREEN_FILE)
    
    # 2. Rename columns for merging
    # Screen columns: ['perturbation', 'spearman_rho', 'gene_symbol', 'ensembl_id']
    if 'gene_symbol' in screen_df.columns:
        screen_df = screen_df.rename(columns={'gene_symbol': 'hgnc_symbol'})
    if 'spearman_rho' in screen_df.columns:
        screen_df = screen_df.rename(columns={'spearman_rho': 'correlation'})

    # 3. Merge Bulk with GWAS
    bulk_gwas = pd.merge(bulk_df, gwas_df, on='hgnc_symbol')
    
    # 4. Merge Screen with GWAS
    screen_gwas = pd.merge(screen_df, gwas_df, on='hgnc_symbol')

    # 5. Analysis: Are GWAS genes rescued by PRG4?
    # PRG4 Rescue = H2O2PRG4_vs_H2O2.log2FoldChange
    rescued_gwas = bulk_gwas[bulk_gwas['H2O2PRG4_vs_H2O2.log2FoldChange'] > 0.3].copy()
    rescued_gwas = rescued_gwas.sort_values('H2O2PRG4_vs_H2O2.log2FoldChange', ascending=False)
    
    print(f"GWAS genes rescued by PRG4: {len(rescued_gwas)}")
    print(rescued_gwas[['hgnc_symbol', 'H2O2PRG4_vs_H2O2.log2FoldChange', 'Locus']])

    # 6. Analysis: Are GWAS genes targets in the virtual screen?
    # Top mimetics (high correlation)
    top_mimetics_gwas = screen_gwas[screen_gwas['correlation'] > 0.05].sort_values('correlation', ascending=False)
    
    # Also find top antagonists (low correlation)
    top_antagonists_gwas = screen_gwas[screen_gwas['correlation'] < -0.05].sort_values('correlation', ascending=True)

    # 7. Save results
    rescued_gwas.to_csv(os.path.join(OUTPUT_DIR, "prg4_rescued_gwas_genes.csv"), index=False)
    top_mimetics_gwas.to_csv(os.path.join(OUTPUT_DIR, "gwas_genes_in_virtual_screen.csv"), index=False)

    # 8. Summary Report
    with open(os.path.join(OUTPUT_DIR, "gwas_integration_summary.md"), "w") as f:
        f.write("# AMD GWAS Integration Summary\n\n")
        f.write(f"Overlap between IAMDGC loci and PRG4 rescue: {len(rescued_gwas)} genes.\n\n")
        f.write("### PRG4-Rescued GWAS Genes\n")
        f.write("These genes are linked to AMD risk and are upregulated by PRG4 treatment in H2O2-stressed RPE.\n\n")
        f.write(rescued_gwas[['hgnc_symbol', 'H2O2PRG4_vs_H2O2.log2FoldChange', 'Evidence']].to_markdown(index=False))
        f.write("\n\n### GWAS Genes as Virtual Screen Hits (Mimetics)\n")
        f.write("Knockdown of these genes mimics the PRG4 rescue signature. High correlation means knockdown ~ PRG4.\n\n")
        if not top_mimetics_gwas.empty:
            f.write(top_mimetics_gwas[['hgnc_symbol', 'correlation', 'Evidence']].head(10).to_markdown(index=False))
        else:
            f.write("No GWAS genes found as mimetics above threshold.")
        
        f.write("\n\n### GWAS Genes as Virtual Screen Hits (Antagonists)\n")
        f.write("Knockdown of these genes opposes the PRG4 rescue signature.\n\n")
        if not top_antagonists_gwas.empty:
            f.write(top_antagonists_gwas[['hgnc_symbol', 'correlation', 'Evidence']].head(10).to_markdown(index=False))
        else:
            f.write("No GWAS genes found as antagonists above threshold.")


if __name__ == "__main__":
    main()
