# GWAS Integration Results

This directory contains the integration of AMD Genome-Wide Association Study (GWAS) loci with the PRG4 rescue analysis.

## Key Findings
1.  **Direct Rescue of Risk Genes:** PRG4 treatment significantly upregulates several key AMD risk genes that are downregulated by H2O2 stress, including:
    *   **CFH** (Complement Factor H): The primary genetic risk factor for AMD.
    *   **C3** (Complement Component 3).
    *   **TRPM1**: Linked to retinal signaling.
2.  **Mimetics in Risk Loci:** Knockdown of **CFH** and **B3GALNT2** show positive correlation with the PRG4 rescue signature in the genome-wide screen.

## Files
*   **`gwas_integration_summary.md`**: Summary of overlapping genes and their statistics.
*   **`prg4_rescued_gwas_genes.csv`**: Detailed list of GWAS genes upregulated by PRG4.
*   **`gwas_genes_in_virtual_screen.csv`**: Correlation of GWAS gene knockdowns with PRG4 rescue.

## Analysis Script
`code/bridge-analysis/validate_gwas_integration.py`
