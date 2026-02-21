
import pandas as pd

# Load DE results (missing symbols)
de_file = "data/RPE_DE_results.csv"
de_df = pd.read_csv(de_file)

# Load TPMS (hopefully has symbols)
tpm_file = "data/RPE_TPMS.csv"
tpm_df = pd.read_csv(tpm_file)

# Check if tpm_df has valid symbols
print("TPM Symbols valid:", tpm_df['hgnc_symbol'].notna().sum(), "/", len(tpm_df))

if tpm_df['hgnc_symbol'].notna().sum() > 0:
    # Create mapping
    mapping = tpm_df[['ensembl_gene_id', 'hgnc_symbol']].dropna().drop_duplicates('ensembl_gene_id')
    
    # Merge
    # Drop old empty symbol column if exists
    if 'hgnc_symbol' in de_df.columns:
        de_df = de_df.drop(columns=['hgnc_symbol'])
        
    merged = pd.merge(de_df, mapping, on='ensembl_gene_id', how='left')
    
    # Save
    merged.to_csv("data/RPE_DE_results.csv", index=False)
    print("Updated RPE_DE_results.csv with symbols from TPMS.")
    print(merged.head())
else:
    print("TPMS also missing symbols!")
