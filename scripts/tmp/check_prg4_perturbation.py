
import anndata
import pandas as pd
import sys

file_path="data/external/perturbseq/K562_gwps_raw_bulk_01.h5ad"
try:
    adata = anndata.read_h5ad(file_path)
    obs_df = adata.obs

    if "gene_transcript" in obs_df.columns:
        found_prg4 = False
        for transcript_id in obs_df["gene_transcript"]:
            parts = transcript_id.split("_")
            if len(parts) > 1 and parts[1] == "PRG4":
                found_prg4 = True
                break
        
        if found_prg4:
            print("Yes, PRG4 is listed as a perturbation target gene in the K562 GWPS dataset.")
        else:
            print("No, PRG4 is not explicitly listed as a perturbation target gene in the K562 GWPS dataset based on the 'gene_transcript' format.")
            
    else:
        print("Error: 'gene_transcript' column not found in the observations.")

except FileNotFoundError:
    print(f"Error: File not found at {file_path}")
except ImportError:
    print("Error: Please ensure 'anndata' and 'pandas' libraries are installed within the activated virtual environment.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
