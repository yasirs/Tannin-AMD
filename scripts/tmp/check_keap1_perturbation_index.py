
import anndata
import pandas as pd

file_path="data/external/perturbseq/K562_gwps_raw_bulk_01.h5ad"
try:
    adata = anndata.read_h5ad(file_path)
    obs_index = adata.obs.index

    found_keap1 = False
    for item in obs_index:
        # The README states the format is typically [ID]_[Target_Gene]_[Target_Site]_[Gene_ID]
        parts = str(item).split("_")
        if len(parts) > 1 and parts[1] == "KEAP1":
            found_keap1 = True
            break
    
    if found_keap1:
        print("Yes, KEAP1 is listed as a perturbation target gene in the K562 GWPS dataset.")
    else:
        print("No, KEAP1 is not explicitly listed as a perturbation target gene in the K562 GWPS dataset based on the index format.")
        
except FileNotFoundError:
    print(f"Error: File not found at {file_path}")
except ImportError:
    print("Error: Please ensure 'anndata' and 'pandas' libraries are installed.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
