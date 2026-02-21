
import anndata
import pandas as pd

file_path="data/external/perturbseq/K562_gwps_raw_bulk_01.h5ad"
try:
    adata = anndata.read_h5ad(file_path)
    print(adata.obs.columns.tolist())
except FileNotFoundError:
    print(f"Error: File not found at {file_path}")
except ImportError:
    print("Error: Please ensure 'anndata' and 'pandas' libraries are installed.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
