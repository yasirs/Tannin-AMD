
import scanpy as sc
try:
    adata = sc.read_h5ad("data/external/perturbseq/rpe1_normalized_bulk_01.h5ad", backed='r')
    print("Obs columns:", adata.obs.columns.tolist())
    print("Obs index sample:", adata.obs_names[:3].tolist())
except Exception as e:
    print(e)
