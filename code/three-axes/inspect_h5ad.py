
import scanpy as sc
try:
    adata = sc.read_h5ad("data/external/perturbseq/rpe1_normalized_bulk_01.h5ad", backed='r')
    print("Var names sample:", adata.var_names[:5].tolist())
    print("Obs names sample:", adata.obs_names[:5].tolist())
except Exception as e:
    print(e)
