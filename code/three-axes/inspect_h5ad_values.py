
import scanpy as sc
import numpy as np
try:
    adata = sc.read_h5ad("data/external/perturbseq/rpe1_normalized_bulk_01.h5ad", backed='r')
    # Read a chunk to check values
    X_chunk = adata.X[:100, :100].toarray() if hasattr(adata.X, 'toarray') else adata.X[:100, :100]
    print("Min:", np.min(X_chunk), "Max:", np.max(X_chunk), "Mean:", np.mean(X_chunk))
except Exception as e:
    print(e)
