
import pandas as pd
df = pd.read_excel("data/RPE_cells/code/RPE_gene pvals.xlsx", nrows=10)
print(df[['ensembl_gene_id', 'hgnc_symbol']])
