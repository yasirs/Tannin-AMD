
import pandas as pd
try:
    df = pd.read_excel("data/RPE_cells/code/RPE_gene pvals.xlsx", nrows=5)
    print(df.columns)
except Exception as e:
    print(e)
