
import pandas as pd
df = pd.read_excel("data/RPE_cells/code/RPE_gene pvals.xlsx")
df.to_csv("data/RPE_DE_results.csv", index=False)
print("Converted to CSV.")
