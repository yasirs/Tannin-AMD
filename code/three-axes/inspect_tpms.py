
import pandas as pd
df = pd.read_csv("data/RPE_TPMS.csv", nrows=5)
print(df.columns.tolist())
