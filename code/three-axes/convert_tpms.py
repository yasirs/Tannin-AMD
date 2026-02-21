
import pandas as pd
try:
    df = pd.read_excel("data/RPE_cells/code/RPE_TPMS.xlsx")
    df.to_csv("data/RPE_TPMS.csv", index=False)
    print("Converted TPMs to CSV.")
except Exception as e:
    print(e)
