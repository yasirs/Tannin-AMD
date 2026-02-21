import gseapy
import pandas as pd

try:
    names = gseapy.get_library_name()
    print("Available libraries sample:", names[:20])
    print("MSigDB related:", [n for n in names if 'MSigDB' in n])
    print("GO related:", [n for n in names if 'GO_Biological' in n])
    print("Pathway related:", [n for n in names if 'Wiki' in n or 'KEGG' in n or 'Reactome' in n])
    print("TF related:", [n for n in names if 'TRRUST' in n or 'ChEA' in n])
except Exception as e:
    print(e)
