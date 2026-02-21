
import gseapy as gp
lib = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')
print(sorted(lib.keys()))
