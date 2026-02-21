
import os
import pandas as pd
import gseapy as gp
import requests

# Constants
OUT_DIR = "results/three-axes/gene-sets"
os.makedirs(OUT_DIR, exist_ok=True)
GMT_FILE = "data/c2.all.v2023.2.Hs.symbols.gmt"

# Manual Lists
KNOWN_AMD_GENES = [
    "CFH", "C3", "ARMS2", "HTRA1", "APOE", "TIMP3", "RPE65", "BEST1", "RLBP1", "TTR",
    "CFI", "C9", "LIPC", "CETP", "COL8A1", "COL10A1", "TNFRSF10A", "IER3", "SLC16A8", "B3GALTL"
]
# Add more from literature if needed
KNOWN_AMD_GENES = sorted(list(set(KNOWN_AMD_GENES)))

def parse_gmt(gmt_path, gene_set_name):
    """Parse a local GMT file for a specific gene set."""
    if not os.path.exists(gmt_path):
        print(f"GMT file not found: {gmt_path}")
        return []
    
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == gene_set_name:
                return parts[2:] # Genes start at index 2
    return []

def fetch_enrichr(library, term):
    """Fetch gene set from Enrichr."""
    try:
        # Using gseapy's get_library_name to verify library exists is optional
        # but fetch_results might fail if term not found.
        # We use gseapy.get_library_name() manually checked before.
        # Actually gseapy.get_library_name() returns libraries, not terms.
        # We need to use gseapy.get_library_name() ? No.
        # We'll valid libraries: MSigDB_Hallmark_2020, GO_Biological_Process_2023
        
        # We can't easily get a single term via gseapy without running enrichment.
        # But we can get the full library as a dict.
        print(f"Fetching library: {library}")
        gs_dict = gp.get_library(name=library, organism='Human')
        
        # Find key matching term (case insensitive or exact)
        # Enrichr terms sometimes differ slightly.
        matched_key = None
        for k in gs_dict.keys():
            if term.lower() in k.lower(): # Loose match
                matched_key = k
                break
        
        if matched_key:
            print(f"  Found term: {matched_key}")
            return gs_dict[matched_key]
        else:
            print(f"  Term '{term}' not found in {library}")
            return []
            
    except Exception as e:
        print(f"Error fetching {library}: {e}")
        return []

def save_genes(gene_list, filename):
    if not gene_list:
        print(f"Warning: Empty gene list for {filename}")
        return
    df = pd.DataFrame({'gene_symbol': gene_list})
    df = df.drop_duplicates().sort_values('gene_symbol')
    path = os.path.join(OUT_DIR, filename)
    df.to_csv(path, index=False)
    print(f"Saved {len(df)} genes to {path}")

def main():
    print("Starting Gene Set Curation...")
    
    # 1. Oxidative Axis
    print("\n--- Oxidative Axis ---")
    ox_hallmark = fetch_enrichr("MSigDB_Hallmark_2020", "Oxidative Phosphorylation")
    ros_hallmark = fetch_enrichr("MSigDB_Hallmark_2020", "Reactive Oxygen Species Pathway") # Added ROS pathway
    ox_go = fetch_enrichr("GO_Biological_Process_2023", "response to oxidative stress (GO:0006979)")
    
    # NRF2 Targets - try "NFE2L2" in ChEA or similar
    # Using 'ChEA_2022' or 'TRRUST_Transcription_Factors_2019'
    nrf2_targets = fetch_enrichr("TRRUST_Transcription_Factors_2019", "NFE2L2")
    if not nrf2_targets:
        nrf2_targets = fetch_enrichr("ChEA_2022", "NFE2L2")
        
    pro_oxidative = sorted(list(set(ox_hallmark + ros_hallmark + ox_go))) # Just a mix for now, will refine
    save_genes(pro_oxidative, "oxidative_pro_disease.csv")
    
    anti_oxidative = nrf2_targets # NRF2 is antioxidant response
    save_genes(anti_oxidative, "oxidative_anti_disease.csv")
    
    
    # 2. Inflammatory Axis
    print("\n--- Inflammatory Axis ---")
    inf_hallmark = fetch_enrichr("MSigDB_Hallmark_2020", "Inflammatory Response")
    comp_hallmark = fetch_enrichr("MSigDB_Hallmark_2020", "Complement")
    nfkb_targets = fetch_enrichr("MSigDB_Hallmark_2020", "TNF-alpha Signaling via NF-kB") # Corrected name
    il6_stat3 = fetch_enrichr("MSigDB_Hallmark_2020", "IL-6/JAK/STAT3 Signaling") # Added related
    
    pro_inflammatory = sorted(list(set(inf_hallmark + comp_hallmark + nfkb_targets + il6_stat3)))
    save_genes(pro_inflammatory, "inflammatory_pro_disease.csv")
    
    # Anti-inflammatory? Literature curated often hard. PRG4 is anti-inflammatory.
    # We will leave anti-inflammatory empty or user-defined for now unless we find "Anti-inflammatory response" GO
    # Or "Resolution of inflammation"
    res_inf = fetch_enrichr("GO_Biological_Process_2023", "regulation of inflammatory response") # Broad
    # Let's save a placeholder
    save_genes(res_inf, "inflammatory_anti_disease.csv") # Saved broad regulation terms
    
    
    # 3. Senescence Axis
    print("\n--- Senescence Axis ---")
    sen_mayo = parse_gmt(GMT_FILE, "SAUL_SEN_MAYO")
    print(f"  Parsed SAUL_SEN_MAYO: {len(sen_mayo)} genes")
    
    sasp_reactome = parse_gmt(GMT_FILE, "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP")
    print(f"  Parsed REACTOME_SASP: {len(sasp_reactome)} genes")
    
    # GO Senescence
    sen_go = fetch_enrichr("GO_Biological_Process_2023", "cellular senescence")
    
    # Gold standards
    gold_sen = ["CDKN1A", "CDKN2A", "TP53", "RB1", "LMNB1", "GLB1"] # p21, p16, p53, Rb, Lamin B1, Beta-gal
    
    pro_senescence = sorted(list(set(sen_mayo + sasp_reactome + sen_go + gold_sen)))
    save_genes(pro_senescence, "senescence_pro_disease.csv")
    
    # Anti-senescence? Pro-proliferation?
    # Maybe "Cell Cycle" hallmark
    cell_cycle = fetch_enrichr("MSigDB_Hallmark_2020", "G2-M Checkpoint") # Corrected Name
    e2f_targets = fetch_enrichr("MSigDB_Hallmark_2020", "E2F Targets") # Added E2F targets (pro-proliferation)
    
    pro_prolif = sorted(list(set(cell_cycle + e2f_targets)))
    save_genes(pro_prolif, "senescence_anti_disease.csv")
    
    
    # 4. Known AMD
    print("\n--- Known AMD Genes ---")
    save_genes(KNOWN_AMD_GENES, "known_amd_genes.csv")

if __name__ == "__main__":
    main()
