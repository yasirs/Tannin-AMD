"""
Configuration for Perturb-seq Gene Set Analysis

Central configuration file containing paths, parameters, and settings.
"""

import os
from pathlib import Path

#%%
# Base paths
PROJECT_ROOT = Path("/home/ysuhail/work/Tannin-AMD")
CODE_DIR = PROJECT_ROOT / "code" / "perturbseq-analysis"
RESULTS_DIR = PROJECT_ROOT / "results" / "perturbseq-analysis"

#%%
# Input data paths
PERTURBSEQ_DATA = {
    "RPE1": PROJECT_ROOT / "data" / "external" / "perturbseq" / "rpe1_normalized_bulk_01.h5ad",
    "K562_GWPS": PROJECT_ROOT / "data" / "external" / "perturbseq" / "K562_gwps_normalized_bulk_01.h5ad",
}

#%%
# Gene set source files
GENE_SET_SOURCES = {
    # Age gene sets
    "age_de": PROJECT_ROOT / "results" / "cohort-GSE29801" / "age_nrf2_analysis" / "age_de_all.csv",
    "are": PROJECT_ROOT / "results" / "cohort-GSE29801" / "age_nrf2_analysis" / "ARE_gene_set.csv",
    
    # AMD gene sets (with gene symbols)
    "amd_macula": PROJECT_ROOT / "results" / "cohort-GSE135092" / "rpe_covariate_de" / "RPE_Macula_all_results_symbols.csv",
    "amd_nonmacula": PROJECT_ROOT / "results" / "cohort-GSE135092" / "rpe_covariate_de" / "RPE_nonMacula_all_results_symbols.csv",
    
    # Axis gene sets
    "senescence_pro": PROJECT_ROOT / "results" / "three-axes" / "gene-sets" / "senescence_pro_disease.csv",
    "senescence_anti": PROJECT_ROOT / "results" / "three-axes" / "gene-sets" / "senescence_anti_disease.csv",
    "redox_pro": PROJECT_ROOT / "results" / "three-axes" / "gene-sets" / "oxidative_pro_disease.csv",
    "redox_anti": PROJECT_ROOT / "results" / "three-axes" / "gene-sets" / "oxidative_anti_disease.csv",
    "inflammation_pro": PROJECT_ROOT / "results" / "three-axes" / "gene-sets" / "inflammatory_pro_disease.csv",
    "inflammation_anti": PROJECT_ROOT / "results" / "three-axes" / "gene-sets" / "inflammatory_anti_disease.csv",
}

#%%
# Output directories
OUTPUT_DIRS = {
    "gene_sets": RESULTS_DIR / "gene-sets",
    "cache": RESULTS_DIR / "cache",
    "Age": RESULTS_DIR / "Age",
    "AMD": RESULTS_DIR / "AMD", 
    "Axes": RESULTS_DIR / "Axes",
    "concordance": RESULTS_DIR / "concordance",
}

#%%
# Analysis parameters
ANALYSIS_PARAMS = {
    # DE thresholds for AMD gene extraction
    "amd_pvalue_threshold": 0.01,
    
    # Mimetic/antagonist selection
    "top_n_mimetics": 200,
    "top_n_antagonists": 200,
    
    # GSEA parameters
    "min_gene_set_size": 10,
    "max_gene_set_size": 2000,
    "permutations": 1000,
}

#%%
# Visualization settings (following coding_preferences.md)
VIZ_SETTINGS = {
    # Font
    "font_family": "serif",  # Will use Palatino-like font
    "font_size": 10,
    
    # Figure size (small physical dimensions)
    "fig_width": 5,
    "fig_height": 4,
    "fig_dpi": 300,
    
    # Colors
    "color_significant": "#1f77b4",  # Blue for significant
    "color_nonsignificant": "#CCCCCC",  # Gray for non-significant
    "cmap_diverging": "BrBG",  # Brown-Blue-Green (brown for high)
    "cmap_sequential": "YlOrBr",  # Yellow-Orange-Brown
    
    # Category colors
    "category_colors": {
        "Age": "#E41A1C",
        "AMD": "#377EB8", 
        "Senescence": "#4DAF4A",
        "Redox": "#984EA3",
        "Inflammation": "#FF7F00",
    },
    
    # Output formats
    "save_pdf": True,
    "save_tiff": True,
    "tiff_compression": "lzw",
}

#%%
# Gene set metadata (will be populated by gene_set_registry.py)
GENE_SET_CATEGORIES = {
    "Age-UP": {"category": "Age", "direction": "UP", "description": "Upregulated with age (GSE29801)"},
    "Age-DOWN": {"category": "Age", "direction": "DOWN", "description": "Downregulated with age (GSE29801)"},
    "ARE": {"category": "Age", "direction": "mixed", "description": "NRF2/ARE targets (MSigDB + literature)"},
    "AMD-Macula-UP": {"category": "AMD", "direction": "UP", "description": "Upregulated in AMD RPE Macula (p<0.01)"},
    "AMD-Macula-DOWN": {"category": "AMD", "direction": "DOWN", "description": "Downregulated in AMD RPE Macula (p<0.01)"},
    "AMD-nonMacula-UP": {"category": "AMD", "direction": "UP", "description": "Upregulated in AMD RPE non-Macula (p<0.01)"},
    "AMD-nonMacula-DOWN": {"category": "AMD", "direction": "DOWN", "description": "Downregulated in AMD RPE non-Macula (p<0.01)"},
    "Senescence-PRO": {"category": "Senescence", "direction": "PRO", "description": "Pro-senescence (SenMayo + GO + SASP)"},
    "Senescence-ANTI": {"category": "Senescence", "direction": "ANTI", "description": "Anti-senescence (Cell cycle + E2F)"},
    "Redox-PRO": {"category": "Redox", "direction": "PRO", "description": "Pro-oxidative stress (HALLMARK + GO)"},
    "Redox-ANTI": {"category": "Redox", "direction": "ANTI", "description": "Antioxidant response (NRF2 targets)"},
    "Inflammation-PRO": {"category": "Inflammation", "direction": "PRO", "description": "Pro-inflammatory (HALLMARK + NFkB)"},
    "Inflammation-ANTI": {"category": "Inflammation", "direction": "ANTI", "description": "Anti-inflammatory (GO regulatory)"},
}

#%%
def ensure_dirs():
    """Create all output directories if they don't exist."""
    for name, path in OUTPUT_DIRS.items():
        path.mkdir(parents=True, exist_ok=True)
        
def get_dataset_output_dir(category: str, dataset: str, analysis_type: str) -> Path:
    """Get the output directory for a specific analysis."""
    # Handle both pair names (Age) and full names (Age-UP, Age-DOWN)
    if category in ["Age", "Age-UP", "Age-DOWN", "ARE"]:
        base = OUTPUT_DIRS["Age"]
    elif category.startswith("AMD"):
        base = OUTPUT_DIRS["AMD"]
    elif "Senescence" in category:
        base = OUTPUT_DIRS["Axes"] / "Senescence"
    elif "Redox" in category:
        base = OUTPUT_DIRS["Axes"] / "Redox"
    elif "Inflammation" in category:
        base = OUTPUT_DIRS["Axes"] / "Inflammation"
    else:
        # Fallback to Axes
        base = OUTPUT_DIRS["Axes"] / category
    
    out_dir = base / dataset / analysis_type
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir
