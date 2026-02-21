# GSE29801 Cohort Metadata Documentation

**GEO Accession:** [GSE29801](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29801)  
**Publication:** Newman et al. 2012, *Genome Medicine*  
**Title:** Systems-level analysis of age-related macular degeneration reveals global biomarkers and phenotype-specific functional networks

---

## 1. Platform & Data Type

| Attribute | Value |
|:----------|:------|
| **Platform** | Affymetrix Human Genome U133 Plus 2.0 Array (GPL4133) |
| **Technology** | Microarray (NOT RNA-seq) |
| **Probes** | ~54,000 probes → ~20,000 unique genes after collapsing |
| **Data Type** | Normalized log2 intensity values |

### What's Included
- Gene expression data only (mRNA levels)
- Clinical/phenotype metadata (see below)

### What's NOT Included
- ❌ Genetic/genotype data (CFH, ARMS2 variants not available)
- ❌ Proteomics, metabolomics, or other omics
- ❌ Single-cell resolution

---

## 2. Sample Summary

| Category | Count |
|:---------|------:|
| **Total Samples** | 293 |
| **Control (Normal)** | 151 |
| **AMD** | 142 |

---

## 3. Tissue Types

Samples are from **2 anatomical regions** × **2 tissue layers** = 4 tissue types:

| Tissue Type | Controls | AMD | Total |
|:------------|:--------:|:---:|------:|
| Macular RPE-choroid | 50 | 41 | 91 |
| Extramacular RPE-choroid | 46 | 38 | 84 |
| Macular Retina | 28 | 32 | 60 |
| Extramacular Retina | 27 | 31 | 58 |
| **Total** | **151** | **142** | **293** |

### Key Points
- **RPE-choroid** (175 samples): Contains RPE cells + choroidal vasculature — **most relevant for RPE biology**
- **Retina** (118 samples): Neural retina (photoreceptors, ganglion cells) — less relevant for RPE-specific analyses
- **Macular**: Central retina where AMD pathology is most severe
- **Extramacular**: Peripheral retina, typically less affected

> **Note for RPE analyses:** Filter to RPE-choroid samples only (N=175).

---

## 4. AMD Clinical Staging

The AMD samples have **detailed clinical classification**:

| AMD Classification | N Samples | Stage Description |
|:-------------------|----------:|:------------------|
| Dry AMD | 61 | Early/intermediate AMD with drusen |
| MD1 (Macular Drusen Stage 1) | 20 | Early AMD |
| MD2 (Macular Drusen Stage 2) | 13 | Intermediate AMD |
| CNV (Choroidal Neovascularization) | 16 | Wet AMD (neovascular) |
| GA (Geographic Atrophy) | 8 | Advanced dry AMD |
| GA/CNV (Mixed) | 12 | Both GA and CNV present |
| Clinical AMD diagnosis | 12 | Unspecified AMD |
| **Normal** | **151** | **Healthy controls** |

### Staging Hierarchy (simplified)
```
Normal → Early (MD1, MD2) → Intermediate (Dry AMD) → Advanced (GA, CNV, GA/CNV)
```

---

## 5. Age Distribution

| Group | Mean Age | Median Age | Range |
|:------|:--------:|:----------:|:-----:|
| Controls | 63.1 years | 68 years | 9–93 |
| AMD | 80.8 years | 83 years | 43–101 |

> ⚠️ **Important:** AMD patients are ~18 years older on average than controls. Age is a significant confounder and should be considered in differential expression analyses.

---

## 6. Other Metadata

| Field | Available | Notes |
|:------|:---------:|:------|
| Gender/Sex | ✅ Yes | Male and Female |
| Age | ✅ Yes | In years |
| Ocular disease status | ✅ Yes | Normal vs AMD |
| AMD classification | ✅ Yes | Detailed staging (see above) |
| Tissue location | ✅ Yes | Macular vs Extramacular |
| Tissue type | ✅ Yes | RPE-choroid vs Retina |
| CFH genotype | ❌ No | Not included |
| ARMS2 genotype | ❌ No | Not included |
| Other SNPs | ❌ No | Not included |

---

## 7. Data Files

| File | Description |
|:-----|:------------|
| `GSE29801_series_matrix.txt.gz` | Full metadata + expression matrix |
| `GSE29801_expression_matrix.csv.gz` | Pre-processed expression matrix |
| `GSE29801_amd_signature.csv` | DE results (AMD vs Control) |
| `GPL4133_old_annotations.txt.gz` | Probe-to-gene symbol mapping |

---

## 8. Recommended Subsets for Analysis

### A. RPE-focused analysis
- Filter to **RPE-choroid** samples only
- N = 96 controls + 79 AMD = **175 samples**

### B. Age-matched analysis
- Filter controls to **age > 70** for better matching with AMD patients
- Or include age as a covariate in statistical models

### C. Disease progression analysis
- Group by staging:
  - **Early**: MD1, MD2
  - **Intermediate**: Dry AMD
  - **Advanced**: GA, CNV, GA/CNV
- Compare gene expression across stages

---

*Documentation generated: January 11, 2026*
