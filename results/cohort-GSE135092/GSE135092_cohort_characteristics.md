# GSE135092 Cohort Characteristics

**GEO Accession:** [GSE135092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135092)  
**Publication:** Orozco et al., Cell Reports (2020)  
**Title:** Integration of eQTL and a Single-Cell Atlas in the Human Eye Identifies Causal Genes for Age-Related Macular Degeneration

---

## 1. Data Type

**Platform:** Bulk RNA-seq (Illumina HiSeq 2500)  
**Sequencing Depth:** 30 million single-end 50bp reads per sample  
**Genome:** GRCh38 (human reference)  
**Quantification:** RPKM → nRPKM (DESeq2 normalized)  
**Genes Measured:** ~58,303 (Ensembl gene IDs)

### Notes on Data Availability
- **Raw reads:** Not provided (patient privacy)
- **Genotype/eQTL data:** Not in GEO (restricted)
- **Available:** Processed expression matrix and differential expression results

---

## 2. Sample Overview

| Metric | Value |
|:-------|------:|
| Total Samples | 537 |
| AMD Samples | 104 |
| Control Samples | 433 |
| Unique Donors | 537 (one sample per donor) |

---

## 3. Disease Status

Samples are classified as:
- **Control:** Donors with no history of ocular diseases (n=433)
- **AMD:** Donors previously diagnosed with age-related macular degeneration (n=104)

> **Note:** AMD staging (e.g., intermediate, neovascular, geographic atrophy) is mentioned in the publication but not available in the GEO metadata.

---

## 4. Tissue Types

Four tissue types are sampled, balanced across conditions:

| Tissue Type | Control | AMD | Total | % Control |
|:------------|--------:|----:|------:|----------:|
| RPE, Macula | 105 | 26 | 131 | 80.2% |
| RPE, non-Macula | 112 | 23 | 135 | 83.0% |
| Retina, Macula | 101 | 27 | 128 | 78.9% |
| Retina, non-Macula | 115 | 28 | 143 | 80.4% |
| **TOTAL** | **433** | **104** | **537** | **80.6%** |

### Tissue Distribution Within Each Group

**AMD (N=104):**
- RPE, Macula: 26 (25.0%)
- RPE, non-Macula: 23 (22.1%)
- Retina, Macula: 27 (26.0%)
- Retina, non-Macula: 28 (26.9%)

**Control (N=433):**
- RPE, Macula: 105 (24.2%)
- RPE, non-Macula: 112 (25.9%)
- Retina, Macula: 101 (23.3%)
- Retina, non-Macula: 115 (26.6%)

---

## 5. Age Distribution

**Age availability:**
- Control: 349 of 433 samples (80.6%)
- AMD: 52 of 104 samples (50.0%)

| Statistic | Control | AMD |
|:----------|--------:|----:|
| Min Age | 59 | 70 |
| Max Age | 95 | 98 |
| Mean Age | 80.3 | 81.8 |

### Age Bins

| Age Bin | Control (N) | Control (%) | AMD (N) | AMD (%) |
|:--------|------------:|------------:|--------:|--------:|
| <70 | 29 | 8.3% | 0 | 0.0% |
| 70-79 | 106 | 30.4% | 25 | 48.1% |
| 80-89 | 173 | 49.6% | 14 | 26.9% |
| 90+ | 41 | 11.7% | 13 | 25.0% |

---

## 6. Donor Matching

- **Within-donor matching:** NOT available (each sample is from a unique donor)
- **Control-case matching:** Group-level age matching only (similar age distributions)
- **Sample independence:** All 537 samples are biologically independent

---

## 7. Key Strengths

1. **Large sample size:** 537 samples (largest AMD RNA-seq cohort)
2. **Tissue resolution:** Macular vs peripheral, RPE vs retina
3. **Balanced design:** Consistent ~4:1 control:AMD ratio across all tissues
4. **Age matching:** Similar age distributions between groups

---

## 8. Limitations

1. **No AMD staging:** Cannot distinguish early vs late disease
2. **No genotype data:** Cannot assess genetic risk variants
3. **No raw reads:** Cannot reprocess from scratch
4. **Missing age data:** ~20% of controls, ~50% of AMD

---

## 9. Files in This Cohort

| File | Description |
|:-----|:------------|
| `GSE135092_DE_results.csv` | Differential expression (AMD vs Control, all genes) |
| `GSE135092_risk_expression.csv` | Expression of known AMD risk genes |
| `GSE135092_cohort_characteristics.md` | This document |

---

*Generated: January 2026*
