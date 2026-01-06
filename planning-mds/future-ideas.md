# Future Analysis Ideas: Clinical Relevance for PRG4 Therapy

This document outlines potential analyses to strengthen the clinical translatability of the PRG4-AMD therapeutic hypothesis, beyond the current GEO disease correlation work.

---

## 1. GWAS Integration: PRG4 Targets × AMD Risk Loci

### Concept
Demonstrate that PRG4-mimetic genes (KEAP1, DICER1, UBR5) or leading edge genes are enriched near AMD GWAS hits, providing human genetic evidence of causality.

### Data Sources
| Resource | URL | Description |
|----------|-----|-------------|
| GWAS Catalog | https://www.ebi.ac.uk/gwas/ | Search "age-related macular degeneration" for summary statistics |
| OpenTargets Genetics | https://genetics.opentargets.org/ | Pre-computed variant-to-gene mappings |
| GTEx Portal | https://gtexportal.org/ | eQTL data for retina (if available) or related tissues |
| AMD GWAS Summary Stats | PMID: 26691988 (Fritsche 2016) | 34 independent loci, largest AMD GWAS |

### Analysis Steps
1. Download AMD GWAS summary statistics from GWAS Catalog or original publication
2. Use **MAGMA** (https://ctg.cncr.nl/software/magma) for gene-set enrichment testing
3. Alternatively, use **SNPsea** or **LDSC** for partitioned heritability
4. Simple overlap: count how many PRG4-rescue genes are within 500kb of a GWAS lead SNP
5. Use OpenTargets to check if KEAP1/DICER1/UBR5 have AMD genetic associations

### Visualization
- Manhattan plot with PRG4-rescue genes annotated as vertical lines
- Venn diagram: GWAS genes ∩ PRG4-mimetics ∩ Leading Edge
- Bar plot of enrichment p-values for different gene sets

---

## 2. Drug Repurposing via Connectivity Map (CMap/L1000)

### Concept
Query the PRG4-rescue signature against large-scale drug perturbation databases to identify FDA-approved compounds that mimic PRG4's effect.

### Data Sources
| Resource | URL | Description |
|----------|-----|-------------|
| CLUE.io | https://clue.io/ | Broad Institute L1000 query interface |
| LINCS L1000 | GEO: GSE92742 | Raw L1000 data if programmatic access needed |
| CREEDS | http://amp.pharm.mssm.edu/CREEDS/ | Curated drug signatures from GEO |
| DrugMatrix | https://ntp.niehs.nih.gov/data/drugmatrix | Toxicogenomics signatures |

### Analysis Steps
1. Prepare PRG4-rescue signature: top 150 UP genes and top 150 DOWN genes (L1000 format)
2. Submit to CLUE.io query tool (free account required)
3. Retrieve ranked drugs by connectivity score (positive = mimics PRG4)
4. Cross-reference top hits with clinical trial databases (ClinicalTrials.gov)
5. Validate hits are biologically plausible (e.g., NRF2 activators like DMF, Bardoxolone)

### Visualization
- Radial "drug wheel" showing top 10-20 mimetic compounds
- Table: Drug name, mechanism, CMap score, clinical status
- Network: Drug → Target → PRG4-mimetic effect

### Expected Hits
- **Sulforaphane** (NRF2 activator, natural compound)
- **Dimethyl fumarate (DMF)** (FDA-approved for MS, NRF2 activator)
- **Bardoxolone methyl** (clinical trials for kidney disease, NRF2)

---

## 3. Target Druggability and Protein-Level Validation

### Concept
Confirm that PRG4-mimetic targets are "druggable" and have existing pharmacological tools.

### Data Sources
| Resource | URL | Description |
|----------|-----|-------------|
| DGIdb | https://www.dgidb.org/ | Drug-Gene Interaction database |
| OpenTargets Platform | https://platform.opentargets.org/ | Tractability scores and known drugs |
| ChEMBL | https://www.ebi.ac.uk/chembl/ | Bioactivity data for compounds |
| UniProt | https://www.uniprot.org/ | Protein annotations |
| Human Protein Atlas | https://www.proteinatlas.org/ | Protein expression in retina/RPE |

### Analysis Steps
1. Query DGIdb API with top 50 PRG4-mimetic genes
2. Extract: number of known drug interactions, drug categories, approval status
3. For KEAP1: identify NRF2 activators that work via KEAP1 inhibition
4. Check Human Protein Atlas for protein expression in eye tissues
5. Use OpenTargets to get "tractability bucket" scores

### Visualization
- Druggability heatmap: genes × druggability metrics
- Sankey diagram: PRG4 rescue → Target genes → Existing drugs
- Bar plot of targets by drug development stage

---

## 4. Clinical Cohort: Disease Staging Analysis

### Concept
Score AMD patient samples by PRG4-rescue signature and test correlation with disease stage or progression.

### Data Sources
| GEO ID | Description | Samples |
|--------|-------------|---------|
| GSE29801 | AMD vs control RPE/choroid (Newman 2012) | 118 samples, macula + periphery |
| GSE50195 | Geographic atrophy RPE | 8 samples |
| GSE135092 | Human retina scRNA-seq with AMD | ~50k cells |
| GSE103186 | RPE from AMD donors | 6 AMD, 6 control |
| GSE181469 | Drusen transcriptomics | Drusen deposits |

### Analysis Steps
1. Download normalized expression matrices from GEO
2. Calculate PRG4-rescue module score for each sample (using ssGSEA or GSVA)
3. Stratify samples by: AMD stage (early, GA, wet), macular vs peripheral, age
4. Statistical test: Kruskal-Wallis or linear regression for score ~ stage
5. For scRNA-seq: score individual cells, compare RPE cells in AMD vs control

### Visualization
- **Boxplot**: PRG4-rescue score across disease stages
- **Kaplan-Meier**: If progression data available, stratify by signature score
- **Violin plots**: Score distribution by tissue region (macula vs periphery)
- **ROC curve**: Can the signature distinguish AMD from control?

---

## 5. Single-Cell Atlas: Cell-Type Specificity

### Concept
Map PRG4 receptor expression and rescue signature onto human retina single-cell atlases to show tissue-specific relevance.

### Data Sources
| Dataset | GEO/Source | Description |
|---------|------------|-------------|
| Lukowski 2019 | GSE137537 | Human retina atlas, 20k cells |
| Voigt 2019 | GSE130636 | Human macula scRNA-seq |
| Menon 2019 | GSE137846 | Fovea vs periphery |
| Sridhar 2020 | GSE142449 | Developing human retina |
| CellxGene | cellxgene.cziscience.com | Curated retina datasets |
| Tabula Sapiens | tabula-sapiens.ds.czbiohub.org | Multi-organ atlas |

### Analysis Steps
1. Download processed AnnData/Seurat objects from CellxGene
2. Identify cell types: RPE, photoreceptors, Müller glia, microglia, endothelial, pericytes
3. Score cells using:
   - PRG4 receptor expression: **CD44, TLR2, TLR4, LGALS3**
   - PRG4-rescue gene set: AUCell or AddModuleScore
4. Compare scores across cell types (ANOVA or Wilcoxon)
5. Identify RPE as the "PRG4-responsive" population

### Visualization
- **UMAP colored by CD44/TLR2 expression**
- **UMAP colored by PRG4-rescue score**
- **Bar plot**: Mean rescue score by cell type
- **Dot plot**: Receptor expression × cell type

---

## 6. Multi-Stressor Meta-Analysis

### Concept
Expand beyond H2O2 and serum starvation to show PRG4 is a universal cytoprotectant across diverse RPE stressors.

### Data Sources
| GEO ID | Stressor | Cell Type |
|--------|----------|-----------|
| GSE129964 | Serum starvation | ARPE-19 (already done) |
| GSE67899 | Cigarette smoke extract | ARPE-19 |
| GSE60436 | A2E (lipofuscin component) | ARPE-19 |
| GSE36139 | Oxidative stress (multiple) | Cell lines |
| GSE148387 | Blue light | RPE |
| GSE75968 | 4-HNE (lipid peroxidation) | ARPE-19 |
| GSE115828 | Complement attack | iPSC-RPE |
| GSE89189 | Sub-RPE deposits | Primary RPE |

### Analysis Steps
1. Download DEG lists or normalized counts from each dataset
2. Define "stress signature" for each: top UP and DOWN genes
3. Correlate each stress signature with PRG4-rescue signature
4. Meta-analysis: combine correlations with random effects model
5. Test heterogeneity (is PRG4 universally protective?)

### Visualization
- **Forest plot**: Correlation coefficient ± 95% CI for each stressor
- **Heatmap**: Stress signatures × PRG4 rescue (showing anti-correlation)
- **Radar/spider plot**: PRG4 protection score across stressor types

---

## 7. Network Medicine: PRG4 as Hub Modulator

### Concept
Build a protein-protein interaction network of PRG4-rescue genes and show that KEAP1/DICER1 are central regulatory hubs.

### Data Sources
| Resource | URL | Description |
|----------|-----|-------------|
| STRING | https://string-db.org/ | PPI database, API available |
| BioGRID | https://thebiogrid.org/ | Curated interactions |
| Reactome | https://reactome.org/ | Pathway-level networks |
| NDEx | https://www.ndexbio.org/ | Network repository |

### Analysis Steps
1. Export PRG4-rescue leading edge genes (top 200-500)
2. Query STRING API for interactions (score > 700)
3. Import network into **Cytoscape** or **igraph**
4. Calculate centrality metrics: betweenness, degree, PageRank
5. Identify modules using Louvain clustering
6. Overlay KEAP1/DICER1/UBR5 and PRG4 receptors on network

### Visualization
- **Network graph**: Nodes colored by function (stress, metabolism, signaling)
- **Hub highlight**: KEAP1/DICER1 as large, central nodes
- **Module layout**: Functional clusters labeled

---

## 8. Complement Pathway Deep Dive

### Concept
CFH (Complement Factor H) is the #1 AMD risk gene and is restored by PRG4. Deep analysis of complement pathway modulation.

### Data Sources
| Resource | Description |
|----------|-------------|
| ImmPort | Curated immune gene lists |
| KEGG hsa04610 | Complement and coagulation cascades |
| Literature | CFH, C3, CFB, CFI alleles in AMD |

### Analysis Steps
1. Extract all complement genes from KEGG pathway
2. Calculate rescue score for each complement gene
3. Test pathway-level enrichment in PRG4 rescue
4. Check if PRG4-mimetics (KEAP1 KD etc.) also modulate complement
5. Literature review: is NRF2 known to regulate complement?

### Visualization
- **Pathway diagram**: KEGG complement cascade with fold-change overlay
- **Bar plot**: Complement genes ranked by rescue magnitude
- **Correlation scatter**: CFH expression vs PRG4-rescue signature

---

## Priority Ranking

| Rank | Analysis | Impact | Effort | Visual Appeal |
|------|----------|--------|--------|---------------|
| 🥇 | CMap Drug Repurposing | ⭐⭐⭐ (Actionable) | Low | ⭐⭐⭐ |
| 🥈 | GWAS Enrichment | ⭐⭐⭐ (Genetic proof) | Medium | ⭐⭐ |
| 🥉 | Single-Cell Atlas | ⭐⭐ (Mechanism) | Medium | ⭐⭐⭐ |
| 4 | Clinical Cohort Staging | ⭐⭐⭐ (Prognostic) | High | ⭐⭐ |
| 5 | Multi-Stressor Meta | ⭐⭐ (Robustness) | Medium | ⭐⭐ |
| 6 | Network Medicine | ⭐ (Supplementary) | Low | ⭐⭐⭐ |
| 7 | Target Druggability | ⭐⭐ (Practical) | Low | ⭐⭐ |
| 8 | Complement Deep Dive | ⭐⭐ (Mechanism) | Low | ⭐⭐ |

---

## Quick Wins (Can do this week)

1. **CLUE.io query** — Just need to format gene list and submit
2. **DGIdb druggability** — API query, 1-hour task
3. **GSE29801 staging** — Classic AMD dataset, straightforward analysis
4. **STRING network** — Export genes, query API, visualize in Cytoscape

---

*Document created: January 6, 2026*
