# Perturb-seq Datasets (Replogle et al. 2022)

This directory contains processed pseudo-bulk datasets from the study:
**"Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq"**
Replogle, J. M., et al. (2022). *Cell*, 185(14), 2559-2575.e28.
[DOI: 10.1016/j.cell.2022.05.013](https://doi.org/10.1016/j.cell.2022.05.013)

## Dataset Overview
The data consists of AnnData (`.h5ad`) files for three Perturb-seq experiments:
1.  **K562 Genome-wide Perturb-seq (GWPS)**: Sampled at day 8 post-transduction.
2.  **K562 Essential-scale Perturb-seq**: Sampled at day 6 post-transduction.
3.  **RPE1 Essential-scale Perturb-seq**: Sampled at day 7 post-transduction.

### File Naming Convention
Files are named following the pattern: `$cell_line_$scale_$type_bulk_01.h5ad`

*   `raw`: Raw pseudo-bulk expression data for genes expressed at >0.01 UMI per cell.
*   `normalized`: gemgroup Z-normalized pseudo-bulk expression data.

## Files in this Directory

| Filename | Description | Size |
| :--- | :--- | :--- |
| `K562_gwps_raw_bulk_01.h5ad` | K562 Genome-wide, Raw pseudo-bulk | 358 MB |
| `K562_gwps_normalized_bulk_01.h5ad` | K562 Genome-wide, Z-normalized pseudo-bulk | 358 MB |
| `K562_essential_raw_bulk_01.h5ad` | K562 Essential-scale, Raw pseudo-bulk | 77 MB |
| `K562_essential_normalized_bulk_01.h5ad` | K562 Essential-scale, Z-normalized pseudo-bulk | 77 MB |
| `rpe1_raw_bulk_01.h5ad` | RPE1 Essential-scale, Raw pseudo-bulk | 91 MB |
| `rpe1_normalized_bulk_01.h5ad` | RPE1 Essential-scale, Z-normalized pseudo-bulk | 91 MB |
| `replogle_20029387.json` | Metadata and download links from Figshare | 12 KB |

## Internal Data Organization
The datasets follow a consistent AnnData structure:

### Observations (`.obs`) - Perturbations
Each row representing a pseudo-bulk population of cells targeted by a specific CRISPR guide (perturbation).
*   **Index (`gene_transcript`)**: Unique identifier for the perturbation, typically formatted as `[ID]_[Target_Gene]_[Target_Site]_[Gene_ID]` (e.g., `10005_ZBTB4_P1_ENSG00000174282`).
*   **Key Metadata Columns**:
    *   `UMI_count_unfiltered`: Total UMIs before filtering.
    *   `num_cells_unfiltered`/`filtered`: Cell counts for the pseudo-bulk aggregate.
    *   `fold_expr` / `pct_expr`: Expression metrics relative to controls.
    *   `cnv_score_z`: Copy number variation score (Z-normalized).
    *   `energy_test_p_value`, `anderson_darling_counts`, `mann_whitney_counts`: Statistical tests for perturbation effects.

### Variables (`.var`) - Measured Genes
Each column represents a gene whose expression was measured.
*   **Index (`gene_id`)**: Ensembl Gene ID.
*   **Key Metadata Columns**:
    *   `gene_name`: Common gene symbol (e.g., `NOC2L`).
    *   `mean`, `std`, `cv`: Basic expression statistics across the dataset.

### Matrix (`.X`)
Contains the actual expression values. For `normalized` files, these are gemgroup Z-normalized pseudo-bulk values.

## Experiment Scales
| Dataset | Number of Perturbations (obs) | Number of Measured Genes (var) |
| :--- | :--- | :--- |
| K562 GWPS (Genome-wide) | 11,258 | 8,248 |
| RPE1 Essential | 2,679 | 8,749 |
| K562 Essential | 2,285 | 8,563 |

## Source
Data was retrieved from Figshare+: [Dataset Link](https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387)
