# Tannin-AMD: Investigating RPE Phenotypes and the Role of PRG4

## Project Overview
The primary goal of this project is to understand how Retinal Pigment Epithelium (RPE) cells develop phenotypes characteristic of Age-related Macular Degeneration (AMD), such as oxidative stress and inflammation (common in dry AMD). 

We are specifically investigating:
1.  **AMD-like Stress**: Identifying the transcriptomic signatures of RPE dysfunction.
2.  **PRG4 Regulation**: Determining if and how Lubricin (PRG4) can control, mitigate, or rescue these AMD-like phenotypes.

## Strategy
We integrate focused experimental data with large-scale functional genomics:
*   **Bulk RNA-seq**: Direct measurement of RPE response to oxidative stress (H2O2) and PRG4 treatment.
*   **Perturb-seq**: Leveraging genome-scale genotype-to-phenotype maps (Replogle et al. 2022) to identify genetic drivers and regulatory nodes involved in these stress pathways.

## Repository Structure
*   `data/RPE_cells/`: Bulk RNA-seq data (CTRL, H2O2, PRG4, and H2O2+PRG4 conditions). Detailed mapping in `planning-mds/RPE_cells_data.md`.
*   `data/external/perturbseq/`: Large-scale Perturb-seq datasets from Replogle et al. 2022.
*   `planning-mds/`: Technical documentation, including GPU environment specs and Figshare download guides.
*   `code/`: Analysis scripts and data processing pipelines.

## Environment
*   **Python**: 3.12.3
*   **Framework**: PyTorch 2.8.0+cu129
*   **Hardware**: 2 × NVIDIA L40 (48GB)
*   **Environment Path**: `~/venvs/torch`
*   Full details in `planning-mds/GPU_ENVIRONMENT.md`.
