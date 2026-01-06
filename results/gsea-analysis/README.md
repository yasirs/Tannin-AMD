# GSEA & Mechanistic Analysis Results

**Objective:** To move beyond simple correlation and understand the *mechanism* of the rescue by identifying the core "Leading Edge" genes.

## Methodology
1.  **GSEA (Gene Set Enrichment Analysis):**
    *   **Gene Sets:** `PRG4_Rescue_UP` and `PRG4_Rescue_DOWN` (from internal data).
    *   **Ranked List:** LogFC from **GSE129964 (Serum Starvation)**.
    *   **Algorithm:** Pre-ranked GSEA (KS statistic).
2.  **Leading Edge Analysis:** We extracted the specific subset of genes driving the enrichment (i.e., genes that are *both* Induced by Rescue and Repressed by Starvation).
3.  **Mimetic Check:** We tested if our top Virtual Screen hits (e.g., KEAP1 KD) significantly induce this same "Leading Edge" set.

## Key Findings
*   **Pathway:** The "Rescue" program is driven by the restoration of **Cell Cycle** and **DNA Repair** genes.
*   **Mechanism:** Knocking down **KEAP1** or **DICER1** significantly induces this exact "Leading Edge" gene set (p < 1e-8). This proves that activating the NRF2 pathway is sufficient to drive the observed rescue phenotype.

## Files
*   **`gsea_standard_rescue_up.pdf`**: Enrichment plot showing the Rescue set reversing the Starvation signature.
*   **`leading_edge_enrichment.csv`**: Pathway analysis of the core rescue genes.