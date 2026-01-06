# Virtual PRG4 Screen Results (Task 7)

**Objective:** To identify molecular targets that, when perturbed, replicate the beneficial "Rescue" phenotype induced by PRG4.

## Methodology
1.  **Query Signature:** The "PRG4 Rescue" signature (Genes DE in Rescue vs Stress, FDR < 0.05).
2.  **Database:** **K562 Genome-Wide Perturb-seq** (~11,000 gene knockdowns).
3.  **Screening Metric:** **Spearman Correlation** between the Query Signature and each Knockdown Profile.
    *   *Positive Correlation:* The KD mimics the drug (potential **agonist**).
    *   *Negative Correlation:* The KD opposes the drug (potential **antagonist**).
4.  **Filtering:** Hits were filtered to include only genes **expressed in RPE1** (Baseline UMI > 0.05) to ensure clinical relevance.

## Key Findings
*   **Top Mimetic:** **KEAP1**.
    *   *Mechanism:* KEAP1 inhibits NRF2. Knocking it down activates NRF2.
    *   *Conclusion:* PRG4 likely acts by activating the NRF2 antioxidant pathway.
*   **Top Antagonists:** **ATR, HSPA5**.
    *   *Mechanism:* DNA Damage Response and ER Stress.
    *   *Conclusion:* Exacerbating cellular stress prevents the rescue.

## Pathway Validation (Task 5)
*   **Method:** KEGG Enrichment on the Top 100 hits.
*   **Result:** Top mimetics are significantly enriched for **Ubiquitin Mediated Proteolysis** (p < 0.03), validating the KEAP1-NRF2 link.

## Files
*   **`prg4_virtual_screen_results.csv`**: The complete ranked list of potential targets.
*   **`enrichment_top_mimetics.csv`**: Biological pathways enriched among the top hits.