# Task 7: Virtual PRG4 Screen Using K562 GWPS

## Objective
Apply the **transfer model** (Task 6) to the full K562 GWPS dataset (~11K KDs) to discover novel PRG4 proxies, antagonists, and AMD-mimicking genes beyond the essential gene set.

## Motivation
Current bridge analysis (RPE1 only) is limited to ~2,679 essential genes. The K562 GWPS has 4× more coverage. By predicting what these KDs would do in RPE cells, we can:
1. Identify non-essential PRG4 proxies (genes whose KD mimics PRG4 rescue)
2. Find non-essential PRG4 antagonists (genes opposing PRG4 effects)
3. Discover broader AMD-mimicking genes (oxidative stress pathways)

## Prerequisites
- **Task 6 must succeed**: Transfer model achieves reasonable accuracy (test ρ > 0.4)
- **Task 1 complete**: Optimal signature threshold determined

## Data Sources

### Bulk RNA-seq Signatures
- **Input**: `data/RPE_cells/code/RPE_gene pvals.xlsx`
- **Signatures** (using Task 1's optimal threshold):
  - H2O2 Stress (H2O2_vs_CTRL): log2FC values
  - PRG4 Baseline (PRG4_vs_CTRL): log2FC values  
  - PRG4 Rescue (H2O2PRG4_vs_H2O2): log2FC values

### Predicted RPE1 Profiles
- **Input**: `results/transfer-model/predicted_RPE1_from_K562_GWPS.h5ad`
  - Contains predicted RPE1 expression profiles for ~11K K562 GWPS KDs
  - Includes confidence scores

### Actual RPE1 Profiles (for comparison)
- **Input**: `data/external/perturbseq/rpe1_normalized_bulk_01.h5ad`
  - Used to compare novel hits with known hits

## Analysis Tasks

### 1. Harmonize Signatures and Predictions
- Match bulk RNA-seq DEGs to genes measured in predicted RPE1 profiles
- Extract common gene set for correlation
- Create signature vectors (log2FC for DEGs)

### 2. Calculate Bridge Correlations
For each of the ~11K predicted RPE1 profiles:

**H2O2 Stress Signature**:
- Spearman correlation between predicted profile and H2O2 log2FC
- Rank KDs by correlation
- **Top positives** = AMD-mimics (induce oxidative stress phenotype)

**PRG4 Baseline Signature**:
- Correlation with PRG4_vs_CTRL log2FC
- **Top positives** = PRG4 baseline mimics

**PRG4 Rescue Signature**:
- Correlation with H2O2PRG4_vs_H2O2 log2FC
- **Top positives** = PRG4 proxies (genes mimicking rescue)
- **Top negatives** = PRG4 antagonists (genes opposing rescue)

### 3. Filter by Confidence
- Use confidence scores from Task 6 predictions
- Flag or filter low-confidence predictions
- Report both "high-confidence hits" and "all hits"

### 4. Comparison with RPE1-Only Results
Compare novel hits (from K562 GWPS) with original RPE1 results:
- Do the top RPE1 hits (SNIP1, DBR1, ARL4D, RFC2/4, etc.) still rank highly?
  - If yes: Validates transfer model
  - If no: Model may be inaccurate; proceed with caution

### 5. Novelty Analysis
Identify genes that are:
- **Novel hits**: High correlation in predicted data, NOT in RPE1 essential set
- **Essential-only hits**: High in RPE1, not captured in K562 GWPS

Characterize each set:
- GO/KEGG enrichment
- Are novel hits biologically plausible for AMD/PRG4?

### 6. Stratification by Gene Type
Break down hits by:
- **Essential vs non-essential** (from DEG2 database or literature)
- **Druggable** (kinases, GPCRs, secreted proteins, etc.)
- **Previously linked to AMD** (GWAS hits, literature mining)

## Expected Outputs

### Correlation Results
- `results/virtual-screen/H2O2_Stress_K562_GWPS_bridge.csv`
  - Columns: KD_gene, Spearman_rho, Rank, Confidence, In_RPE1_dataset
  
- `results/virtual-screen/PRG4_Baseline_K562_GWPS_bridge.csv`
- `results/virtual-screen/PRG4_Rescue_K562_GWPS_bridge.csv`

### Novelty Analysis
- `results/virtual-screen/novel_AMD_mimics.csv` (top hits not in RPE1)
- `results/virtual-screen/novel_PRG4_proxies.csv`
- `results/virtual-screen/novel_PRG4_antagonists.csv`

### Enrichment
- `results/virtual-screen/novel_hits_enrichment.csv`
  - GO/KEGG enrichment for each category

### Visualizations
1. **Correlation distribution histogram**:
   - Compare RPE1-only vs K562 GWPS-predicted distributions

2. **Rank comparison scatter**:
   - X-axis: Rank in RPE1 (for overlapping genes)
   - Y-axis: Rank in K562 GWPS predictions
   - Validate transfer model quality

3. **Enrichment bar charts**:
   - Top pathways in novel vs essential hits

4. **Confidence vs Correlation scatter**:
   - Identify high-confidence + high-correlation hits

5. **Druggability/AMD GWAS annotation**:
   - Highlight clinically relevant hits

### Summary Report
- `results/virtual-screen/virtual_screen_summary.md`
  - Top 10-20 novel hits for each signature
  - Biological interpretation
  - Prioritization for experimental validation

## Success Criteria
- **Validation**: Top RPE1 hits also rank highly in predicted data (Spearman ρ > 0.7)
- **Novelty**: At least 10-20 high-confidence novel hits per signature
- **Biological plausibility**: Novel hits enrich for relevant pathways (oxidative stress, DNA repair, inflammation, etc.)

## Potential Follow-Up
- **Literature mining**: Check if novel hits have prior links to AMD, retinal disease, or neuroprotection
- **Drug repurposing**: Identify if any novel hits are druggable targets with existing compounds
- **Experimental validation**: Prioritize top novel hits for wet-lab validation (siRNA KD in RPE cells)

## Dependencies
- **Task 1 (Signature Robustness)**: Provides optimal signature threshold
- **Task 6 (Transfer Model)**: Must succeed to generate predictions

## Next Steps
- **Task 8 (Multi-Signature Integration)**: Plot novel hits in 2D/3D correlation space
- **Task 9 (Network Analysis)**: Build interaction networks around novel hits
- Discuss experimental validation with biology team
