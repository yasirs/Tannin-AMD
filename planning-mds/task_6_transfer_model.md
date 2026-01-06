# Task 6: Transfer Model Development

## Objective
Build a computational model to translate perturbation effects from K562 to RPE1, enabling use of the broader K562 GWPS dataset (~11K KDs) to infer RPE biology.

## Motivation
- **RPE1** has only ~2,679 essential gene KDs
- **K562 GWPS** has ~11,258 KDs (genome-wide coverage)
- **Gap**: ~8,500 additional perturbations not in RPE1

If we can learn a mapping function **K562 → RPE1**, we can:
1. Predict what these 8,500 KDs would do in RPE cells
2. Find PRG4 proxies/antagonists not in the essential gene set
3. Discover novel AMD-relevant biology

## Prerequisites
- **Task 4 must show sufficient concordance** (median ρ > 0.3-0.4)
- **Task 5 must validate** that conserved processes exist
- If concordance is too low, this task may not be viable

## Data Sources

### Training Data (Overlapping KDs)
- **RPE1**: `data/external/perturbseq/rpe1_normalized_bulk_01.h5ad`
- **K562 Essential**: `data/external/perturbseq/K562_essential_normalized_bulk_01.h5ad`
- **Overlap**: ~2,000-2,500 perturbations present in both

### Target Data (K562 GWPS)
- **K562 GWPS**: `data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad`
- Contains the additional ~8,500 KDs not in RPE1

### Filtered Gene Set
- **Task 3 output**: `results/baseline-expression/commonly_expressed_genes.csv`
- Model should operate only on genes expressed in both cell types

## Modeling Approaches

### Option 1: Gene-Level Linear Transfer
**Simplest approach**: For each gene, fit a linear model.

For gene *g*:
```
RPE1_expression[g] = α[g] · K562_expression[g] + β[g]
```

- Train on overlapping KDs (N~2000)
- Learn α and β for each gene
- Apply to K562 GWPS to predict RPE1 effects

**Pros**: Simple, interpretable, fast
**Cons**: Assumes linear relationship, ignores gene-gene interactions

### Option 2: Ridge/Lasso Regression
**Multi-gene model**: Predict each RPE1 gene's response using multiple K562 genes.

For RPE1 gene *g*:
```
RPE1[g] = Σ w[g,k] · K562[k] + intercept
```

- Regularization prevents overfitting
- Can capture some cross-talk between genes
- Interpretable coefficients

**Pros**: More expressive than linear, still interpretable
**Cons**: Computationally heavier, may need feature selection

### Option 3: Neural Network (Multi-Task Learning)
**Deep learning approach**: Shared encoder + cell-type-specific decoders.

Architecture:
- **Input**: K562 expression profile (N_genes dimensions)
- **Shared encoder**: Learns cell-type-invariant perturbation features
- **RPE1 decoder**: Predicts RPE1 expression from shared features
- **K562 decoder**: Reconstructs K562 expression (for regularization)

Train on overlapping KDs, apply encoder → RPE1 decoder to K562 GWPS.

**Pros**: Can capture complex non-linear relationships
**Cons**: Black-box, requires hyperparameter tuning, prone to overfitting

### Option 4: Pathway-Level Transfer
**Biological approach**: Map to functional modules first, then transfer.

1. Decompose expression profiles into **pathway scores** (e.g., PCA, NMF, or curated gene sets)
2. Learn pathway-level transfer: `RPE1_pathway[i] = f(K562_pathway[i])`
3. More robust to cell-type differences if pathways are conserved

**Pros**: Biologically interpretable, robust
**Cons**: Requires good pathway definitions, may lose gene-level resolution

## Recommended Initial Approach
Start with **Option 1 + Option 2**:
1. Fit simple gene-level linear models (baseline)
2. Fit Ridge regression with cross-validation (improved baseline)
3. Evaluate both on held-out test set
4. If performance is insufficient, proceed to Option 3 or 4

## Analysis Tasks

### 1. Data Preparation
- Merge RPE1 and K562 Essential datasets by overlapping KDs
- Filter to commonly expressed genes (from Task 3)
- Split overlapping KDs into **train (80%) / test (20%)** sets
  - Stratify by concordance level (ensure test set has both high and low concordance KDs)

### 2. Model Training
For each approach (Linear, Ridge, potentially NN):
- Train on training set
- Tune hyperparameters (e.g., Ridge alpha) via cross-validation
- Document model coefficients/parameters

### 3. Model Evaluation
On held-out test set:
- **Per-KD correlation**: For each test KD, correlate predicted RPE1 vs actual RPE1 profile
  - Metric: Spearman ρ
  - Plot distribution of test correlations

- **Gene-level accuracy**: For each gene, how well are its responses predicted?
  - Metric: R² or Spearman ρ across KDs

- **Pathway-level accuracy**: Do predicted profiles enrich for correct pathways?

### 4. Comparison to Naive Baseline
- **Baseline 1**: Just use K562 expression as-is (no transfer)
- **Baseline 2**: Use average RPE1 response (ignore K562)

Model should outperform both.

### 5. Apply to K562 GWPS
For each of the ~8,500 novel KDs in K562 GWPS:
- Extract K562 expression profile
- Apply trained model to predict RPE1 expression profile
- Save predictions

### 6. Uncertainty Quantification
For predictions:
- Calculate confidence intervals (bootstrap or ensemble methods)
- Flag low-confidence predictions (may exclude from downstream analysis)

## Expected Outputs

### Model Files
- `models/transfer_model_linear.pkl` (or .pt for PyTorch)
- `models/transfer_model_ridge.pkl`
- `models/transfer_model_best.pkl` (selected model)

### Evaluation Results
- `results/transfer-model/model_performance.csv`
  - Columns: Model, Train_Corr, Test_Corr, Test_R2, etc.
  
- `results/transfer-model/per_gene_accuracy.csv`
  - Gene-wise performance metrics

- `results/transfer-model/test_predictions.csv`
  - Actual vs predicted for test set KDs

### Predictions
- `results/transfer-model/predicted_RPE1_from_K562_GWPS.h5ad`
  - AnnData with predicted RPE1 profiles for all K562 GWPS KDs
  - Include confidence scores in `.obs`

### Visualizations
1. **Model comparison bar chart**: Test correlation for each approach
2. **Scatter plot**: Actual vs predicted RPE1 expression (select KDs)
3. **Gene-wise accuracy heatmap**: Which genes are well-predicted?
4. **Concordance vs predictability**: Do high-concordance KDs predict better?

### Summary Report
- `results/transfer-model/transfer_model_summary.md`

## Success Criteria
- **Test correlation ρ > 0.4**: Model captures meaningful signal
- **Outperforms naive baseline**: Transfer adds value over direct use
- **High-concordance KDs predict well**: Validates consistency with Task 4

**If success criteria not met**: May need more sophisticated models (Option 3/4) or reconsider approach.

## Dependencies
- **Task 3 (Baseline Expression)**: Provides gene filtering
- **Task 4 (Cell-Type Concordance)**: Must show sufficient conservation
- **Task 5 (Validation)**: Should pass before investing in complex models

## Next Steps
- **Task 7 (Virtual PRG4 Screen)**: Apply model to find PRG4 proxies in K562 GWPS
- **Task 8 (Multi-Signature Integration)**: Use predicted profiles for 2D/3D analysis
- **Task 9 (Network Analysis)**: Build networks using predicted RPE1 responses
