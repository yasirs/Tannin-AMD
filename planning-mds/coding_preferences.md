# Coding Preferences

## Language & Environment
- **Primary Language**: R
- **R Path**: `/usr/bin/R`
- **Visualization Library**: ggplot2
- **Python Environment**: `~/venvs/torch`
- **Package Management**: 
  - Install new packages if needed.
  - **Do not** uninstall packages belonging to others.

## Visualization Standards
- **Font**: Palatino Linotype
- **Figure Size**: Aim for small physical dimensions (inches).
- **File Formats**: 
  - PDF
  - TIFF (Compression = LZW+P, Resolution = 300 dpi)
- **Theme**: Use `theme_classic()` by default unless another theme is specifically required.
- **Axes Text**: Always explicitly set `axis.text = element_text(color = "black")`.

### Heatmaps (ggplot2)
- **Dendrograms**: Do not plot unless specifically needed.
- **Gene Heatmaps**: Input data should be TPM converted to **gene-wise z-scores**.
- **Color Palettes**:
  - **Diverging Data**: BrBG (Brown for high values).
  - **Positive 'One-Sided' Data**: YlOrB (Brown for high values).

## Workflow & Coding Style
- **Code Chunks**: Separated by `#%%` (no space between `#` and `%%`).
- **Comments**: Use standard `# comment` format. Avoid lines with many repeated hashes (e.g. `#######`).
- **Caching**:
  - For computationally heavy analyses, store cached results in a dedicated subfolder.
  - Prevent recomputation for minor downstream changes.
  - **Invalidation**: Delete cache to force recomputation if logic changes.