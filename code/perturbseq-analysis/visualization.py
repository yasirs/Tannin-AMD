"""
Visualization Module

Generate publication-ready figures following coding_preferences.md.
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import sys

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))
from config import VIZ_SETTINGS, OUTPUT_DIRS

#%%
def setup_style():
    """Set up matplotlib style to match coding preferences."""
    plt.rcParams.update({
        # Font settings
        "font.family": "serif",
        "font.serif": ["Palatino", "P052", "DejaVu Serif", "Times New Roman"],
        "font.size": VIZ_SETTINGS["font_size"],
        
        # Figure settings
        "figure.figsize": (VIZ_SETTINGS["fig_width"], VIZ_SETTINGS["fig_height"]),
        "figure.dpi": VIZ_SETTINGS["fig_dpi"],
        "figure.facecolor": "white",
        "figure.edgecolor": "white",
        
        # Axes settings
        "axes.facecolor": "white",
        "axes.edgecolor": "black",
        "axes.linewidth": 0.8,
        "axes.labelcolor": "black",
        "axes.labelsize": VIZ_SETTINGS["font_size"],
        
        # Tick settings
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.labelsize": VIZ_SETTINGS["font_size"] - 1,
        "ytick.labelsize": VIZ_SETTINGS["font_size"] - 1,
        
        # Legend
        "legend.frameon": False,
        "legend.fontsize": VIZ_SETTINGS["font_size"] - 1,
        
        # Grid
        "axes.grid": False,
        
        # Spines
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

# Apply style on import
setup_style()

#%%
def save_figure(fig, filepath: Path, save_pdf: bool = True, save_tiff: bool = True):
    """Save figure in PDF and TIFF formats."""
    filepath = Path(filepath)
    
    if save_pdf:
        pdf_path = filepath.with_suffix('.pdf')
        fig.savefig(pdf_path, format='pdf', bbox_inches='tight', dpi=VIZ_SETTINGS["fig_dpi"])
        print(f"    Saved: {pdf_path.name}")
    
    if save_tiff:
        tiff_path = filepath.with_suffix('.tiff')
        fig.savefig(tiff_path, format='tiff', bbox_inches='tight', 
                   dpi=VIZ_SETTINGS["fig_dpi"], pil_kwargs={'compression': 'tiff_lzw'})
        print(f"    Saved: {tiff_path.name}")
    
    plt.close(fig)

#%%
def plot_mimetic_distribution(
    df: pd.DataFrame,
    score_col: str = "mimetic_score",
    title: str = "Mimetic Score Distribution",
    output_path: Path = None,
    color: str = None
):
    """
    Plot distribution of mimetic scores.
    
    Shows histogram with density overlay and marks for top/bottom.
    """
    scores = df[score_col].values
    
    fig, ax = plt.subplots(figsize=(5, 4))
    
    # Histogram
    color = color or "#1f77b4"
    ax.hist(scores, bins=50, density=True, alpha=0.7, color=color, edgecolor='white', linewidth=0.5)
    
    # Add vertical lines for mean and zero
    ax.axvline(0, color='gray', linestyle='--', linewidth=1, label='Zero')
    ax.axvline(np.mean(scores), color='red', linestyle='-', linewidth=1, label=f'Mean: {np.mean(scores):.3f}')
    
    # Labels
    ax.set_xlabel('Mimetic Score')
    ax.set_ylabel('Density')
    ax.set_title(title)
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def plot_convergence_scatter(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    significant_mask: np.ndarray = None,
    highlight_genes: list = None,
    title: str = "Convergence Scatter",
    xlabel: str = None,
    ylabel: str = None,
    output_path: Path = None,
    category_color: str = None
):
    """
    Plot convergence scatter - simple colored dots only.
    
    Per user request: significant points colored, non-significant gray, NO shapes.
    """
    x = df[x_col].values
    y = df[y_col].values
    
    fig, ax = plt.subplots(figsize=(5, 5))
    
    # Color scheme
    color_sig = category_color or VIZ_SETTINGS["color_significant"]
    color_nonsig = VIZ_SETTINGS["color_nonsignificant"]
    
    if significant_mask is not None:
        # Non-significant points (gray)
        ax.scatter(x[~significant_mask], y[~significant_mask], 
                  c=color_nonsig, s=10, alpha=0.3, rasterized=True)
        
        # Significant points (colored)
        ax.scatter(x[significant_mask], y[significant_mask],
                  c=color_sig, s=15, alpha=0.8, rasterized=True)
    else:
        # All points with same color
        ax.scatter(x, y, c=color_sig, s=10, alpha=0.5, rasterized=True)
    
    # Highlight specific genes if provided
    if highlight_genes:
        for gene in highlight_genes:
            mask = df["target_gene"].str.upper() == gene.upper()
            if mask.any():
                ax.scatter(x[mask], y[mask], c='red', s=50, edgecolors='black', 
                          linewidths=0.5, zorder=10)
                # Add label
                for idx in np.where(mask)[0]:
                    ax.annotate(gene, (x[idx], y[idx]), fontsize=7,
                               xytext=(5, 5), textcoords='offset points')
    
    # Add reference lines
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # Labels
    ax.set_xlabel(xlabel or x_col)
    ax.set_ylabel(ylabel or y_col)
    ax.set_title(title)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def plot_enrichment_distribution(
    df: pd.DataFrame,
    score_col: str = "enrichment_score",
    group_col: str = None,
    title: str = "Enrichment Score Distribution",
    output_path: Path = None
):
    """
    Plot distribution of enrichment scores, optionally grouped.
    """
    fig, ax = plt.subplots(figsize=(5, 4))
    
    if group_col and group_col in df.columns:
        groups = df[group_col].unique()
        colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3']
        
        for i, group in enumerate(groups):
            mask = df[group_col] == group
            scores = df.loc[mask, score_col].values
            color = colors[i % len(colors)]
            ax.hist(scores, bins=30, density=True, alpha=0.5, 
                   label=f"{group} (n={mask.sum()}, μ={np.mean(scores):.3f})",
                   color=color)
    else:
        scores = df[score_col].values
        ax.hist(scores, bins=50, density=True, alpha=0.7, color='#1f77b4')
    
    ax.axvline(0, color='gray', linestyle='--', linewidth=1)
    ax.set_xlabel('Enrichment Score')
    ax.set_ylabel('Density')
    ax.set_title(title)
    
    if group_col:
        ax.legend(loc='upper right', fontsize=8)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
def plot_top_genes_bar(
    df: pd.DataFrame,
    n_top: int = 20,
    score_col: str = "mimetic_score",
    gene_col: str = "target_gene",
    title: str = "Top Genes",
    output_path: Path = None,
    color: str = None
):
    """Plot bar chart of top N genes by score."""
    # Get top and bottom genes
    df_sorted = df.sort_values(score_col, ascending=False)
    
    top_genes = df_sorted.head(n_top)
    bottom_genes = df_sorted.tail(n_top)
    
    fig, axes = plt.subplots(1, 2, figsize=(8, 5))
    
    # Top mimetics
    ax1 = axes[0]
    colors_top = ['#4DAF4A'] * n_top  # Green for positive
    ax1.barh(range(n_top), top_genes[score_col].values[::-1], color=colors_top, edgecolor='white')
    ax1.set_yticks(range(n_top))
    ax1.set_yticklabels(top_genes[gene_col].values[::-1])
    ax1.set_xlabel('Mimetic Score')
    ax1.set_title('Top Mimetics')
    
    # Top antagonists
    ax2 = axes[1]
    colors_bottom = ['#E41A1C'] * n_top  # Red for negative
    ax2.barh(range(n_top), bottom_genes[score_col].values[::-1], color=colors_bottom, edgecolor='white')
    ax2.set_yticks(range(n_top))
    ax2.set_yticklabels(bottom_genes[gene_col].values[::-1])
    ax2.set_xlabel('Mimetic Score')
    ax2.set_title('Top Antagonists')
    
    plt.suptitle(title, fontsize=12)
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, axes

#%%
def plot_overlap_bar(
    overlaps: dict,
    title: str = "Gene Set Overlaps",
    output_path: Path = None
):
    """Plot bar chart showing overlap counts between sets."""
    fig, ax = plt.subplots(figsize=(6, 4))
    
    names = list(overlaps.keys())
    counts = list(overlaps.values())
    
    colors = ['#1f77b4' if c > 0 else '#CCCCCC' for c in counts]
    bars = ax.bar(range(len(names)), counts, color=colors, edgecolor='white')
    
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Overlap Count')
    ax.set_title(title)
    
    # Add count labels on bars
    for bar, count in zip(bars, counts):
        if count > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                   str(count), ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig, ax

#%%
if __name__ == "__main__":
    print("Testing visualization module...")
    
    # Create test data
    np.random.seed(42)
    test_df = pd.DataFrame({
        "target_gene": [f"GENE{i}" for i in range(100)],
        "mimetic_score": np.random.randn(100) * 0.1,
        "x_score": np.random.randn(100) * 0.2,
        "y_score": np.random.randn(100) * 0.2,
    })
    
    # Test distribution plot
    output_dir = OUTPUT_DIRS["cache"]
    fig, ax = plot_mimetic_distribution(test_df, output_path=output_dir / "test_distribution")
    
    # Test scatter plot
    sig_mask = np.abs(test_df["x_score"]) > 0.3
    fig, ax = plot_convergence_scatter(
        test_df, "x_score", "y_score",
        significant_mask=sig_mask,
        output_path=output_dir / "test_scatter"
    )
    
    print("Visualization tests passed!")
