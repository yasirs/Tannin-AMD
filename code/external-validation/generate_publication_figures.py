#!/usr/bin/env python3
"""
Publication-Ready Figures for Age-NRF2-ARE Analysis
Generates both static (PDF/TIFF) and interactive (HTML) versions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from pathlib import Path
from PIL import Image

# Try plotly for Sankey
try:
    import plotly.graph_objects as go
    import kaleido
    HAS_PLOTLY = True
except:
    HAS_PLOTLY = False
    print("Warning: plotly/kaleido not available, skipping Sankey diagrams")

# Paths
RESULTS_DIR = Path("/home/ysuhail/work/Tannin-AMD/results/cohort-GSE29801/age_nrf2_analysis")

# Load data
print("Loading data...")
forward_df = pd.read_csv(RESULTS_DIR / "forward_age_de_are_gsea.csv")
aging_mim_df = pd.read_csv(RESULTS_DIR / "backward_aging_mimetics.csv")
prg4_mim_df = pd.read_csv(RESULTS_DIR / "backward_prg4_mimetics.csv")
convergence = pd.read_csv(RESULTS_DIR / "convergence_overlaps.csv")
are_genes = pd.read_csv(RESULTS_DIR / "ARE_gene_set.csv")['gene'].tolist()

# Define convergent genes
aging_antagonists = set(aging_mim_df.sort_values('mimetic_score').head(200)['target_gene'])
prg4_mimetics = set(prg4_mim_df.sort_values('prg4_mimetic_score', ascending=False).head(200)['target_gene'])
convergent_genes = aging_antagonists & prg4_mimetics
mirna_genes = {'DICER1', 'DROSHA', 'XPO5'}

# Helper function to save both PDF and TIFF
def save_figure(fig, name, dpi=300):
    """Save figure in both PDF and LZW-compressed TIFF"""
    pdf_path = RESULTS_DIR / f"{name}.pdf"
    tiff_path = RESULTS_DIR / f"{name}.tiff"
    
    fig.savefig(pdf_path, dpi=dpi, bbox_inches='tight')
    
    # Save as PNG first, then convert to TIFF with compression
    png_temp = RESULTS_DIR / f"{name}_temp.png"
    fig.savefig(png_temp, dpi=dpi, bbox_inches='tight')
    
    # Convert to TIFF with LZW compression
    img = Image.open(png_temp)
    img.save(tiff_path, compression='tiff_lzw')
    png_temp.unlink()  # Remove temp file
    
    print(f"  Saved: {name}.pdf, {name}.tiff")

#%% ==========================================================================
# FIGURE 1: SCHEMATIC - Analysis Overview
# ==========================================================================
print("\n[Figure 1] Creating analysis overview schematic...")

fig, ax = plt.subplots(figsize=(12, 8))
ax.set_xlim(0, 12)
ax.set_ylim(0, 8)
ax.axis('off')

# Title
ax.text(6, 7.5, 'Age-NRF2-ARE Analysis Overview', fontsize=16, fontweight='bold', 
        ha='center', va='center')

# Box style
box_style = dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', linewidth=2)

# Data sources (left)
ax.add_patch(FancyBboxPatch((0.5, 5.5), 2.5, 1.2, boxstyle='round,pad=0.1', 
                             facecolor='#E8F4F8', edgecolor='#2C3E50', linewidth=2))
ax.text(1.75, 6.1, 'GSE29801\nAge-DE Genes', fontsize=10, ha='center', va='center', fontweight='bold')
ax.text(1.75, 5.7, 'UP: 322  DOWN: 150', fontsize=8, ha='center', va='center')

ax.add_patch(FancyBboxPatch((0.5, 3.8), 2.5, 1.2, boxstyle='round,pad=0.1', 
                             facecolor='#FEF9E7', edgecolor='#B7950B', linewidth=2))
ax.text(1.75, 4.4, 'PRG4 Rescue\nSignature', fontsize=10, ha='center', va='center', fontweight='bold')
ax.text(1.75, 4.0, 'UP: 4012  DOWN: 2257', fontsize=8, ha='center', va='center')

ax.add_patch(FancyBboxPatch((0.5, 2.1), 2.5, 1.2, boxstyle='round,pad=0.1', 
                             facecolor='#FDEDEC', edgecolor='#C0392B', linewidth=2))
ax.text(1.75, 2.7, 'ARE/NRF2\nGene Set', fontsize=10, ha='center', va='center', fontweight='bold')
ax.text(1.75, 2.3, '242 genes', fontsize=8, ha='center', va='center')

# Perturb-seq (center)
ax.add_patch(FancyBboxPatch((4, 3.5), 3.5, 2.5, boxstyle='round,pad=0.1', 
                             facecolor='#E8DAEF', edgecolor='#8E44AD', linewidth=2))
ax.text(5.75, 5.5, 'K562 GWPS Perturb-seq', fontsize=11, ha='center', va='center', fontweight='bold')
ax.text(5.75, 4.9, '11,258 knockdowns', fontsize=9, ha='center', va='center')
ax.text(5.75, 4.4, '8,248 genes profiled', fontsize=9, ha='center', va='center')
ax.text(5.75, 3.8, 'Compute mimetic scores\nusing GSEA-based methods', fontsize=8, ha='center', va='center', style='italic')

# Results (right)
ax.add_patch(FancyBboxPatch((8.5, 5.0), 3, 1.5, boxstyle='round,pad=0.1', 
                             facecolor='#D5F5E3', edgecolor='#1E8449', linewidth=2))
ax.text(10, 5.8, 'Aging Antagonists', fontsize=10, ha='center', va='center', fontweight='bold')
ax.text(10, 5.3, 'Top: RBM25, FGFR1OP\nCNOT3, MED9', fontsize=8, ha='center', va='center')

ax.add_patch(FancyBboxPatch((8.5, 3.0), 3, 1.5, boxstyle='round,pad=0.1', 
                             facecolor='#D4E6F1', edgecolor='#2874A6', linewidth=2))
ax.text(10, 3.8, 'PRG4 Mimetics', fontsize=10, ha='center', va='center', fontweight='bold')
ax.text(10, 3.3, 'Top: SHOC2, DICER1\nAMBRA1, IMPDH2', fontsize=8, ha='center', va='center')

# Convergence box (bottom center)
ax.add_patch(FancyBboxPatch((4, 0.5), 4, 1.5, boxstyle='round,pad=0.1', 
                             facecolor='#FCF3CF', edgecolor='#D68910', linewidth=3))
ax.text(6, 1.5, '14 Convergent Genes', fontsize=12, ha='center', va='center', fontweight='bold')
ax.text(6, 1.0, 'DICER1, DROSHA, XPO5 (miRNA)\nMED14, RBBP6, TANGO6, UBA3...', 
        fontsize=9, ha='center', va='center')

# Arrows
arrow_style = dict(arrowstyle='->', color='#34495E', lw=2, mutation_scale=15)
# Data to Perturb-seq
ax.annotate('', xy=(4, 4.75), xytext=(3, 6.1), arrowprops=arrow_style)
ax.annotate('', xy=(4, 4.75), xytext=(3, 4.4), arrowprops=arrow_style)
ax.annotate('', xy=(4, 4.75), xytext=(3, 2.7), arrowprops=arrow_style)
# Perturb-seq to results
ax.annotate('', xy=(8.5, 5.5), xytext=(7.5, 4.75), arrowprops=arrow_style)
ax.annotate('', xy=(8.5, 3.5), xytext=(7.5, 4.25), arrowprops=arrow_style)
# Results to convergence
ax.annotate('', xy=(7, 2), xytext=(9, 3), arrowprops=arrow_style)
ax.annotate('', xy=(5, 2), xytext=(9, 5), arrowprops=arrow_style)

plt.tight_layout()
save_figure(fig, "fig1_analysis_overview_schematic")
plt.close()

#%% ==========================================================================
# FIGURE 2: FORWARD ANALYSIS
# ==========================================================================
print("\n[Figure 2] Creating forward analysis figure...")

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

# Panel A: Histograms
for i, direction in enumerate(['UP', 'DOWN']):
    subset = forward_df[forward_df['age_direction'] == direction]
    color = '#E74C3C' if direction == 'UP' else '#3498DB'
    axes[0].hist(subset['ARE_enrichment'], bins=25, alpha=0.6, 
                 label=f'Age-{direction} (n={len(subset)})', color=color, edgecolor='black')
axes[0].axvline(0, color='black', linestyle='--', linewidth=1.5)
axes[0].set_xlabel('ARE Enrichment Score', fontsize=11)
axes[0].set_ylabel('Number of Knockdowns', fontsize=11)
axes[0].set_title('A. Distribution of ARE Effects', fontsize=12, fontweight='bold')
axes[0].legend(fontsize=9)

# Panel B: Top activators and suppressors (age-DOWN genes)
age_down = forward_df[forward_df['age_direction'] == 'DOWN'].copy()
top_activators = age_down.nlargest(8, 'ARE_enrichment')
top_suppressors = age_down.nsmallest(8, 'ARE_enrichment')
combined = pd.concat([top_activators, top_suppressors])
colors = ['#27AE60' if x > 0 else '#C0392B' for x in combined['ARE_enrichment']]

y_pos = range(len(combined))
axes[1].barh(y_pos, combined['ARE_enrichment'], color=colors, edgecolor='black')
axes[1].set_yticks(y_pos)
axes[1].set_yticklabels(combined['target_gene'], fontsize=9)
axes[1].axvline(0, color='black', linestyle='-', linewidth=1)
axes[1].set_xlabel('ARE Enrichment Score', fontsize=11)
axes[1].set_title('B. Top Age-DOWN Gene Effects on ARE', fontsize=12, fontweight='bold')
axes[1].invert_yaxis()

# Panel C: Summary statistics
summary_data = forward_df.groupby('age_direction')['ARE_enrichment'].agg(['mean', 'std', 'count'])
x_pos = [0, 1]
means = [summary_data.loc['UP', 'mean'], summary_data.loc['DOWN', 'mean']]
stds = [summary_data.loc['UP', 'std'], summary_data.loc['DOWN', 'std']]
colors = ['#E74C3C', '#3498DB']
bars = axes[2].bar(x_pos, means, yerr=stds, color=colors, edgecolor='black', capsize=5)
axes[2].set_xticks(x_pos)
axes[2].set_xticklabels(['Age-UP\n(n=153)', 'Age-DOWN\n(n=92)'], fontsize=11)
axes[2].axhline(0, color='black', linestyle='--', linewidth=1)
axes[2].set_ylabel('Mean ARE Enrichment ± SD', fontsize=11)
axes[2].set_title('C. Summary by Direction', fontsize=12, fontweight='bold')

plt.tight_layout()
save_figure(fig, "fig2_forward_analysis")
plt.close()

#%% ==========================================================================
# FIGURE 3: SANKEY DIAGRAM - Backward Analysis Flow
# ==========================================================================
if HAS_PLOTLY:
    print("\n[Figure 3] Creating Sankey diagram...")
    
    # Define nodes
    labels = [
        "All Knockdowns\n(11,258)",  # 0
        "Aging Mimetics\n(Top 200)",  # 1
        "Neutral",  # 2
        "Aging Antagonists\n(Top 200)",  # 3
        "PRG4 Mimetics\nOverlap (14)",  # 4
        "No PRG4 Overlap",  # 5
        "miRNA genes (3)",  # 6
        "Other (11)"  # 7
    ]
    
    colors = [
        "#34495E",  # All
        "#E74C3C",  # Mimetics
        "#95A5A6",  # Neutral
        "#27AE60",  # Antagonists  
        "#3498DB",  # PRG4 overlap
        "#BDC3C7",  # No PRG4
        "#9B59B6",  # miRNA
        "#F39C12"   # Other
    ]
    
    # Define flows (source, target, value)
    source = [0, 0, 0, 3, 3, 4, 4]
    target = [1, 2, 3, 4, 5, 6, 7]
    values = [200, 10858, 200, 14, 186, 3, 11]
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=20,
            thickness=30,
            line=dict(color="black", width=1),
            label=labels,
            color=colors
        ),
        link=dict(
            source=source,
            target=target,
            value=values,
            color=['rgba(231,76,60,0.4)', 'rgba(149,165,166,0.2)', 'rgba(39,174,96,0.4)',
                   'rgba(52,152,219,0.5)', 'rgba(189,195,199,0.3)', 
                   'rgba(155,89,182,0.5)', 'rgba(243,156,18,0.4)']
        )
    )])
    
    fig.update_layout(
        title_text="Backward Analysis: Flow from Knockdowns to Convergent Targets",
        title_font_size=16,
        font_size=12,
        width=900,
        height=500
    )
    
    # Save interactive
    fig.write_html(str(RESULTS_DIR / "fig3_sankey_interactive.html"))
    print("  Saved: fig3_sankey_interactive.html")
    
    # Try to save static (requires Chrome)
    try:
        fig.write_image(str(RESULTS_DIR / "fig3_sankey_flow.pdf"), scale=2)
        fig.write_image(str(RESULTS_DIR / "fig3_sankey_flow_temp.png"), scale=3)
        img = Image.open(RESULTS_DIR / "fig3_sankey_flow_temp.png")
        img.save(RESULTS_DIR / "fig3_sankey_flow.tiff", compression='tiff_lzw')
        (RESULTS_DIR / "fig3_sankey_flow_temp.png").unlink()
        print("  Saved: fig3_sankey_flow.pdf, fig3_sankey_flow.tiff")
    except Exception as e:
        print(f"  Warning: Could not save static Sankey (Chrome required): {e}")
        print("  Creating matplotlib flow diagram as fallback...")
        
        # Matplotlib fallback flow diagram
        fig_fb, ax_fb = plt.subplots(figsize=(12, 6))
        ax_fb.set_xlim(0, 12)
        ax_fb.set_ylim(0, 6)
        ax_fb.axis('off')
        
        # Title
        ax_fb.text(6, 5.7, 'Backward Analysis: Flow from Knockdowns to Therapeutic Targets', 
                  fontsize=14, fontweight='bold', ha='center')
        
        # Boxes
        boxes_data = [
            (0.5, 2.5, 2, 2, "All Knockdowns\n11,258", "#34495E", "white"),
            (3.5, 4, 2, 1.2, "Aging\nMimetics\n(200)", "#E74C3C", "white"),
            (3.5, 2.5, 2, 0.8, "Neutral", "#95A5A6", "black"),
            (3.5, 1, 2, 1.2, "Aging\nAntagonists\n(200)", "#27AE60", "white"),
            (7, 1.5, 2, 1.5, "PRG4 Mimetic\nOverlap\n(14 genes)", "#3498DB", "white"),
            (10, 1.5, 1.8, 1.5, "miRNA genes\n(3)\nDICER1\nDROSHA\nXPO5", "#9B59B6", "white"),
        ]
        
        for x, y, w, h, text, color, tcolor in boxes_data:
            box = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.1',
                                  facecolor=color, edgecolor='black', linewidth=2)
            ax_fb.add_patch(box)
            ax_fb.text(x + w/2, y + h/2, text, ha='center', va='center', 
                      fontsize=9, color=tcolor, fontweight='bold')
        
        # Arrows
        arrow_s = dict(arrowstyle='->', color='#2C3E50', lw=2, mutation_scale=15)
        ax_fb.annotate('', xy=(3.5, 4.6), xytext=(2.5, 3.5), arrowprops=arrow_s)
        ax_fb.annotate('', xy=(3.5, 2.9), xytext=(2.5, 3.5), arrowprops=arrow_s)
        ax_fb.annotate('', xy=(3.5, 1.6), xytext=(2.5, 3.5), arrowprops=arrow_s)
        ax_fb.annotate('', xy=(7, 2.25), xytext=(5.5, 1.6), arrowprops=arrow_s)
        ax_fb.annotate('', xy=(10, 2.25), xytext=(9, 2.25), arrowprops=arrow_s)
        
        plt.tight_layout()
        save_figure(fig_fb, "fig3_flow_diagram")
        plt.close(fig_fb)

#%% ==========================================================================
# FIGURE 4: CONVERGENCE SCATTER WITH HIGHLIGHTS
# ==========================================================================
print("\n[Figure 4] Creating convergence scatter plot...")

merged = aging_mim_df.merge(prg4_mim_df[['target_gene', 'prg4_mimetic_score']], on='target_gene')
merged = merged.dropna(subset=['mimetic_score', 'prg4_mimetic_score'])

fig, ax = plt.subplots(figsize=(10, 8))

# Background points
ax.scatter(-merged['mimetic_score'], merged['prg4_mimetic_score'], 
           alpha=0.2, s=10, c='#95A5A6', label='Other genes')

# Highlight ARE genes
are_mask = merged['target_gene'].isin(are_genes)
ax.scatter(-merged.loc[are_mask, 'mimetic_score'], 
           merged.loc[are_mask, 'prg4_mimetic_score'],
           alpha=0.7, s=40, c='#E74C3C', marker='s', label=f'ARE genes (n={are_mask.sum()})')

# Highlight convergent genes
conv_mask = merged['target_gene'].isin(convergent_genes)
conv_data = merged[conv_mask]
ax.scatter(-conv_data['mimetic_score'], conv_data['prg4_mimetic_score'],
           alpha=1, s=100, c='#3498DB', marker='D', edgecolors='black', linewidth=1.5,
           label=f'Convergent (n={len(conv_data)})')

# Highlight miRNA genes specially
mirna_mask = merged['target_gene'].isin(mirna_genes)
mirna_data = merged[mirna_mask]
ax.scatter(-mirna_data['mimetic_score'], mirna_data['prg4_mimetic_score'],
           alpha=1, s=200, c='#9B59B6', marker='*', edgecolors='black', linewidth=1.5,
           label='miRNA genes (DICER1, DROSHA, XPO5)', zorder=5)

# Label miRNA genes
for _, row in mirna_data.iterrows():
    ax.annotate(row['target_gene'], 
                (-row['mimetic_score'], row['prg4_mimetic_score']),
                xytext=(10, 10), textcoords='offset points',
                fontsize=10, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#9B59B6', lw=1.5))

# Reference lines
ax.axhline(0, color='black', linestyle=':', alpha=0.5)
ax.axvline(0, color='black', linestyle=':', alpha=0.5)

# Quadrant labels
ax.text(0.35, 0.18, 'Aging Antagonist\n& PRG4 Mimetic', fontsize=9, 
        ha='center', transform=ax.transAxes, color='#27AE60', fontweight='bold')
ax.text(0.15, 0.82, 'Aging Mimetic\n& PRG4 Mimetic', fontsize=9,
        ha='center', transform=ax.transAxes, color='#F39C12')

ax.set_xlabel('Aging Antagonist Score\n(reverses aging phenotype)', fontsize=12)
ax.set_ylabel('PRG4 Mimetic Score\n(mimics PRG4 effect)', fontsize=12)
ax.set_title('Convergence Analysis: Identifying Therapeutic Targets', fontsize=14, fontweight='bold')
ax.legend(loc='upper left', fontsize=10)

plt.tight_layout()
save_figure(fig, "fig4_convergence_scatter")
plt.close()

#%% ==========================================================================
# FIGURE 5: CONVERGENT GENES - miRNA PATHWAY SCHEMATIC
# ==========================================================================
print("\n[Figure 5] Creating miRNA pathway schematic...")

fig, ax = plt.subplots(figsize=(12, 7))
ax.set_xlim(0, 12)
ax.set_ylim(0, 7)
ax.axis('off')

# Title
ax.text(6, 6.7, 'Convergent Therapeutic Targets: miRNA Biogenesis Pathway', 
        fontsize=14, fontweight='bold', ha='center')

# Legend box
legend_box = FancyBboxPatch((0.2, 5.5), 3, 1, boxstyle='round,pad=0.1',
                             facecolor='white', edgecolor='black', linewidth=1)
ax.add_patch(legend_box)
ax.plot([0.4], [6.2], 'o', color='#27AE60', markersize=10)
ax.text(0.7, 6.2, 'Aging Antagonist', fontsize=9, va='center')
ax.plot([0.4], [5.8], 'o', color='#3498DB', markersize=10)
ax.text(0.7, 5.8, 'PRG4 Mimetic', fontsize=9, va='center')

# miRNA biosynthesis pathway boxes
boxes = [
    (1, 3.5, 'Pri-miRNA\n(nucleus)', '#F5EEF8'),
    (3.5, 3.5, 'DROSHA\n⬇️ KD reverses aging', '#D5F5E3'),
    (6, 3.5, 'Pre-miRNA', '#F5EEF8'),
    (8, 3.5, 'XPO5\n⬇️ KD reverses aging', '#D5F5E3'),
    (10, 3.5, 'Pre-miRNA\n(cytoplasm)', '#F5EEF8'),
]

for x, y, text, color in boxes:
    box = FancyBboxPatch((x-0.7, y-0.5), 1.6, 1.2, boxstyle='round,pad=0.1',
                          facecolor=color, edgecolor='black', linewidth=2)
    ax.add_patch(box)
    ax.text(x+0.1, y+0.1, text, ha='center', va='center', fontsize=9, fontweight='bold')

# DICER1 box (larger, emphasized)
dicer_box = FancyBboxPatch((4.5, 1.3), 3, 1.5, boxstyle='round,pad=0.1',
                            facecolor='#D4E6F1', edgecolor='#2874A6', linewidth=3)
ax.add_patch(dicer_box)
ax.text(6, 2.3, 'DICER1', ha='center', va='center', fontsize=14, fontweight='bold')
ax.text(6, 1.8, 'Knockdown:', ha='center', va='center', fontsize=10)
ax.text(6, 1.5, '✓ Reverses aging  ✓ Mimics PRG4', ha='center', va='center', 
        fontsize=10, color='#27AE60', fontweight='bold')

# Mature miRNA
ax.add_patch(FancyBboxPatch((8.5, 1.3), 2.5, 1.2, boxstyle='round,pad=0.1',
                             facecolor='#FCF3CF', edgecolor='#D68910', linewidth=2))
ax.text(9.75, 1.9, 'Mature miRNA', ha='center', va='center', fontsize=11, fontweight='bold')

# Arrows
arrow_style = dict(arrowstyle='->', color='#34495E', lw=2.5, mutation_scale=20)
ax.annotate('', xy=(2.7, 4.0), xytext=(1.8, 4.0), arrowprops=arrow_style)
ax.annotate('', xy=(5.2, 4.0), xytext=(4.4, 4.0), arrowprops=arrow_style)
ax.annotate('', xy=(7.2, 4.0), xytext=(6.8, 4.0), arrowprops=arrow_style)
ax.annotate('', xy=(9.2, 4.0), xytext=(8.8, 4.0), arrowprops=arrow_style)
ax.annotate('', xy=(10, 2.7), xytext=(10, 2.9), arrowprops=dict(arrowstyle='->', color='#34495E', lw=2, mutation_scale=15))
ax.annotate('', xy=(6, 2.8), xytext=(6, 2.9), arrowprops=dict(arrowstyle='->', color='#34495E', lw=2, mutation_scale=15))
ax.annotate('', xy=(8.5, 2.0), xytext=(7.5, 2.0), arrowprops=arrow_style)

# Interpretation text
ax.text(6, 0.5, 'Interpretation: Reducing miRNA biogenesis both reverses aging phenotype AND mimics PRG4 therapeutic effect',
        ha='center', va='center', fontsize=10, style='italic',
        bbox=dict(boxstyle='round', facecolor='#FDEBD0', edgecolor='#D68910'))

plt.tight_layout()
save_figure(fig, "fig5_mirna_pathway_schematic")
plt.close()

#%% ==========================================================================
# FIGURE 6: SUMMARY BAR CHART - Interpretable Overview
# ==========================================================================
print("\n[Figure 6] Creating summary bar chart...")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Set sizes and overlaps
labels = ['Aging\nAntagonists', 'PRG4\nMimetics', 'ARE\nGenes', 
          'Antag. ∩\nPRG4 Mim.', 'miRNA in\nConvergent']
values = [199, 194, 155, 14, 3]
colors = ['#27AE60', '#3498DB', '#E74C3C', '#9B59B6', '#F39C12']

bars = axes[0].bar(range(len(labels)), values, color=colors, edgecolor='black', linewidth=1.5)
axes[0].set_xticks(range(len(labels)))
axes[0].set_xticklabels(labels, fontsize=10)
axes[0].set_ylabel('Number of Genes', fontsize=12)
axes[0].set_title('A. Gene Set Sizes and Key Overlaps', fontsize=12, fontweight='bold')

for bar, val in zip(bars, values):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3, 
                 str(val), ha='center', va='bottom', fontsize=11, fontweight='bold')

# Panel B: Top convergent genes with scores
top_conv = []
for gene in convergent_genes:
    age_score = aging_mim_df[aging_mim_df['target_gene'] == gene]['mimetic_score'].values
    prg4_score = prg4_mim_df[prg4_mim_df['target_gene'] == gene]['prg4_mimetic_score'].values
    if len(age_score) > 0 and len(prg4_score) > 0:
        top_conv.append({
            'gene': gene,
            'age_antagonist': -age_score[0],
            'prg4_mimetic': prg4_score[0],
            'is_mirna': gene in mirna_genes
        })

top_conv_df = pd.DataFrame(top_conv)
top_conv_df = top_conv_df.sort_values('age_antagonist', ascending=True)

y_pos = range(len(top_conv_df))
colors = ['#9B59B6' if x else '#3498DB' for x in top_conv_df['is_mirna']]
bars = axes[1].barh(y_pos, top_conv_df['age_antagonist'], color=colors, edgecolor='black')
axes[1].set_yticks(y_pos)
axes[1].set_yticklabels(top_conv_df['gene'], fontsize=10)
axes[1].set_xlabel('Aging Antagonist Score', fontsize=12)
axes[1].set_title('B. Convergent Genes (sorted by antagonist strength)', fontsize=12, fontweight='bold')

# Legend
purple_patch = mpatches.Patch(color='#9B59B6', label='miRNA biogenesis')
blue_patch = mpatches.Patch(color='#3498DB', label='Other function')
axes[1].legend(handles=[purple_patch, blue_patch], loc='lower right', fontsize=9)

plt.tight_layout()
save_figure(fig, "fig6_summary_overview")
plt.close()

#%% ==========================================================================
# Create summary table for report
# ==========================================================================
print("\n[Table] Creating summary statistics table...")

summary_stats = pd.DataFrame({
    'Metric': [
        'Age-DE genes (extramacular)',
        'Age-UP genes',
        'Age-DOWN genes',
        'PRG4-UP genes (lfc > 0.5)',
        'PRG4-DOWN genes (lfc < -0.5)',
        'ARE/NRF2 genes',
        'Total Perturb-seq knockdowns',
        'Age-DE with knockdowns',
        'Top 200 Aging Antagonists',
        'Top 200 PRG4 Mimetics',
        'Convergent (Antagonist ∩ PRG4 Mim.)',
        'Convergent miRNA genes'
    ],
    'Value': [
        '472',
        '322',
        '150',
        '4,012',
        '2,257',
        '242',
        '11,258',
        '245',
        '199 unique',
        '194 unique',
        '14',
        '3 (DICER1, DROSHA, XPO5)'
    ]
})
summary_stats.to_csv(RESULTS_DIR / "summary_statistics_table.csv", index=False)
print("  Saved: summary_statistics_table.csv")

print("\n" + "="*60)
print("ALL FIGURES GENERATED SUCCESSFULLY")
print("="*60)

# List all output files
print("\nOutput files:")
for f in sorted(RESULTS_DIR.glob("fig*.pdf")):
    print(f"  {f.name}")
for f in sorted(RESULTS_DIR.glob("fig*.tiff")):
    print(f"  {f.name}")
for f in sorted(RESULTS_DIR.glob("*.html")):
    print(f"  {f.name}")
