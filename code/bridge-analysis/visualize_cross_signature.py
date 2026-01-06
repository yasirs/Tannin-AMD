import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Paths
BASE_DIR = "/home/ysuhail/work/Tannin-AMD"
DATA_DIR = os.path.join(BASE_DIR, "results/bridge-results")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/bridge-results")

sns.set_theme(style="white", context="talk")

def plot_cross_signature_dotplot():
    # Load the three result files
    df_h2o2 = pd.read_csv(os.path.join(DATA_DIR, "H2O2_Stress_bridge_correlations.csv"))
    df_prg4 = pd.read_csv(os.path.join(DATA_DIR, "PRG4_Baseline_bridge_correlations.csv"))
    df_rescue = pd.read_csv(os.path.join(DATA_DIR, "PRG4_Rescue_bridge_correlations.csv"))

    # Rename columns for merging
    df_h2o2 = df_h2o2.rename(columns={"spearman_rho": "H2O2_Stress", "pvalue": "H2O2_p"})
    df_prg4 = df_prg4.rename(columns={"spearman_rho": "PRG4_Baseline", "pvalue": "PRG4_p"})
    df_rescue = df_rescue.rename(columns={"spearman_rho": "PRG4_Rescue", "pvalue": "Rescue_p"})

    # Merge
    merged = df_h2o2.merge(df_prg4, on=["perturbation", "target_gene"])
    merged = merged.merge(df_rescue, on=["perturbation", "target_gene"])

    # Select top 30 interesting genes (top in any of the categories)
    top_stress = merged.sort_values("H2O2_Stress", ascending=False).head(10)
    top_baseline = merged.sort_values("PRG4_Baseline", ascending=False).head(10)
    top_rescue = merged.sort_values("PRG4_Rescue", ascending=False).head(10)
    
    selected_genes = pd.concat([top_stress, top_baseline, top_rescue]).drop_duplicates()
    
    # Reshape for plotting
    plot_data = selected_genes.melt(
        id_vars=["target_gene"], 
        value_vars=["H2O2_Stress", "PRG4_Baseline", "PRG4_Rescue"],
        var_name="Signature", value_name="Spearman_Rho"
    )

    plt.figure(figsize=(10, 12))
    
    # Create the dot plot
    # Color by Rho, size by Rho absolute value
    plot_data["Abs_Rho"] = plot_data["Spearman_Rho"].abs()
    
    g = sns.scatterplot(
        data=plot_data,
        x="Signature",
        y="target_gene",
        size="Abs_Rho",
        hue="Spearman_Rho",
        palette="RdBu_r",
        sizes=(20, 300),
        hue_norm=(-0.4, 0.4),
        edgecolor="black",
        linewidth=0.5
    )
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Spearman ρ")
    plt.title("Cross-Signature Comparison of Top Hits", pad=20)
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Cross_Signature_Comparison_DotPlot.png"), dpi=300)
    plt.close()

if __name__ == "__main__":
    print("Generating Cross-Signature Dot Plot...")
    plot_cross_signature_dotplot()
    print(f"Comparison plot saved to {OUTPUT_DIR}")
