# %%
import glob
import math
import os
import string

import matplotlib.font_manager as fm
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

# Configuration
cm = 1/2.54
path = "./Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 6
plot_linewidth = 0.5

figsize = (10*cm, 10*cm)
fig = plt.figure(figsize=figsize)

pad = 11
# Define GridSpec layout
row = 99+ pad*2
col = 100
gsfig = gridspec.GridSpec(
    row, col,
    left=0,
    right=1,
    bottom=0,
    top=1,
    wspace=1,
    hspace=1
)

# Create subplots within GridSpec
gs1 = gsfig[0:33, 0:50]
gs2 = gsfig[33+pad:66+pad, 0:50]
gs3 = gsfig[66+pad*2:99+pad*2, 0:50]
gs4 = gsfig[0:row, 70:100]
ax1 = fig.add_subplot(gs1)
ax2 = fig.add_subplot(gs2)
ax3 = fig.add_subplot(gs3)
ax4 = fig.add_subplot(gs4)

axes = [ax1, ax2, ax3]

# Define data directories
DEG = "./DEG"
META = "./Metadata"
TARGET = "./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below"
df_reg_res = pd.read_csv(os.path.join(TARGET, "regression_results.txt"), sep="\t")
list_genesym_target = df_reg_res["Selected_Features"].to_list()

list_deg_files = sorted(glob.glob(f"{DEG}/**/*healthy.txt", recursive=True))
list_deg_samples = sorted(glob.glob(f"{META}/**/*long.txt", recursive=True))
list_deg_title = ["Acute Phase", "Mid Phase", "Late Phase"]

# Colors for each category
colors = {
    "COVID19": "tab:blue",
    "AgePred": "tab:orange",
    "None": "tab:grey"
}

linecolors = {
    "COVID19": "dodgerblue",
    "AgePred": "gold",
    "None": "silver"
}

list_df_deg_phase_target = list()
for i, (file, sample, title) in enumerate(zip(list_deg_files, list_deg_samples, list_deg_title)):
    fig_letter = list(string.ascii_uppercase)[i]
    df_deg = pd.read_csv(file, sep='\t')
    df_deg_phase = df_deg.copy()
    df_deg_phase["phase"] = title
    df_deg_phase["genesym"] = df_deg_phase["ID"].apply(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]))
    df_deg_phase_target = df_deg_phase[df_deg_phase["genesym"].isin(list_genesym_target)]
    list_df_deg_phase_target.append(df_deg_phase_target)

df_deg_phase_target = pd.concat(list_df_deg_phase_target, axis=0)
df_deg_phase_target.to_csv("./Supplementary_Tables/Table_S5.txt", sep="\t", index=False)

list_fc_sig_visits = []
list_fc_nosig_visits = []
list_fc_target_visits = []
for i, (file, sample, title) in enumerate(zip(list_deg_files, list_deg_samples, list_deg_title)):
    fig_letter = list(string.ascii_uppercase)[i]
    df_deg = pd.read_csv(file, sep='\t')
    df_deg_phase = df_deg.copy()
    df_deg_phase = df_deg_phase[df_deg_phase["baseMean"] > 1]
    df_deg_phase = df_deg_phase[df_deg_phase["log2FoldChange"].notna()]
    df_deg_phase = df_deg_phase[df_deg_phase["padj"].notna()]
    df_deg_phase["significance"] = df_deg_phase["padj"].apply(lambda x: -math.log10(float(x) + 1e-300))
    df_deg_phase["magnitude"] = df_deg_phase["log2FoldChange"].apply(abs)

    df_deg_phase_sig_only = df_deg_phase[np.logical_and(df_deg_phase["padj"] < 0.05, df_deg_phase["magnitude"] >= 1)]
    df_deg_phase_no_sig_only = df_deg_phase[np.logical_and(df_deg_phase["padj"] >= 0.05, df_deg_phase["magnitude"] < 1)]

    list_genes = df_deg_phase["ID"].to_list()
    list_fc = df_deg_phase["magnitude"].to_list()
    list_padj = df_deg_phase["significance"].to_list()

    list_genes_sig = df_deg_phase_sig_only["ID"].to_list()
    list_genesym_sig = list(map(lambda x: x.split(".")[0] + "_" + "_".join(str(x).split("_")[1:]), list_genes_sig))
    list_fc_sig = df_deg_phase_sig_only["magnitude"].to_list()
    list_fc_sig_visits.append(list_fc_sig)
    list_padj_sig = df_deg_phase_sig_only["significance"].to_list()

    list_genes_nosig = df_deg_phase_no_sig_only["ID"].to_list()
    list_genesym_nosig = list(map(lambda x: x.split(".")[0] + "_" + "_".join(str(x).split("_")[1:]), list_genes_nosig))
    list_fc_nosig = df_deg_phase_no_sig_only["magnitude"].to_list()
    list_fc_nosig_visits.append(list_fc_nosig)
    list_padj_nosig = df_deg_phase_no_sig_only["significance"].to_list()

    def get_dict_fc_padj(list_genes, list_fc, list_padj, list_genesym):
        dict_gene_fc_padj = {}
        for gene, fc, padj in zip(list_genes, list_fc, list_padj):
            gene = str(gene).split(".")[0] + "_" + "_".join(str(gene).split("_")[1:])
            fc = float(fc)
            padj = float(padj)
            if gene in list_genesym:
                dict_gene_fc_padj[gene] = {"fc": fc, "padj": padj}
        
        return dict_gene_fc_padj

    dict_gene_fc_padj = get_dict_fc_padj(list_genes, list_fc, list_padj, list_genesym_target)
    list_fc_target = list(map(lambda x: dict_gene_fc_padj[x]["fc"] if x in dict_gene_fc_padj else 0, list_genesym_target))
    list_padj_target = list(map(lambda x: dict_gene_fc_padj[x]["padj"] if x in dict_gene_fc_padj else 1, list_genesym_target))
    list_genesymbol = list(map(lambda x: x.split("_")[-1], list_genesym_target))
    list_fc_target_visits.append(list_fc_target)

    axes[i].scatter(x=list_fc, y=list_padj, s=2, edgecolors=colors["None"], linewidths=0.3, facecolor="None", alpha=1, zorder=0.5, label="None")
    axes[i].scatter(x=list_fc_sig, y=list_padj_sig, s=2, edgecolors=colors["COVID19"], linewidths=0.3, facecolor="None", alpha=1, zorder=0.5, label="COVID19")
    axes[i].scatter(x=list_fc_target, y=list_padj_target, s=3, edgecolors=colors["AgePred"], linewidths=0.6, facecolor="None", zorder=900, label="AgePred")

    for fc, padj, genesymbol in zip(list_fc_target, list_padj_target, list_genesymbol):
        if genesymbol == "VSIG4":
            axes[i].plot([fc, 4.0], [padj, 100], linewidth=0.4, color="k", zorder=999)
            axes[i].hlines(y=100, xmin=4.0, xmax=4.9, linewidth=0.4, color="k", zorder=999)
            axes[i].text(5, 100, genesymbol, fontstyle="italic", color=colors["AgePred"], ha="left", va="center", weight="bold", fontsize=plt.rcParams["font.size"], zorder=999)

    df_meta = pd.read_csv(sample, sep='\t')
    axes[i].set_title(f"{title} (N={len(df_meta)})", fontsize=plt.rcParams["font.size"] + 1, pad=0.2, weight="bold")
    axes[i].set_ylabel(f"Significance\n{r'-log$_{10}$(FDR)'}", fontsize=plt.rcParams["font.size"] + 1)
    if i == 2:
        axes[i].set_xlabel(f"Absolute Relative Expression\n|{r'log$_{2}$(COVID19/Healthy)'}|", fontsize=plt.rcParams["font.size"] + 1)
    axes[i].axvline(x=(np.mean(list_fc_target)), linewidth=plot_linewidth, linestyle="solid", color=linecolors["AgePred"], zorder=10)
    axes[i].axhline(y=(np.mean(list_padj_target)), linewidth=plot_linewidth, linestyle="solid", color=linecolors["AgePred"], zorder=10)
    axes[i].axvline(x=(np.mean(list_fc_sig)), linewidth=plot_linewidth, linestyle="solid", color=linecolors["COVID19"], zorder=20)
    axes[i].axhline(y=(np.mean(list_padj_sig)), linewidth=plot_linewidth, linestyle="solid", color=linecolors["COVID19"], zorder=20)
    axes[i].axvline(x=(np.mean(list_fc_nosig)), linewidth=plot_linewidth, linestyle="solid", color=linecolors["None"], zorder=5)
    axes[i].axhline(y=(np.mean(list_padj_nosig)), linewidth=plot_linewidth, linestyle="solid", color=linecolors["None"], zorder=5)
    axes[i].text(-0.33, 1.0, fig_letter, weight="bold", fontsize=plt.rcParams["font.size"] + 2, zorder=999, transform=axes[i].transAxes)
    axes[i].legend(markerscale=2.5, edgecolor="k", frameon=False)
    axes[i].set_xlim(-1.0, 8.0)
    axes[i].set_ylim(-10, 310)

fig_letter = list(string.ascii_uppercase)[len(list_deg_files)]

# Provided data
data = {
    "Acute": {
        "COVID19": list_fc_sig_visits[0],
        "AgePred": list_fc_target_visits[0],
        "None": list_fc_nosig_visits[0]
    },
    "Mid": {
        "COVID19": list_fc_sig_visits[1],
        "AgePred": list_fc_target_visits[1],
        "None": list_fc_nosig_visits[1]
    },
    "Late": {
        "COVID19": list_fc_sig_visits[2],
        "AgePred": list_fc_target_visits[2],
        "None": list_fc_nosig_visits[2]
    }
}

# Position adjustment for boxplots
positions = np.array([0, 1, 2])  # Adjust the x-positions of boxplots
width = 0.25  # Width of each boxplot

from scipy.stats.mstats import kruskal

_, pval_acute = kruskal(list_fc_sig_visits[0], list_fc_target_visits[0], list_fc_nosig_visits[0])
_, pval_mid = kruskal(list_fc_sig_visits[1], list_fc_target_visits[1], list_fc_nosig_visits[1])
_, pval_late = kruskal(list_fc_sig_visits[2], list_fc_target_visits[2], list_fc_nosig_visits[2])

from scikit_posthocs import posthoc_dunn

df_acute = posthoc_dunn([list_fc_sig_visits[0], list_fc_target_visits[0], list_fc_nosig_visits[0]], p_adjust="bonferroni")
df_acute.columns = list(colors.keys())
df_acute.index = list(colors.keys())
df_mid = posthoc_dunn([list_fc_sig_visits[1], list_fc_target_visits[1], list_fc_nosig_visits[1]], p_adjust="bonferroni")
df_mid.columns = list(colors.keys())
df_mid.index = list(colors.keys())
df_late = posthoc_dunn([list_fc_sig_visits[2], list_fc_target_visits[2], list_fc_nosig_visits[2]], p_adjust="bonferroni")
df_late.columns = list(colors.keys())
df_late.index = list(colors.keys())

mean_values = {"COVID19": [], "AgePred": [], "None": []}
categories = list(mean_values.keys())
# Plotting the boxplots
for i, (phase, phase_data) in enumerate(data.items()):
    for j, (category, values) in enumerate(phase_data.items()):
        x = positions[i] + (j - 1) * width
        bp = ax4.boxplot([values], 
                              positions=[x], 
                              widths=width, 
                              patch_artist=True, 
                              boxprops=dict(facecolor=colors[category], 
                                            linewidth=0.2, 
                                            edgecolor='k'),
                              whiskerprops=dict(linewidth=0.2),
                              capprops=dict(linewidth=0.2),
                              medianprops=dict(color="k", linewidth=0.2),
                              meanprops=dict(marker="o", markersize=4, markeredgecolor="k", markeredgewidth=0.2, markerfacecolor=linecolors[category]),
                              showmeans=True,
                              showfliers=False)
        mean_values[category].append(np.mean(values))

        ax4.text(x=(0+0+width)/2, 
                      y=3.34, 
                      s=f"{df_acute.loc['AgePred', 'None']:.1e}", 
                      color= "k", 
                      fontsize=plt.rcParams["font.size"],
                      ha="center")
        ax4.hlines(y=3.3, xmin=0, xmax=0+width, color="k", linewidth=0.4)
        ax4.vlines(x=0, ymin=3.05, ymax=3.3, color="k", linewidth=0.4)
        ax4.vlines(x=0+width, ymin=3.05, ymax=3.3, color="k", linewidth=0.4)
    
        ax4.text(x=(1+1+width)/2, 
                      y=3.24, 
                      s=f"{df_mid.loc['AgePred', 'None']:.1e}", 
                      color= "k",
                      fontsize=plt.rcParams["font.size"],
                      ha="center")
        ax4.hlines(y=3.2, xmin=1, xmax=1+width, color="k", linewidth=0.4)
        ax4.vlines(x=1, ymin=2.95, ymax=3.2, color="k", linewidth=0.4)
        ax4.vlines(x=1+width, ymin=2.95, ymax=3.2, color="k", linewidth=0.4)
       
        ax4.text(x=(2+2+width)/2, 
                      y=3.14,
                      s=round(df_late.loc['AgePred', 'None'], 2),  
                      color= "red",
                      fontsize=plt.rcParams["font.size"],
                      ha="center")
        ax4.hlines(y=3.1, xmin=2, xmax=2+width, color="red", linewidth=0.4)
        ax4.vlines(x=2, ymin=2.85, ymax=3.1, color="red", linewidth=0.4)
        ax4.vlines(x=2+width, ymin=2.85, ymax=3.1, color="red", linewidth=0.4)

for category in mean_values.keys():
    adjusted_positions = positions + (categories.index(category) - 1) * width
    ax4.plot(adjusted_positions, mean_values[category], color=linecolors[category], linestyle='-', linewidth=2.0, zorder=2)

# Set x-tick labels and adjust plot properties
ax4.set_xticks(positions)
ax4.set_xticklabels(data.keys())
ax4.set_xlabel("Infection Phases", fontsize=plt.rcParams["font.size"]+1)
ax4.set_ylabel(f"Absolute Relative Expression", fontsize=plt.rcParams["font.size"]+1)
ax4.set_ylim(-0.1, 3.5)
# Add text and adjust layout
ax4.text(-0.43, 1.0, 
        fig_letter, 
        weight="bold", 
        fontsize=plt.rcParams["font.size"]+2, 
        zorder=999, 
        transform=ax4.transAxes)
ax4.spines[['right', 'top']].set_visible(False)

plt.tight_layout()
plt.savefig("./Figures/Supplementary_Figure_8.tiff", dpi=600, bbox_inches="tight")
plt.savefig("./Figures/Supplementary_Figure_8.png", dpi=600, bbox_inches="tight")
plt.show()
plt.close()
# %%
