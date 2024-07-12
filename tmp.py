# %%
import glob
import math
import os
import string

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

# Configuration
cm = 1/2.54
path = "/BiO/Access/kyungwhan1998/miniconda3/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
plot_linewidth = 0.5

figsize = (10*cm, 5*cm)
fig = plt.figure(figsize=figsize)

# Define GridSpec layout
col = 180
row = 210
gsfig = gridspec.GridSpec(
    row, col,
    left=0,
    right=1,
    bottom=0,
    top=1,
    wspace=1,
    hspace=1
)

# Define data directories
DEG = "./DEG"
META = "./Metadata"
TARGET = "./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below"
df_reg_res = pd.read_csv(os.path.join(TARGET, "regression_results.txt"), sep="\t")
list_deg_files = sorted(glob.glob(f"{DEG}/**/*healthy.txt", recursive=True))
list_deg_samples = sorted(glob.glob(f"{META}/**/*long.txt", recursive=True))
list_deg_title = ["Acute Phase", "Mid Phase", "Late Phase"]

# Colors for each category
colors = {
    "COVID19": "tab:green",
    "AgePred": "tab:red",
    "None": "silver"
}

list_fc_sig_visits = []
list_fc_nosig_visits = []
list_fc_target_visits = []

for i, (file, sample, title) in enumerate(zip(list_deg_files, list_deg_samples, list_deg_title)):
    fig_letter = list(string.ascii_lowercase)[i]
    df_deg = pd.read_csv(file, sep='\t')
    df_deg_phase = df_deg.copy()
    df_deg_phase = df_deg_phase[df_deg_phase["baseMean"] > 1]
    df_deg_phase = df_deg_phase[df_deg_phase["log2FoldChange"].notna()]
    df_deg_phase = df_deg_phase[df_deg_phase["padj"].notna()]
    df_deg_phase["significance"] = df_deg_phase["padj"].apply(lambda x: -math.log10(float(x) + 1e-300))
    df_deg_phase["magnitude"] = df_deg_phase["log2FoldChange"].apply(abs)

    df_deg_phase_sig_only = df_deg_phase[np.logical_and(df_deg_phase["padj"] < 0.05, df_deg_phase["magnitude"] >= 1)]
    df_deg_phase_no_sig_only = df_deg_phase[np.logical_and(df_deg_phase["padj"] >= 0.05, df_deg_phase["magnitude"] < 1)]

    list_genesym_target = df_reg_res["Selected_Features"].to_list()
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

    ax_main = fig.add_subplot(gsfig[60 * i:60 * (i + 1), 0:100])  # Main subplot
    ax_inset = ax_main.inset_axes([0.6, 0.6, 0.35, 0.35])  # Inset axis

    # Main scatter plot
    ax_main.scatter(x=list_fc, y=list_padj, s=2, edgecolors="darkgrey", linewidths=0.3, facecolor="None", alpha=1, zorder=0.5, label="None")
    ax_main.scatter(x=list_fc_sig, y=list_padj_sig, s=2, edgecolors="forestgreen", linewidths=0.3, facecolor="None", alpha=1, zorder=0.5, label="COVID19")
    ax_main.scatter(x=list_fc_target, y=list_padj_target, s=4, edgecolors="firebrick", linewidths=0.3, facecolor="None", zorder=900, label="AgePred")

    for fc, padj, genesymbol in zip(list_fc_target, list_padj_target, list_genesymbol):
        if genesymbol == "VSIG4":
            ax_main.plot([fc, 4.0], [padj, 100], linewidth=0.4, color="k", zorder=999)
            ax_main.hlines(y=100, xmin=4.0, xmax=4.9, linewidth=0.4, color="k", zorder=999)
            ax_main.text(5, 100, genesymbol, fontstyle="italic", color="firebrick", ha="left", va="center", weight="bold", fontsize=plt.rcParams["font.size"], zorder=999)

    df_meta = pd.read_csv(sample, sep='\t')
    ax_main.set_title(f"{title} (N={len(df_meta)})", fontsize=plt.rcParams["font.size"] + 1, pad=0.2, weight="bold")
    ax_main.set_xlabel(f"Absolute Relative Expression\n|{r'log$_{2}$(COVID19/Healthy)'}|", fontsize=plt.rcParams["font.size"] + 1)
    ax_main.set_ylabel(f"Significance\n{r'-log$_{10}$(FDR)'}", fontsize=plt.rcParams["font.size"] + 1)
    ax_main.axvline(x=(np.mean(list_fc_target)), linewidth=plot_linewidth, linestyle="solid", color=colors["AgePred"], zorder=10)
    ax_main.axhline(y=(np.mean(list_padj_target)), linewidth=plot_linewidth, linestyle="solid", color=colors["AgePred"], zorder=10)
    ax_main.axvline(x=(np.mean(list_fc_sig)), linewidth=plot_linewidth, linestyle="solid", color=colors["COVID19"], zorder=20)
    ax_main.axhline(y=(np.mean(list_padj_sig)), linewidth=plot_linewidth, linestyle="solid", color=colors["COVID19"], zorder=20)
    ax_main.axvline(x=(np.mean(list_fc_nosig)), linewidth=plot_linewidth, linestyle="solid", color=colors["None"], zorder=5)
    ax_main.axhline(y=(np.mean(list_padj_nosig)), linewidth=plot_linewidth, linestyle="solid", color=colors["None"], zorder=5)
    ax_main.text(-0.5, 1.05, fig_letter, weight="bold", fontsize=plt.rcParams["font.size"] + 2, zorder=999, transform=ax_main.transAxes)
    ax_main.legend(markerscale=2.5, edgecolor="k", frameon=False)
    ax_main.set_xlim(-1.0, 8.0)
    ax_main.set_ylim(-10, 310)

    # Inset scatter plot (same data for demonstration)
    ax_inset.scatter(x=list_fc, y=list_padj, s=2, edgecolors="darkgrey", linewidths=0.3, facecolor="None", alpha=1, zorder=0.5, label="None")
    ax_inset.scatter(x=list_fc_sig, y=list_padj_sig, s=2, edgecolors="forestgreen", linewidths=0.3, facecolor="None", alpha=1, zorder=0.5, label="COVID19")
    ax_inset.scatter(x=list_fc_target, y=list_padj_target, s=4, edgecolors="firebrick", linewidths=0.3, facecolor="None", zorder=900, label="AgePred")
    ax_inset.set_xlim(-1.0, 8.0)
    ax_inset.set_ylim(-10, 310)
    ax_inset.set_title("Inset Plot", fontsize=plt.rcParams["font.size"], pad=0.2, weight="bold")

plt.show()