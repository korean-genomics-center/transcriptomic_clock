# %%
import os
import json
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib.lines import Line2D
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from transcriptomic_clock import Clock
import draw_enrichment_plot as enrich

# %% 
WORKDIR = ".."
cm = 1/2.54
path = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 7
width = 16*cm
height = 16*cm
figsize = (width, height)
plot_linewidth = 3

fig = plt.figure(figsize = figsize)
col = 180
row = 210
gsfig = gridspec.GridSpec(
                        row, col,
                        left = 0, 
                        right = 1, 
                        bottom = 0,
                        top = 1,
                        wspace = 1, 
                        hspace = 1)

gs1 = gsfig[10:60, 0:45]
ax1 = fig.add_subplot(gs1)

gs2 = gsfig[80:130, 0:45]
ax2 = fig.add_subplot(gs2)

gs3 = gsfig[150:200, 0:45] 
ax3 = fig.add_subplot(gs3)

gs4 = gsfig[10:200, 55:65]
ax4 = fig.add_subplot(gs4) 

gs5 = gsfig[15:85, 155:180]
ax5 = fig.add_subplot(gs5)

gs6 = gsfig[110:160, 155:180]
ax6 = fig.add_subplot(gs6)

gs7 = gsfig[180:185, 110:160]
ax7 = fig.add_subplot(gs7)

gs8 = gsfig[190:200, 110:158]
ax8 = fig.add_subplot(gs8)

train_group = "healthy_cohort_train_prev_added_back_10_80_90s_median_filter"
list_test_group = [train_group, "healthy_cohort_valid_prev", "healthy_cohort2"]
dir_trained_model = f"{WORKDIR}/LASSO_INFO/{train_group}/corr_0.35_above_pval_0.05_below"
path_dict_name_conversion_json = f"{WORKDIR}/Data/dictionary_name_conversion.json"
with open (path_dict_name_conversion_json, mode="rb") as fb1:
    dict_name_conversion = json.load(fb1)

list_fig_letters = ["A", "B", "C"]
axes = [ax1, ax2, ax3]
for i, fig_letter in enumerate(list_fig_letters):
    file_testing_data = f"{dir_trained_model}/{list_test_group[i]}_std_scaled.txt"
    os.makedirs(dir_trained_model, exist_ok=True)
    clock = Clock(dir_trained_model, file_testing_data)
    df_stat = clock.make_dataframe_error_statistics_per_sample()
    axes[i] = clock.draw_age_prediction_scatterplot(ax=axes[i], fontsize=plt.rcParams["font.size"])
    axes[i].set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+1)
    axes[i].set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+1)
    groupname = dict_name_conversion[list_test_group[i]]
    sample_size = len(df_stat)
    title = f"{groupname} (N={sample_size})"
    axes[i].text(-0.27, 1.1, 
            s=fig_letter, 
            fontdict={"fontsize":plt.rcParams["font.size"]+2, 
                        "fontweight":"bold"}, 
            ha='left', 
            va='top', 
            transform=axes[i].transAxes)
    axes[i].set_title(title, fontsize=plt.rcParams["font.size"]+1, weight="bold", pad=1.0)
    axes[i].set_xlim(10, 90)
    axes[i].set_ylim(10, 90)

file_testing_data = f"{dir_trained_model}/{list_test_group[0]}_std_scaled.txt"
clock = Clock(dir_trained_model, file_testing_data)
clock.draw_feature_importance_barplot(ax4, fontsize=plt.rcParams["font.size"])

ax4.text(-0.47, 
        1.025, 
        s="D", 
        fontdict={"fontsize":plt.rcParams["font.size"]+2, 
                "fontweight":"bold"}, 
        ha='left', 
        va='top', 
        transform=ax4.transAxes)

ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("left")
for label in ax4.get_yticklabels()[:10]:
    label.set_fontweight("bold")

ax4.set_xlabel("Regression\nCoefficient", fontsize=plt.rcParams["font.size"]+1)
ax4.set_ylabel(f"{len(ax4.get_yticklabels())} Significant Age Predictors", fontsize=plt.rcParams["font.size"]+1, rotation=90, weight="bold")
ax4.set_xlim(-5, 5)
ax4.axhline(y=9.5, linewidth=0.5, linestyle="dashed", color="grey")
ax4.legend(loc='upper left', 
            bbox_to_anchor=(1.5, -0.005), 
            ncol=1, 
            frameon=True, 
            borderaxespad=0,
            fontsize=plt.rcParams["font.size"])

name_db1 = "Gene Ontology\n(Biological Process)"
name_db2 = "Molecular Signatures\nDatabase (Cell Type)"
path_db1 = f"{WORKDIR}/Enrichment/Enrichment_GO_BP_216_gene_sets.csv"
path_db2 = f"{WORKDIR}/Enrichment/Enrichment_Celltype_MSigDB_216_gene_sets.csv"
n_max_show = 10
num_minimum_gene = 10
Upper_first = True
no_pathway_id = True
adjust_showing_pathway_length = True
len_showing_adjust_db1 = 30
len_showing_adjust_db2 = 50

size_factor = "nGenes"
size_minmax = (10, 50)
size_show_level = 9

hue_factor = "Neg_Log10_FDR"
palette = "OrRd"

color_label_significant = "black"
color_label_nonsignificant = "gray"

table_db1 = pd.read_csv(path_db1)
table_db2 = pd.read_csv(path_db2)

table_db1["Pathway_id"] = table_db1["Pathway"].apply(lambda x : x.split()[0])
table_db2["Pathway_id"] = table_db2["Pathway"].apply(lambda x : x.split()[0])
table_db1["Pathway_name"] = table_db1["Pathway"].apply(lambda x : ' '.join(x.split()[1:]))
table_db2["Pathway_name"] = table_db2["Pathway"].apply(lambda x : ' '.join(x.split()[1:]))

if Upper_first:
    table_db1["Pathway_name"] = table_db1["Pathway_name"].apply(lambda x : x[0].upper() + x[1:])
    table_db2["Pathway_name"] = table_db2["Pathway_name"].apply(lambda x : x[0].upper() + x[1:])

if no_pathway_id:
    cols_show_pathway = ["Pathway_name"]
else:
    cols_show_pathway = ["Pathway_id", "Pathway_name"]        

table_db1["Pathway_show"] = table_db1[cols_show_pathway].apply(lambda row : ':'.join(row.to_list()), axis = 1)
table_db2["Pathway_show"] = table_db2[cols_show_pathway].apply(lambda row : ':'.join(row.to_list()), axis = 1)

if adjust_showing_pathway_length:
    table_db1["Pathway_show"] = table_db1["Pathway_show"].apply(lambda val : enrich.adjust_phrase_length(val, len_showing_adjust_db1))
    table_db2["Pathway_show"] = table_db2["Pathway_show"].apply(lambda val : enrich.adjust_phrase_length(val, len_showing_adjust_db2))

table_db1_filtered = enrich.filter_table_with_fdr_fe(table_db1, n_max_show, num_minimum_gene)
table_db2_filtered = enrich.filter_table_with_fdr_fe(table_db2, n_max_show, num_minimum_gene)

hue_min = min(table_db1_filtered[hue_factor].to_list() + table_db2_filtered[hue_factor].to_list())
hue_max = max(table_db1_filtered[hue_factor].to_list() + table_db2_filtered[hue_factor].to_list())
size_min = min(table_db1_filtered[size_factor].to_list() + table_db2_filtered[size_factor].to_list())
size_max = max(table_db1_filtered[size_factor].to_list() + table_db2_filtered[size_factor].to_list())

import math

hue_min_stepped = max(math.floor(hue_min / 0.5) * 0.5, 0)
hue_max_stepped = min(math.ceil(hue_max / 0.5) * 0.5, 10)

if hue_min_stepped > hue_max_stepped:
    hue_min_stepped = 0 

size_min_stepped = math.floor(size_min / 0.5) * 0.5
size_max_stepped = math.ceil(size_max / 0.5) * 0.5

# db1
sns.scatterplot(data = table_db1_filtered,
                x = "Fold Enrichment",
                y = "Pathway_show",
                hue = "Neg_Log10_FDR",
                hue_norm = (hue_min_stepped, hue_max_stepped),
                palette = palette,
                edgecolor = "k",
                size = "nGenes",
                size_norm = (size_min_stepped, size_max_stepped),
                sizes = size_minmax,
                ax = ax5,
                legend = False,
                zorder = 3)
tick_colors = list(map({1:color_label_significant, 0:color_label_nonsignificant}.__getitem__, table_db1_filtered["FDR_sig"]))
for ticklabel, tickcolor in zip(ax5.get_yticklabels(), tick_colors):
    ticklabel.set_color(tickcolor)
ax5.grid(linewidth=0.5, zorder = 1)
ax5.set_ylabel("")
ax5.set_xlabel("Fold Enrichment", fontsize=plt.rcParams["font.size"]+1)
ax5.set_xlim(0, 4)
ax5.set_title(name_db1, fontsize = plt.rcParams["font.size"]+1, fontweight="bold", pad=1.0)
ax5.text(-2.0, 1.09, "E", transform = ax5.transAxes, size = plt.rcParams["font.size"]+2, weight = "bold")

# db2
sns.scatterplot(data = table_db2_filtered,
                x = "Fold Enrichment",
                y = "Pathway_show",
                hue = "Neg_Log10_FDR",
                hue_norm = (hue_min_stepped, hue_max_stepped),
                palette = palette,
                edgecolor = "k",
                size = "nGenes",
                size_norm = (size_min_stepped, size_max_stepped),
                sizes = size_minmax,
                ax = ax6,
                legend = False,
                zorder = 3)
tick_colors = list(map({1:color_label_significant, 0:color_label_nonsignificant}.__getitem__, table_db2_filtered["FDR_sig"]))
for ticklabel, tickcolor in zip(ax2.get_yticklabels(), tick_colors):
    ticklabel.set_color(tickcolor)
ax6.grid(linewidth=0.5, zorder = 1)
ax6.set_ylabel("")
ax6.set_xlabel("Fold Enrichment", fontsize=plt.rcParams["font.size"]+1)
ax6.set_xlim(0, 10)
ax6.set_title(name_db2, fontsize = plt.rcParams["font.size"]+1, fontweight="bold", pad=1.0)
ax6.text(-2.0, 1.0, "F", transform = ax6.transAxes, size = plt.rcParams["font.size"]+2, weight = "bold")

# Legend (FDR)
import matplotlib as mpl

cmap = sns.color_palette(palette, as_cmap = True)
norm = mpl.colors.Normalize(vmin=0, vmax=10)
plt.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap),
             cax = ax7,
             orientation = "horizontal",
             label = r"-log$_{10}$FDR",
             )
ax7.xaxis.set_label_position("top")

# Legend (Size)
list_handles = []
list_labels = []
size_step_show = (size_max_stepped - size_min_stepped) / (size_show_level - 1)

for step_multiplier in range(size_show_level):
    size_val = size_min_stepped + size_step_show * step_multiplier
    if size_val > size_max:
        break
    size_val_norm = ((size_val - size_min_stepped) / (size_max_stepped - size_min_stepped)) * \
                    (size_minmax[1] - size_minmax[0]) + size_minmax[0]
    list_handles.append(plt.scatter([],  [], s=size_val_norm, c="gray"))
    list_labels.append(str(int(size_val)))

ax8.legend(
    list_handles,
    list_labels,
    scatterpoints=1,
    ncol=3,  # Reduced column number for better spacing
    title="Number of Genes",
    loc="upper right",
    bbox_to_anchor=(1.08, 0.9),  # Adjusted position to avoid overlap
    frameon=True,
    handletextpad=1.2  # Added padding for readability
)
ax8.set_axis_off()

plt.savefig(f"{WORKDIR}/Figures/Revision_Figure_3.tiff", dpi=300, bbox_inches="tight")
plt.savefig(f"{WORKDIR}/Figures/Revision_Figure_3.png", dpi=300, bbox_inches="tight")
# plt.savefig(f"{WORKDIR}/Figures/Figure_1.tiff", dpi=300, bbox_inches="tight")
# plt.savefig(f"{WORKDIR}/Figures/Figure_1.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()
# %%
