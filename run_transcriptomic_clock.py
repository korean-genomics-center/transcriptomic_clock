# %%
import os

import draw_enrichment_plot as enrich
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib.lines import Line2D

from transcriptomic_clock import Clock

# %%
dict_name_conversion = {"healthy_train": "Healthy Training Set (N=370)",
                        "healthy": "Healthy Validation Set (N=93)", 
                        "healthy_bgi": "Healthy Test Set (N=96)", 
                        "anxiety": "Anxiety", 
                        "depression": "Major Depression", 
                        "suicide_attempt": "Suicide Attempt",
                        "viral_v1_long": "Acute Phase",
                        "viral_v2_long": "Mid Phase",
                        "viral_v3_long": "Late Phase",
                        "viral_recov": "Convalescent"}


# %% [draw_feature_importance_barplot]
cm = 1/2.54
path = "./Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
width = 12*cm
height = 12*cm
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

gs4 = gsfig[10:200, 55:70]
ax4 = fig.add_subplot(gs4) 

gs5 = gsfig[15:85, 140:170]
ax5 = fig.add_subplot(gs5)

gs6 = gsfig[110:160, 140:170]
ax6 = fig.add_subplot(gs6)

gs7 = gsfig[180:185, 110:160]
ax7 = fig.add_subplot(gs7)

gs8 = gsfig[190:200, 110:158]
ax8 = fig.add_subplot(gs8)

# gs9 = gsfig[165:210, 95:175]
# ax9 = fig.add_subplot(gs9)


list_test_group = ["healthy_train", "healthy", "healthy_bgi"]
dir_trained_model = f"./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below"
file_testing_data = f"{dir_trained_model}/testing_dataset_selected_{list_test_group[0]}.tsv"
os.makedirs(dir_trained_model, exist_ok=True)
clock = Clock(dir_trained_model, file_testing_data)
clock.make_dataframe_error_statistics_per_sample()
ax1 = clock.draw_age_prediction_scatterplot(ax=ax1, fontsize=plt.rcParams["font.size"])
ax1.set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+1)
ax1.set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+1)
# Set title for each subplot
title = dict_name_conversion[list_test_group[0]]
ax1.text(-0.27, 1.1, 
        s="a", 
        fontdict={"fontsize":plt.rcParams["font.size"]+2, 
                    "fontweight":"bold"}, 
        ha='left', 
        va='top', 
        transform=ax1.transAxes)
ax1.set_title(title, fontsize=plt.rcParams["font.size"]+1, weight="bold", pad=1.0)
ax1.set_xlim(10, 90)
ax1.set_ylim(10, 90)

file_testing_data = f"{dir_trained_model}/testing_dataset_selected_{list_test_group[1]}.tsv"
os.makedirs(dir_trained_model, exist_ok=True)
clock = Clock(dir_trained_model, file_testing_data)
clock.make_dataframe_error_statistics_per_sample()
ax2 = clock.draw_age_prediction_scatterplot(ax=ax2, fontsize=plt.rcParams["font.size"])
ax2.set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+1)
ax2.set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+1)
# Set title for each subplot
title = dict_name_conversion[list_test_group[1]]
ax2.text(-0.27, 1.1, 
        s="b", 
        fontdict={"fontsize":plt.rcParams["font.size"]+2, 
                    "fontweight":"bold"}, 
        ha='left', 
        va='top', 
        transform=ax2.transAxes)
ax2.set_title(title, fontsize=plt.rcParams["font.size"]+1, weight="bold", pad=1.0)
ax2.set_xlim(10, 90)
ax2.set_ylim(10, 90)


file_testing_data = f"{dir_trained_model}/testing_dataset_selected_{list_test_group[2]}.tsv"
os.makedirs(dir_trained_model, exist_ok=True)
clock = Clock(dir_trained_model, file_testing_data)
clock.make_dataframe_error_statistics_per_sample()
ax3 = clock.draw_age_prediction_scatterplot(ax=ax3, fontsize=plt.rcParams["font.size"])
ax3.set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+1)
ax3.set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+1)
# Set title for each subplot
title = dict_name_conversion[list_test_group[2]]
ax3.text(-0.27, 1.1, 
        s="c", 
        fontdict={"fontsize":plt.rcParams["font.size"]+2, 
                    "fontweight":"bold"}, 
        ha='left', 
        va='top', 
        transform=ax3.transAxes)
ax3.set_title(title, fontsize=plt.rcParams["font.size"]+1, weight="bold", pad=1.0)
ax3.set_xlim(10, 90)
ax3.set_ylim(10, 90)

file_testing_data = f"{dir_trained_model}/testing_dataset_selected_healthy.tsv"
clock = Clock(dir_trained_model, file_testing_data)
clock.draw_feature_importance_barplot(ax4, fontsize=plt.rcParams["font.size"])

# Adding text
ax4.text(-0.47, 
        1.025, 
        s="d", 
        fontdict={"fontsize":plt.rcParams["font.size"]+2, 
                "fontweight":"bold"}, 
        ha='left', 
        va='top', 
        transform=ax4.transAxes)

# Move y-tick labels to the right side
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("left")
for label in ax4.get_yticklabels()[:10]:
    label.set_fontweight("bold")
# Set labels
ax4.set_xlabel("Regression\nCoefficient", fontsize=plt.rcParams["font.size"]+1)
ax4.set_ylabel(f"{len(ax4.get_yticklabels())} Significant Age Predictors", fontsize=plt.rcParams["font.size"]+1, rotation=90, weight="bold")

# Setting limits
ax4.set_xlim(-5, 5)

ax4.axhline(y=9.5, linewidth=0.5, linestyle="dashed", color="grey")
ax4.legend(loc='upper left', 
            bbox_to_anchor=(1.25, -0.005), 
            ncol=1, 
            frameon=True, 
            borderaxespad=0,
            fontsize=plt.rcParams["font.size"])

name_db1 = "Gene Ontology\n(Biological Process)"
name_db2 = "Molecular Signatures Database\n(Cell Type)"
path_db1 = "Enrichment/282_GOBP_10MinPathway_001FDR.csv"
path_db2 = "Enrichment/282_MSigDB_10MinPathway_001FDR.csv"
n_max_show = 10
Upper_first = True
no_pathway_id = True
adjust_showing_pathway_length = True
len_showing_adjust_db1 = 30
len_showing_adjust_db2 = 50

size_factor = "nGenes"
size_minmax = (10, 50)
size_show_level = 6

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
    table_db2["Pathway_name"] = table_db2["Pathway_name"].apply(lambda x : x[0].upper() + x[1:] if len(x) < 30 else ' '.join(x.split(' ')[-2:]))

if no_pathway_id:
    cols_show_pathway = ["Pathway_name"]
else:
    cols_show_pathway = ["Pathway_id", "Pathway_name"]        

table_db1["Pathway_show"] = table_db1[cols_show_pathway].apply(lambda row : ':'.join(row.to_list()), axis = 1)
table_db2["Pathway_show"] = table_db2[cols_show_pathway].apply(lambda row : ':'.join(row.to_list()), axis = 1)

if adjust_showing_pathway_length:
    table_db1["Pathway_show"] = table_db1["Pathway_show"].apply(lambda val : enrich.adjust_phrase_length(val, len_showing_adjust_db1))
    table_db2["Pathway_show"] = table_db2["Pathway_show"].apply(lambda val : enrich.adjust_phrase_length(val, len_showing_adjust_db2))

table_db1_filtered = enrich.filter_table_with_fdr_fe(table_db1, n_max_show)
table_db2_filtered = enrich.filter_table_with_fdr_fe(table_db2, n_max_show)

hue_min = min(table_db1_filtered[hue_factor].to_list() + table_db2_filtered[hue_factor].to_list())
hue_max = max(table_db1_filtered[hue_factor].to_list() + table_db2_filtered[hue_factor].to_list())
size_min = min(table_db1_filtered[size_factor].to_list() + table_db2_filtered[size_factor].to_list())
size_max = max(table_db1_filtered[size_factor].to_list() + table_db2_filtered[size_factor].to_list())

import math

# Calculate hue_min_stepped
hue_min_stepped = max(math.floor(hue_min / 0.5) * 0.5, 0)
# Ensure hue_max_stepped is set appropriately based on your data range
# Example: Adjust to cover your expected range if necessary
hue_max_stepped = min(math.ceil(hue_max / 0.5) * 0.5, 10)  # Example range adjustment

# Ensure hue_min_stepped is less than or equal to hue_max_stepped
if hue_min_stepped > hue_max_stepped:
    hue_min_stepped = 0  # or any other default value you prefer

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
ax5.set_xlim(0, 40)
ax5.set_title(name_db1, fontsize = plt.rcParams["font.size"]+1, fontweight="bold", pad=1.0)
ax5.text(-1.5, 1.09, "e", transform = ax5.transAxes, size = plt.rcParams["font.size"]+2, weight = "bold")

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
ax6.set_xlim(0, 20)
ax6.set_title(name_db2, fontsize = plt.rcParams["font.size"]+1, fontweight="bold", pad=1.0)
ax6.text(-1.5, 1.0, "f", transform = ax6.transAxes, size = plt.rcParams["font.size"]+2, weight = "bold")

# Legend (FDR)
import matplotlib as mpl

cmap = sns.color_palette(palette, as_cmap = True)
norm = mpl.colors.Normalize(vmin=0, vmax=hue_max_stepped)
plt.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap),
             cax = ax7,
             orientation = "horizontal",
             label = r"statistcal signficance (-log$_{10}$FDR)",
             anchor = (0.5, 1),
             panchor = (1, 1))
ax7.xaxis.set_label_position("top")

# Legend (Size)
list_handles = list()
list_labels = list()
size_step_show = math.ceil((size_max_stepped - size_min_stepped) / (size_show_level-1))
for step_multiplier in range(size_show_level):
    size_val = size_min_stepped + size_step_show * step_multiplier
    if size_val > size_max:
        break
    size_val_norm = size_val * ((size_minmax[1] - size_minmax[0]) / (size_max_stepped - size_min_stepped))
    list_handles.append(plt.scatter([],[], s = size_val_norm, c = "gray"))
    list_labels.append(str(int(size_val)))
ax8.legend(list_handles,
           list_labels,
           scatterpoints = 1,
           ncol = 3 if len(list_handles) > 4 else 2,
           title = "Number of Genes",
           loc = "upper right",
           bbox_to_anchor = (1, 0.8),
           labelspacing = 1,
           frameon = True)
ax8.set_axis_off()

# Draw Rectangle
# rec = mpl.patches.Rectangle((0,0), 1, 1, linewidth = 1, edgecolor = "gray", facecolor = "none")
# ax9.add_patch(rec)
# ax9.set_axis_off()

plt.savefig("Figures/Fig_1.pdf", dpi=300, bbox_inches="tight")
plt.savefig("Figures/Fig_1.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# %%
cm = 1/2.54
path = "/BiO/Access/kyungwhan1998/miniconda3/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
width = 8*cm
height = 16*cm

list_test_group = ["viral_v1_long", 
                   "viral_v2_long", 
                   "viral_v3_long", 
                   "viral_recov", 
                   "anxiety",
                   "suicide_attempt", 
                   "depression"]

# Create a figure to hold all subplots
fig, axs = plt.subplots(4, 2, figsize=(width, height), constrained_layout=True)

# Flatten the axs array to iterate over it easily
axs = axs.flatten()

# Iterate over each test group
for i, testname in enumerate(list_test_group):
    try:
        dir_trained_model = f"./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below"
        file_testing_data = f"{dir_trained_model}/testing_dataset_selected_{testname}.tsv"
        os.makedirs(dir_trained_model, exist_ok=True)
        
        # Assuming Clock class and its methods (make_dataframe_error_statistics_per_sample, draw_age_prediction_scatterplot) exist
        clock = Clock(dir_trained_model, file_testing_data)
        df_clock = clock.make_dataframe_error_statistics_per_sample()
        df_clock.to_csv(os.path.join(dir_trained_model, f"sample_statistics_{testname}.tsv"), sep="\t", index=False)
        print(os.path.join(dir_trained_model, f"sample_statistics_{testname}.tsv"))

        
        # Plot age prediction scatterplot for the current test group
        ax = clock.draw_age_prediction_scatterplot(ax=axs[i], fontsize=plt.rcParams["font.size"])
        ax.set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+1)
        ax.set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+1)
        # Set title for each subplot
        title = dict_name_conversion[testname]
        figletter=chr(ord('a') + i)
        ax.text(-0.32, 1.2, 
                s=figletter, 
                fontdict={"fontsize":plt.rcParams["font.size"]+2, 
                          "fontweight":"bold"}, 
                ha='left', 
                va='top', 
                transform=ax.transAxes)
        ax.legend().set_visible(False)
        axs[i].set_title(f"{title} (N={len(df_clock)})", fontsize=plt.rcParams["font.size"]+1, weight="bold", pad=0.5)
        ax.set_xlim(10, 90)
        ax.set_ylim(10, int(ax.get_yticklabels()[-1].get_text())+5)
    except Exception as e:
        print(f"Error processing {testname}: {str(e)}")

# Hide unused subplots
for j in range(len(list_test_group), len(axs)):
    axs[j].axis('off')

# Extract legend handles and labels from one of the subplots
handles, labels = axs[0].get_legend_handles_labels()
# Create a new legend in the last subplot
axs[-1].legend(handles, labels, loc='center', fontsize=plt.rcParams["font.size"]+1)

# Adjust layout and show plot
plt.savefig("Figures/Extended_Data_Fig_6.png", dpi=300)
plt.savefig("Figures/Extended_Data_Fig_6.pdf", dpi=300)
plt.show()
plt.close()

# %%
def make_df_perfomance_evaluation(workdir, testing_data, list_corr_thres):
    dict_annot_age_eval_plot = dict()
    for corr_thres in list_corr_thres:
        dir_trained_model = f"{workdir}/corr_{corr_thres}_above_pval_0.05_below"
        file_testing_data = f"{dir_trained_model}/{testing_data}"
        clock = Clock(dir_trained_model, file_testing_data)
        print(f"----{corr_thres}----")
        try:
            clock.make_dataframe_clock_features()
            dict_annot_corr_thres = clock.get_annotation_prediction_evaluation_plot()
            dict_annot_age_eval_plot[corr_thres] = dict_annot_corr_thres
        except Exception as e:
            print(e)
            print(f"No data for {corr_thres}... Skipping...")
    df_annot_age_eval_plot = pd.DataFrame(dict_annot_age_eval_plot)

    return df_annot_age_eval_plot

# %% [draw_plot_performance_profile]
# cm = 1/2.54
# path = "/BiO/Access/kyungwhan1998/miniconda3/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf"
# prop = fm.FontProperties(fname=path)
# plt.rcParams['font.family'] = prop.get_name()
# plt.rcParams["font.size"] = 5
# width = 10*cm
# height = 5*cm
# figsize = (width, height)
# plot_linewidth = 1

# WORKDIR = "./LASSO_INFO/standardized/healthy_illumina"
# dir_trained_model = f"{WORKDIR}/corr_0.3_above_pval_0.05_below"
# list_corr_thres = list(map(lambda x: round(float(x), 2), list(np.arange(0.3, 0.51, 0.01))))
# clock = Clock(dir_trained_model, "")
# df_annot_age_eval_plot_train = make_df_perfomance_evaluation(WORKDIR, "training_dataset.tsv", list_corr_thres)
# df_annot_age_eval_plot_test1 = make_df_perfomance_evaluation(WORKDIR, f"testing_dataset_selected_healthy.tsv", list_corr_thres)
# col_train, mae_train, corr_train, list_num_input_genes_train, list_num_features_train = clock.summarize_performance_statistics(df_annot_age_eval_plot_train)
# col_test, mae_test, corr_test, list_num_input_genes_test, list_num_features_test = clock.summarize_performance_statistics(df_annot_age_eval_plot_test1)
# fig, ax1, ax2 = clock.draw_plot_performance_profile(col_train, 
#                                     col_test, 
#                                     mae_train, 
#                                     mae_test, 
#                                     corr_train, 
#                                     corr_test, 
#                                     list_num_features_train, 
#                                     list_num_input_genes_train,
#                                     plot_linewidth=plot_linewidth, 
#                                     figsize=figsize)
# # Annotate
# ax1.annotate('Thresholds\n\nNo. LASSO Features\n\nNo. Input Genes', 
#              xy=(0.01, 1.3), 
#              xycoords=ax1.get_yaxis_transform(),
#              xytext=(0, 0), 
#              textcoords="offset points", 
#              ha="right", 
#              va="center", 
#              fontsize=plt.rcParams["font.size"], 
#              weight="bold")

# # Grid and labels
# ax1.grid(axis="both", linestyle="dashed", linewidth=plot_linewidth-0.5)
# ax1.set_ylabel('MAE', color='darkred', fontsize=plt.rcParams["font.size"]+2, labelpad=10)
# ax2.set_ylabel("Pearson's r", color='darkblue', fontsize=plt.rcParams["font.size"]+2, labelpad=10)  # we already handled the x-label with ax1

# # Add vertical lines
# plt.axvline(x=0.30, color="grey", linewidth=plot_linewidth)
# plt.axvline(x=0.37, color="grey", linewidth=plot_linewidth)

# # Legend
# legend_elements = [
#                 Line2D([0], [0], color='k', linestyle="dotted", linewidth=plot_linewidth, label='Training'),
#                 Line2D([0], [0], color='k', linestyle="solid", linewidth=plot_linewidth, label='Validation')
#                 ]
# plt.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(1, 1), prop={'size': plt.rcParams["font.size"]}, frameon=False)

# # Final adjustments
# plt.margins(x=0.1)
# plt.tight_layout()

# # Show the plot
# plt.show()
# plt.close()

# %%

cm = 1/2.54
path = "/BiO/Access/kyungwhan1998/miniconda3/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
width = 10/1.4*cm
height = 14/1.4*cm
figsize = (width, height)
plot_linewidth = 3

fig = plt.figure(figsize = figsize)
col = 100
row = 140
gsfig = gridspec.GridSpec(
    row, col,
    left = 0, right = 1, bottom = 0,
    top = 1, wspace = 1, hspace = 1
)

gs1 = gsfig[0:95, 0:30]
ax1 = fig.add_subplot(gs1)

gs2 = gsfig[0:95, 70:100]
ax2 = fig.add_subplot(gs2)

gs3 = gsfig[115:120, 0:30] 
ax3 = fig.add_subplot(gs3)

gs4 = gsfig[105:120, 60:90]
ax4 = fig.add_subplot(gs4) 


name_db1 = "Gene Ontology\n(Cellular Component)"
name_db2 = "KEGG"
path_db1 = "Enrichment/3407_GOCC_100MinPathway_001FDR.csv"
path_db2 = "Enrichment/3407_kegg_100MinPathway_001FDR.csv"
n_max_show = 15
Upper_first = True
no_pathway_id = True
adjust_showing_pathway_length = True
len_showing_adjust_db1 = 30
len_showing_adjust_db2 = 50

size_factor = "nGenes"
size_minmax = (10, 50)
size_show_level = 10

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

table_db1_filtered = enrich.filter_table_with_fdr_fe(table_db1, n_max_show)
table_db2_filtered = enrich.filter_table_with_fdr_fe(table_db2, n_max_show)

hue_min = min(table_db1_filtered[hue_factor].to_list() + table_db2_filtered[hue_factor].to_list())
hue_max = max(table_db1_filtered[hue_factor].to_list() + table_db2_filtered[hue_factor].to_list())
size_min = min(table_db1_filtered[size_factor].to_list() + table_db2_filtered[size_factor].to_list())
size_max = max(table_db1_filtered[size_factor].to_list() + table_db2_filtered[size_factor].to_list())

import math

# Calculate hue_min_stepped
hue_min_stepped = max(math.floor(hue_min / 0.5) * 0.5, 0)
# Ensure hue_max_stepped is set appropriately based on your data range
# Example: Adjust to cover your expected range if necessary
hue_max_stepped = min(math.ceil(hue_max / 0.5) * 0.5, 10)  # Example range adjustment

# Ensure hue_min_stepped is less than or equal to hue_max_stepped
if hue_min_stepped > hue_max_stepped:
    hue_min_stepped = 0  # or any other default value you prefer

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
                ax = ax1,
                legend = False,
                zorder = 3)
tick_colors = list(map({1:color_label_significant, 0:color_label_nonsignificant}.__getitem__, table_db1_filtered["FDR_sig"]))
for ticklabel, tickcolor in zip(ax1.get_yticklabels(), tick_colors):
    ticklabel.set_color(tickcolor)
ax1.grid(linewidth=0.5, zorder = 1)
ax1.set_ylabel("")
ax1.set_xlabel("Fold Enrichment", fontsize=plt.rcParams["font.size"]+1)
ax1.set_xlim(0, 4)
ax1.set_title(name_db1, fontsize = plt.rcParams["font.size"]+1, fontweight="bold", pad=1)
ax1.text(-1.2, 1.02, "a", transform = ax1.transAxes, size = plt.rcParams["font.size"]+2, weight = "bold")

# db2
import textwrap


# Define a function to wrap text
def wrap_labels(labels, width):
    return ['\n'.join(textwrap.wrap(label, width)) for label in labels]

# Wrap the 'Pathway_show' labels
wrapped_labels = wrap_labels(table_db2_filtered['Pathway_show'], width=20)
table_db2_filtered['Wrapped_Pathway_show'] = wrapped_labels

# Create the scatter plot with wrapped labels
sns.scatterplot(data=table_db2_filtered,
                x="Fold Enrichment",
                y="Wrapped_Pathway_show",
                hue="Neg_Log10_FDR",
                hue_norm=(hue_min_stepped, hue_max_stepped),
                palette=palette,
                edgecolor="k",
                size="nGenes",
                size_norm=(size_min_stepped, size_max_stepped),
                sizes=size_minmax,
                ax=ax2,
                legend=False,
                zorder=3)
tick_colors = list(map({1:color_label_significant, 0:color_label_nonsignificant}.__getitem__, table_db2_filtered["FDR_sig"]))
for ticklabel, tickcolor in zip(ax2.get_yticklabels(), tick_colors):
    ticklabel.set_color(tickcolor)
ax2.grid(linewidth=0.5, zorder = 1)
ax2.set_ylabel("")
ax2.set_xlabel("Fold Enrichment", fontsize=plt.rcParams["font.size"]+1)
ax2.set_xlim(0, 4)
ax2.set_title(name_db2, fontsize = plt.rcParams["font.size"]+1, fontweight="bold", pad=1)
ax2.text(-1.2, 1.02, "b", transform = ax2.transAxes, size = plt.rcParams["font.size"]+2, weight = "bold")

# Legend (FDR)
import matplotlib as mpl

cmap = sns.color_palette(palette, as_cmap = True)
norm = mpl.colors.Normalize(vmin=0, vmax=hue_max_stepped)
plt.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmap),
             cax = ax3,
             orientation = "horizontal",
             label = r"-log$_{10}$FDR",
             anchor = (0.5, 1),
             panchor = (1, 1))
ax3.xaxis.set_label_position("top")

# Legend (Size)
list_handles = list()
list_labels = list()
size_step_show = math.ceil((size_max_stepped - size_min_stepped) / (size_show_level-1))
for step_multiplier in range(size_show_level):
    size_val = size_min_stepped + size_step_show * step_multiplier
    if size_val > size_max:
        break
    size_val_norm = size_val * ((size_minmax[1] - size_minmax[0]) / (size_max_stepped - size_min_stepped))
    list_handles.append(plt.scatter([],[], s = size_val_norm, c = "gray"))
    list_labels.append(str(int(size_val)))
ax4.legend(list_handles,
           list_labels,
           scatterpoints = 1,
           ncol = 3 if len(list_handles) > 4 else 2,
           title = "Number of Genes",
           loc = "upper right",
           bbox_to_anchor = (1, 0.8),
           frameon = True)
ax4.set_axis_off()

# # Draw Rectangle
# rec = mpl.patches.Rectangle((0,0), 1, 1, linewidth = 1, edgecolor = "gray", facecolor = "none")
# ax5.add_patch(rec)
# ax5.set_axis_off()

plt.savefig("Figures/Extended_Data_Fig_2.pdf", dpi=300, bbox_inches="tight")
plt.savefig("Figures/Extended_Data_Fig_2.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# %%
