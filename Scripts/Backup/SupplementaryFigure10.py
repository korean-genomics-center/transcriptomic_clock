# %%
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
import draw_enrichment_plot as enrich

# %%
WORKDIR = ".."
cm = 1/2.54
path = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
width = 10*cm
height = 20*cm
figsize = (width, height)
plot_linewidth = 3

fig = plt.figure(figsize = figsize)
col = 100
row = 200
gsfig = gridspec.GridSpec(
    row, col,
    left = 0, right = 1, bottom = 0,
    top = 1, wspace = 1, hspace = 1
)

gs1 = gsfig[0:95, 0:30]
ax1 = fig.add_subplot(gs1)

gs2 = gsfig[0:95, 60:90]
ax2 = fig.add_subplot(gs2)

gs3 = gsfig[110:115, 0:30] 
ax3 = fig.add_subplot(gs3)

gs4 = gsfig[105:120, 60:90]
ax4 = fig.add_subplot(gs4) 


name_db1 = "Gene Ontology\n(Cellular Component)"
name_db2 = "KEGG"
path_db1 = f"{WORKDIR}/Enrichment/Enrichment_Celltype_MSigDB_7608_gene_sets.csv"
path_db2 = f"{WORKDIR}/Enrichment/Enrichment_GO_BP_7608_gene_sets.csv"
n_max_show = 15
num_minimum_gene = 10
Upper_first = True
no_pathway_id = True
adjust_showing_pathway_length = True
len_showing_adjust_db1 = 30
len_showing_adjust_db2 = 30

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
ax1.text(-1.0, 1.02, "A", transform = ax1.transAxes, size = plt.rcParams["font.size"]+2, weight = "bold")

# db2
import textwrap


# Define a function to wrap text
def wrap_labels(labels, width):
    return ['\n'.join(textwrap.wrap(label, width)) for label in labels]

# Wrap the 'Pathway_show' labels
wrapped_labels = wrap_labels(table_db2_filtered['Pathway_show'], width=30)
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
ax2.text(-0.9, 1.02, "B", transform = ax2.transAxes, size = plt.rcParams["font.size"]+2, weight = "bold")

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

plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_10.tiff", dpi=300, bbox_inches="tight")
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_10.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# %%
