# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
path_reg = "./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below/regression_results.txt"

df_reg = pd.read_csv(path_reg,sep="\t")
df_reg_pos = df_reg[df_reg["Regression_Coefficient"] > 0]
df_reg_neg = df_reg[df_reg["Regression_Coefficient"] < 0]
list_features = df_reg["Selected_Features"].to_list()
list_pos_features = df_reg_pos["Selected_Features"].to_list()
list_neg_features = df_reg_neg["Selected_Features"].to_list()
# %%
path_norm_exp = "./FeatureTable/healthy_illumina.txt"

def get_end_idx_gene_list(df_feature: pd.DataFrame, gene_tag: str) -> int:
    list_header = list(df_feature.columns)
    list_gene_idx = list()
    for idx, column in enumerate(list_header):
        if str(column).startswith(gene_tag):
            list_gene_idx.append(idx)
    end_idx_gene = int(max(list_gene_idx))

    return end_idx_gene
df_norm_exp = pd.read_csv(path_norm_exp, sep="\t")
df_norm_exp_set_idx = df_norm_exp.set_index("SUBJID")
end_tag = get_end_idx_gene_list(df_norm_exp_set_idx, "ENSG")
df_norm_exp_set_idx_gene_only = df_norm_exp_set_idx.iloc[:,:end_tag+1]
df_meta = df_norm_exp_set_idx.iloc[:,end_tag+1:]
corr_matrix = df_norm_exp_set_idx_gene_only.corr(method="pearson")

# %%
from collections import Counter
dict_gene_to_highcorr_cnt = corr_matrix.apply(lambda row: sum(list(map(lambda val: val > 0.68, row)))-1, axis = 1).to_dict()
sorted(dict_gene_to_highcorr_cnt.items(), key = lambda x : x[1], reverse = True)
plt.hist(dict_gene_to_highcorr_cnt.values(), color="grey")
plt.ylabel("Count", fontsize=16)
plt.xlabel("Co-Expressed Genes", fontsize=16)
plt.show()
plt.close()

# %%
list_col = list(map(lambda x: x.split(".")[0] + "_" + x.split("_")[-1], corr_matrix.columns))
corr_matrix_selected = corr_matrix.copy()
corr_matrix_selected.columns = list_col
corr_matrix_selected.index = list_col
corr_matrix_selected =corr_matrix_selected[list_features]

# %%
from matplotlib.ticker import MaxNLocator
import math
from collections import Counter

range_thres = np.arange(0.50, 0.70, 0.01)
fig, axes = plt.subplots(nrows=math.ceil(np.sqrt(len(range_thres))), ncols=math.floor(np.sqrt(len(range_thres))),figsize=(10, 10), sharex=True, sharey=True)
for thres, ax in zip(range_thres, axes.flatten()):
    
    dict_marker_highcorr = dict()
    for featname in corr_matrix_selected.columns:
        dict_gene_to_corr = corr_matrix_selected.loc[:, featname].to_dict()
        list_gene_highcorr = list(filter(lambda genename: dict_gene_to_corr[genename] > thres, dict_gene_to_corr.keys()))
        dict_marker_highcorr[featname] = list_gene_highcorr

    list_diff = list()
    for featname, list_highcorr in dict_marker_highcorr.items():
        # print(featname, len(list_highcorr), len(set(list_highcorr) - set(list_features)))  
        diff = (len(list_highcorr) - len(set(list_highcorr) - set(list_features)))
        list_diff.append(diff)
    
    dict_cnt = Counter(list_diff)

    ax.bar(list(map(lambda x: int(x)-1,dict_cnt.keys())), list(map(int, dict_cnt.values())), color="grey")
    ax.set_title(f'Threshold: {thres:.2f}')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

fig.supylabel("Count", fontsize=16)
fig.supxlabel("Number of Overlaps", fontsize=16)
fig.tight_layout()
plt.show()
plt.close()

# %%
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Assuming corr_matrix_selected is already defined
# Create the clustermap with the default settings
cg = sns.clustermap(
    corr_matrix_selected.T,
    cmap="coolwarm",
    xticklabels=False,
    figsize=(5, 15),
    vmin=-1,
    vmax=1
)

# Update yticklabels
cg.ax_heatmap.set_yticklabels(
    [r'$\it{' + ticklabel.get_text().split("_")[-1] + '}$' for ticklabel in cg.ax_heatmap.get_yticklabels()], weight="bold"
)
cg.ax_heatmap.set_xlabel("12,546 genes", fontsize=12)

# Hide the column dendrogram
cg.ax_col_dendrogram.set_visible(False)

# Adjust the layout of the row dendrogram to fit the heatmap
heatmap_pos = cg.ax_heatmap.get_position()
dendrogram_pos = cg.ax_row_dendrogram.get_position()

# Adjust the row dendrogram to match the heatmap height
cg.ax_row_dendrogram.set_position([dendrogram_pos.x0, heatmap_pos.y0, dendrogram_pos.width, heatmap_pos.height])

# Remove the default colorbar
cg.cax.set_visible(False)

# Create a new horizontal colorbar at the bottom
cbar_ax = cg.fig.add_axes([heatmap_pos.x0, heatmap_pos.y0 - 0.05, heatmap_pos.width, 0.015])
cbar = plt.colorbar(cg.ax_heatmap.collections[0], cax=cbar_ax, orientation='horizontal')
cbar.set_label('Gene Co-Expression')

# Add ticklabels as the value of the correlation
cbar.set_ticks(np.linspace(-1, 1, 5))
cbar.set_ticklabels(['-1', '-0.5', '0', '0.5', '1'])

# Customize colorbar appearance to move ticks to the bottom
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.ax.xaxis.set_label_position('bottom')
cbar.ax.tick_params(axis='x', length=5)
for spine in cbar.ax.spines:
    cbar.ax.spines[spine].set_color('black')
    cbar.ax.spines[spine].set_linewidth(0.1)

plt.show()

# %%
from itertools import chain
dict_feat_module = dict()
for feature in list_features:
    corr_feature = corr_matrix_selected[feature]
    corr_feature_abs = corr_feature.apply(abs)
    corr_feature_sort = corr_feature_abs.sort_values(ascending=False)
    corr_feature_final = corr_feature_sort.head(6)
    dict_feat_module["_".join(feature.split("_")[1:])] = list(map(lambda x: x.split("_")[-1], list(corr_feature_final.keys())))

list(chain(*list(dict_feat_module.values())))


# %%
from itertools import chain
dict_feat_module = dict()
for feature in list_features:
    corr_feature = corr_matrix_selected[feature]
    corr_feature_abs = corr_feature.apply(abs)
    corr_feature_filt = corr_feature_abs[corr_feature_abs > 0.50]
    dict_feat_module[feature] = list(map(lambda x: x.split("_")[-1], list(corr_feature_filt.keys())))

list(chain(*list(dict_feat_module.values())))
# %%
