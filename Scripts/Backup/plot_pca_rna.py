# %%
import glob
import os
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pcarna import PCARna
from utility import Utility
from matplotlib.pyplot import get_cmap
import math
import string
import numpy as np

# %%
list_meta_columns = ["Sample_Age_Group",
                    "Sample_Phenotype",
                    "Sample_Sex",
                    "Sampling_Year",
                    "Sequencing_Year",
                    "Sequencing_Platform",
                    "Read_Length"]

dict_conversion = {"readLength": "Read_Length", 
 "Project_Year": "Sampling_Year", 
 "Sequencing_Date": "Sequencing_Year",
 "Sequencing_Type_Platform": "Sequencing_Platform", 
 "Sample_Trait": "Sample_Phenotype"}

# %%
WORKDIR = f"/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_20241024_protein_coding_only"
list_drop_samples = []    
list_select_samples = []
pca = PCARna(list_drop_samples=list_drop_samples, list_select_samples=list_select_samples)
util = Utility()
pcx_num = 1
pcy_num = 2
n_components = 20
colsample = "Project-ID"

# %%
list_study_groups = ["healthy_cohort1", "depression", "GSE273149", "GSE152641"]
list_pca_input = glob.glob(f"{WORKDIR}/**/*.txt", recursive=True)
list_pca_input_filt = list(filter(lambda x: os.path.basename(x).replace(".txt", "") in list_study_groups, list_pca_input))
list_pca_exp_input = sorted(list(filter(lambda x: "FeatureTable" in x, list_pca_input_filt)))
path_target_gene = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_20241024_protein_coding_only/LASSO_INFO/corr_0.32_above_pval_0.05_below/regression_results.txt"
df_target_gene = pd.read_csv(path_target_gene, sep="\t")
list_target_gene = df_target_gene["Selected_Features"].to_list()
# %%
list_df_exp_select = list()
for pca_exp in list_pca_exp_input:
    df_exp = pd.read_csv(pca_exp, sep="\t")
    list_colheader_selected = ["Project-ID"] + list_target_gene
    df_exp_select = df_exp[list_colheader_selected]
    list_df_exp_select.append(df_exp_select)

# %%
df_exp_merged = pd.concat(list_df_exp_select, axis=1)
# %%
df_exp_drop_row = df_exp[~df_exp[colsample].isin(list_drop_samples)]
df_exp_transposed = df_exp_drop_row.set_index(colsample).T
exp_values = df_exp_transposed.values
scaler = StandardScaler().fit(df_exp_transposed)
exp_values_std = scaler.transform(df_exp_transposed)
pcaobj = PCA(n_components=n_components).fit(exp_values_std)
pca_exp = pcaobj.transform(exp_values_std)

# df_meta = pd.read_csv(meta, sep="\t")
# df_meta = df_meta.rename(columns=dict_conversion)
# df_meta.loc[df_meta["Sample_Phenotype"].str.startswith("Mental"), "Sample_Phenotype"] = "Mental illness"
# df_meta.loc[df_meta[colsample].str.startswith("C19-"), "Sample_Phenotype"] = "Infection"
# df_meta_drop = util.drop_samples_meta(df_meta, list_drop_samples, colsample=colsample)
# df_meta_shortened = pca.shorten_meta_info(df_meta_drop)

# df_pca_meta = pca.merge_pcadata_and_metadata(pca_exp, df_meta, colsampleid=colsample)

# %%
figsize= (10, 10)
ncol = 3
plot_linewidth = 1

pcx = f"PC{pcx_num}"
pcy = f"PC{pcy_num}"
nrow = math.ceil(len(list_meta_columns) / ncol)
fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=figsize)  # Increased figure size

# Flatten axes to handle the case when it's a 2D array
if nrow > 1 or ncol > 1:
    axes = axes.flatten()
else:
    axes = np.array([axes]).flatten()

dict_legend_ncol = dict()
for i, (meta, ax) in enumerate(zip(list_meta_columns, axes)):
    fig_letter = list(string.ascii_uppercase)[i]
    val_groupby_meta_pcx = df_pca_meta.groupby(meta)[pcx].apply(np.array).tolist()
    val_groupby_meta_pcy = df_pca_meta.groupby(meta)[pcy].apply(np.array).tolist()
    cmap = get_cmap("Reds", len(val_groupby_meta_pcx))
    list_colors = [cmap(ind) for ind in range(len(val_groupby_meta_pcx))]
    dim_level = 10
    list_colors_cmyk = [pca.rgb_to_cmyk(rgb[0], rgb[1], rgb[2]) for rgb in list_colors]
    list_colors_cmyk_dim = [(cmyk[0], cmyk[1], cmyk[2], cmyk[3] + dim_level) for cmyk in list_colors_cmyk]
    list_colors_rgb_dim = [pca.cmyk_to_rgb(cmyk[0], cmyk[1], cmyk[2], cmyk[3]) for cmyk in list_colors_cmyk_dim]
    list_legends = df_pca_meta.groupby(meta)[pcy].apply(np.array).index.tolist()
    list_legends = [int(float(x)) if str(x).endswith(".0") else x for x in list_legends]
    legend_ncol = dict_legend_ncol.get(meta, 1)

    for (x, y, color, legend) in zip(val_groupby_meta_pcx, val_groupby_meta_pcy, list_colors_rgb_dim, list_legends):
        ax.scatter(x, y, c=color, s=plot_linewidth, label=legend, alpha=0.3, zorder=-11) 
        pca.confidence_ellipse(x, y, ax, n_std=3.0, edgecolor=color, facecolor='none', linewidth=plot_linewidth*2.0, alpha=0.5)
        ax.scatter(np.mean(x), np.mean(y), c=color, edgecolor="k", linewidth=plot_linewidth*0.5, marker="*", s=plot_linewidth*70, alpha=1, zorder=999)  # Increased mean marker size

    ax.set_title(meta, 
                fontsize=plt.rcParams["font.size"]+1,
                fontweight="bold",
                pad=2)  # Adjusted padding to avoid overlap with y-tick labels
    ax.set_ylabel(f"principal component {pcy_num}\n({round(pcaobj.explained_variance_ratio_[pcy_num - 1] * 100, 1)}%)", 
                fontsize=plt.rcParams["font.size"]+1, 
                labelpad=2)  # Added label padding
    ax.set_xlabel(f"principal component {pcx_num}\n({round(pcaobj.explained_variance_ratio_[pcx_num - 1] * 100, 1)}%)", 
                fontsize=plt.rcParams["font.size"]+1, 
                labelpad=2)  # Added label padding
            
    legend = ax.legend(loc="upper right", 
                    fontsize=plt.rcParams["font.size"],
                    ncol=legend_ncol, 
                    markerscale=plot_linewidth,
                    handletextpad=1.0, 
                    bbox_to_anchor=(1.0, 1.0), 
                    frameon=True)

    plt.setp(legend.get_texts(), weight='bold')

    for handle in legend.legendHandles:
        handle.set_alpha(1)
    
    ax.annotate(fig_letter,
                xy=(-0.4, 1.04),
                xycoords='axes fraction',
                xytext=(0, 0),
                textcoords='offset points',
                size=plt.rcParams["font.size"]+2, 
                ha='left', 
                va='center',
                fontweight="bold",
                color="black")
    
for extra_ind in range(len(list_meta_columns), len(axes)):
    axes[extra_ind].axis("off")
    
plt.tight_layout()
plt.show()
plt.close()

# %%
