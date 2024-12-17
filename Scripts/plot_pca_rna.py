# %%
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pcarna import PCARna
from matplotlib.pyplot import get_cmap
import math
import string
import numpy as np

# %%
dict_conversion = {"readLength": "Read_Length", 
                    "Project_Year": "Sampling_Year", 
                    "Sequencing_Date": "Sequencing_Year",
                    "Sequencing_Type_Platform": "Sequencing_Platform", 
                    "Sample_Trait": "Sample_Phenotype"}
list_meta_columns = list(dict_conversion.values())
figsize= (7, 5)
ncol = 1
plot_linewidth = 1
n_components = 20
list_drop_samples = []    
list_select_samples = []
pca = PCARna(list_drop_samples=list_drop_samples, list_select_samples=list_select_samples)

# %%
def get_pcaobj_and_transformed_data(df_train, df_test, n_components):
    # Transpose input data for PCA
    df_train_transposed = df_train.T
    df_test_transposed = df_test.T
    
    # Standardize the training data and apply the same scaler to the test data
    scaler = StandardScaler().fit(df_train_transposed)
    train_scaled = scaler.transform(df_train_transposed)
    test_scaled = scaler.transform(df_test_transposed)
    
    # Fit PCA on the training data
    pcaobj = PCA(n_components=n_components).fit(train_scaled)
    
    # Transform both training and test data using the trained PCA
    train_pca = pcaobj.transform(train_scaled)
    test_pca = pcaobj.transform(test_scaled)
    
    # Create PCA result DataFrames
    pcs = [f"PC{i+1}" for i in range(n_components)]
    df_train_pca = pd.DataFrame(train_pca, index=df_train_transposed.index, columns=pcs)
    df_test_pca = pd.DataFrame(test_pca, index=df_test_transposed.index, columns=pcs)


    return pcaobj, df_train_pca, df_test_pca

def get_pca_exp_meta_merged(df_pca_exp, df_meta_pca_input):
    df_pca_exp = df_pca_exp.sort_index(axis=0)
    df_meta_pca_input = df_meta_pca_input.sort_index(axis=0)
    if list(df_pca_exp.index) == list(df_meta_pca_input.index):
        df_pca_exp_meta_merged = pd.concat([df_pca_exp, df_meta_pca_input], axis=1)
    else:
        print("Different Order of Sample Names")

    return df_pca_exp_meta_merged

# %%
import os
group = "viral_v1_prev"
WORKDIR = ".."
path_exp_pca_input = os.path.join(WORKDIR, "Expression", f"{group}_pca_input.txt")
path_meta_pca_input = os.path.join(WORKDIR, "Metadata", f"{group}_pca_input.txt")
df_exp_pca_input = pd.read_csv(path_exp_pca_input, sep="\t", index_col=[0])
df_meta_pca_input = pd.read_csv(path_meta_pca_input, sep="\t", index_col=[0])
pcaobj, df_pca_exp, _ = get_pcaobj_and_transformed_data(df_exp_pca_input, df_exp_pca_input, n_components)
df_pca_meta = get_pca_exp_meta_merged(df_pca_exp, df_meta_pca_input)
df_pca_meta["Sample_Project"] = list(map(lambda x: str(x)[0], list(df_pca_meta.index)))

# %%
pcx_num = 1
pcy_num = 2
pcx = f"PC{pcx_num}"
pcy = f"PC{pcy_num}"
nrow = math.ceil(len(list_meta_columns) / ncol)
fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=figsize)

if nrow > 1 or ncol > 1:
    axes = axes.flatten()
else:
    axes = np.array([axes]).flatten()

dict_legend_ncol = dict()
for i, (meta, ax) in enumerate(zip(list_meta_columns, axes)):
    fig_letter = list(string.ascii_uppercase)[i]
    val_groupby_meta_pcx = df_pca_meta.groupby(meta)[pcx].apply(np.array).tolist()
    val_groupby_meta_pcy = df_pca_meta.groupby(meta)[pcy].apply(np.array).tolist()
    cmap = get_cmap("RdBu", len(val_groupby_meta_pcx))
    list_colors = [cmap(ind) for ind in range(len(val_groupby_meta_pcx))]
    dim_level = 10
    list_colors_cmyk = [pca.rgb_to_cmyk(rgb[0], rgb[1], rgb[2]) for rgb in list_colors]
    list_colors_cmyk_dim = [(cmyk[0], cmyk[1], cmyk[2], cmyk[3] + dim_level) for cmyk in list_colors_cmyk]
    list_colors_rgb_dim = [pca.cmyk_to_rgb(cmyk[0], cmyk[1], cmyk[2], cmyk[3]) for cmyk in list_colors_cmyk_dim]
    list_legends = df_pca_meta.groupby(meta)[pcy].apply(np.array).index.tolist()
    list_legends = [int(float(x)) if str(x).endswith(".0") else x for x in list_legends]
    legend_ncol = dict_legend_ncol.get(meta, 1)

    for (x, y, color, legend) in zip(val_groupby_meta_pcx, val_groupby_meta_pcy, list_colors_rgb_dim, list_legends):
        ax.scatter(x, y, c=color, s=plot_linewidth*50, label=legend, alpha=0.3, zorder=-11) 
        pca.confidence_ellipse(x, y, ax, n_std=3.0, edgecolor=color, facecolor='none', linewidth=plot_linewidth*2.0, alpha=0.5)
        ax.scatter(np.mean(x), np.mean(y), c=color, edgecolor="k", linewidth=plot_linewidth*0.5, marker="*", s=plot_linewidth*500, alpha=1, zorder=999)  # Increased mean marker size

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
                    bbox_to_anchor=(1.5, 1.0), 
                    frameon=True)

    plt.setp(legend.get_texts(), weight='bold')

    for handle in legend.legendHandles:
        handle.set_alpha(1)
    
    ax.annotate(fig_letter,
                xy=(-0.15, 1.04),
                xycoords='axes fraction',
                xytext=(0, 0),
                textcoords='offset points',
                size=plt.rcParams["font.size"]+10, 
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
