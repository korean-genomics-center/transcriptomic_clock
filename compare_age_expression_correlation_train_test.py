# %%
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from scipy.stats import pearsonr

# Load and prepare data
train_corr = "./LinearRegression/standardized/healthy_illumina/corr_0.0_above_pval_1.0_below/linear_regression_expression_Sample_Age_train.tsv"
valid_corr = "./LinearRegression/standardized/healthy_illumina/corr_0.0_above_pval_1.0_below/linear_regression_expression_Sample_Age_test.tsv"
test_corr = "./LinearRegression/standardized/healthy_illumina/corr_0.0_above_pval_1.0_below/linear_regression_expression_Sample_Age_bgi.tsv"

df_train_corr = pd.read_csv(train_corr, sep="\t")
df_valid_corr = pd.read_csv(valid_corr, sep="\t")
df_test_corr = pd.read_csv(test_corr, sep="\t")
df_train_corr = df_train_corr[["gene_id", "corr"]]
df_valid_corr = df_valid_corr[["gene_id", "corr"]]
df_test_corr = df_test_corr[["gene_id", "corr"]]
df_train_valid = pd.merge(df_train_corr, df_valid_corr, on="gene_id", how="inner", suffixes=["_Train", "_Valid"])
data = pd.merge(df_train_valid, df_test_corr, on="gene_id", how="inner")
data = data.rename(columns={"corr": "corr_Test"})

# %%
# Set font properties and figure size
cm = 1/2.54
path = "./Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
width = 16*cm
height = 10*cm
figsize = (width, height)
plot_linewidth = 0.3

# Create figure and GridSpec layout
fig = plt.figure(figsize=figsize)
col = 160
row = 100
gsfig = gridspec.GridSpec(row, col, left=0, right=1, bottom=0, top=1, wspace=2, hspace=2)

gs1 = gsfig[10:45, 10:50]
ax1 = fig.add_subplot(gs1)

gs2 = gsfig[10:45, 60:100]
ax2 = fig.add_subplot(gs2)

gs3 = gsfig[10:45, 110:150]
ax3 = fig.add_subplot(gs3)

gs4 = gsfig[55:90, 10:150]
ax4 = fig.add_subplot(gs4)

# Plot in the first subplot
sns.scatterplot(ax=ax1, data=data, x="corr_Train", y="corr_Valid", color="white", edgecolor="k")
sns.kdeplot(ax=ax1, data=data, x="corr_Train", y="corr_Valid", fill=True, cmap="inferno")
sns.lineplot(ax=ax1, x=np.arange(-0.7, 0.7, 0.1), y=np.arange(-0.7, 0.7, 0.1), color="grey", linestyle="dashed", linewidth=1)
corr, pval = pearsonr(data["corr_Train"], data["corr_Valid"])
ax1.set_xlabel("Age Correlations in train data", fontsize=plt.rcParams["font.size"]+1)
ax1.set_ylabel("Age Correlations in validation data", fontsize=plt.rcParams["font.size"]+1)
ax1.text(-0.65, 0.60, f"r = {round(corr, 2)}", weight="bold", fontsize=plt.rcParams["font.size"] + 1)
ax1.set_xlim(-0.7, 0.7)
ax1.set_ylim(-0.7, 0.7)

# Plot in the second subplot
sns.scatterplot(ax=ax2, data=data, x="corr_Train", y="corr_Test", color="white", edgecolor="k")
sns.kdeplot(ax=ax2, data=data, x="corr_Train", y="corr_Test", fill=True, cmap="inferno")
sns.lineplot(ax=ax2, x=np.arange(-0.7, 0.7, 0.1), y=np.arange(-0.7, 0.7, 0.1), color="grey", linestyle="dashed", linewidth=1)
corr, pval = pearsonr(data["corr_Train"], data["corr_Test"])
ax2.set_xlabel("Age Correlations in train data", fontsize=plt.rcParams["font.size"]+1)
ax2.set_ylabel("Age Correlations in test data", fontsize=plt.rcParams["font.size"]+1)
ax2.text(-0.65, 0.60, f"r = {round(corr, 2)}", weight="bold", fontsize=plt.rcParams["font.size"] + 1)
ax2.set_xlim(-0.7, 0.7)
ax2.set_ylim(-0.7, 0.7)

# Plot in the third subplot
sns.scatterplot(ax=ax3, data=data, x="corr_Valid", y="corr_Test", color="white", edgecolor="k")
sns.kdeplot(ax=ax3, data=data, x="corr_Valid", y="corr_Test", fill=True, cmap="inferno")
sns.lineplot(ax=ax3, x=np.arange(-0.7, 0.7, 0.1), y=np.arange(-0.7, 0.7, 0.1), color="grey", linestyle="dashed", linewidth=1)
corr, pval = pearsonr(data["corr_Valid"], data["corr_Test"])
ax3.set_xlabel("Age Correlations in validation data", fontsize=plt.rcParams["font.size"]+1)
ax3.set_ylabel("Age Correlations in test data", fontsize=plt.rcParams["font.size"]+1)
ax3.text(-0.65, 0.60, f"r = {round(corr, 2)}", weight="bold", fontsize=plt.rcParams["font.size"] + 1)
ax3.set_xlim(-0.7, 0.7)
ax3.set_ylim(-0.7, 0.7)

# Prepare data for stacked bar plot in the fourth subplot
path_selected_features = "./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below/regression_results.txt"
df_features = pd.read_csv(path_selected_features, sep="\t")
list_selected_features = df_features["Selected_Features"].to_list()
data_set_idx = data.set_index("gene_id")
data_selec = data_set_idx.loc[list_selected_features]
label = data_selec.index
label_short = list(map(lambda x: x.split("_")[-1], label))
x = np.arange(len(label))

# Plot stacked bar plot
width = 0.2
ax4.bar(x - width, abs(data_selec["corr_Train"]), width, label='Train', color="blue", edgecolor="black", linewidth=plot_linewidth)
ax4.bar(x, abs(data_selec["corr_Valid"]), width, label='Validation', color="orange", edgecolor="black", linewidth=plot_linewidth)
ax4.bar(x + width, abs(data_selec["corr_Test"]), width, label='Test', color="green", edgecolor="black", linewidth=plot_linewidth)

ax4.legend(fontsize=plt.rcParams["font.size"], title="Dataset", title_fontsize=plt.rcParams["font.size"]+1)
ax4.set_xticks(x)
ax4.set_xticklabels(label_short, rotation=45, rotation_mode="anchor", ha="right", fontstyle="italic", weight="bold", fontsize=plt.rcParams["font.size"])
ax4.set_yticklabels(ax4.get_yticklabels(minor=False), fontsize=plt.rcParams["font.size"])
ax4.set_xlabel("Gene Symbol", fontsize=plt.rcParams["font.size"]+1)
ax4.set_ylabel("Absolute Correlation", fontsize=plt.rcParams["font.size"]+1)
ax4.margins(0.01)

# Add figure letters
fig.text(0.01, 0.95, 'A', ha='center', va='center', fontsize=plt.rcParams["font.size"]+2, weight='bold')
fig.text(0.32, 0.95, 'B', ha='center', va='center', fontsize=plt.rcParams["font.size"]+2, weight='bold')
fig.text(0.63, 0.95, 'C', ha='center', va='center', fontsize=plt.rcParams["font.size"]+2, weight='bold')
fig.text(0.01, 0.50, 'D', ha='center', va='center', fontsize=plt.rcParams["font.size"]+2, weight='bold')

fig.tight_layout()
plt.savefig("./Figures/Supplementary_Figure_3.tiff", dpi=600, bbox_inches="tight")
plt.savefig("./Figures/Supplementary_Figure_3.png", dpi=600, bbox_inches="tight")
plt.show()
plt.close()
# %%
