# %%
import os
import matplotlib.font_manager as fm
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from transcriptomic_clock import Clock

# %%
dict_name_conversion = {"GTEx": "GTEx (Whole Blood)"}
WORKDIR = ".."
cm = 1/2.54
path = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 9
width = 15*cm
height = 15*cm
figsize = (width, height)
plot_linewidth = 3
tick_spacing = 10

fig = plt.figure(figsize=(width, height))
gsfig = gridspec.GridSpec(2, 2, left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.2)

WORKDIR = "../"
dir_trained_model = f"{WORKDIR}/LASSO_INFO/healthy_cohort_train_prev_median_filter/corr_0.35_above_pval_0.05_below"
file_testing_data = f"{dir_trained_model}/{list(dict_name_conversion.keys())[0]}_std_scaled.txt"
os.makedirs(dir_trained_model, exist_ok=True)
clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
clock.make_dataframe_error_statistics_per_sample()
df_testing_dataset = clock.get_testing_dataset()
df_metadata = pd.read_csv("/BiO/Research/GeroExpressome/Resources/Data/GTEx/GTEx_Analysis_v8_Annotations_Metadata_Whole_Blood.txt", sep="\t")
df_metadata_select = df_metadata[["RACE", "DTHHRDY", "SMRIN", "Sample_Trait"]]
df_metadata_select["RNAQuality"] = df_metadata_select["SMRIN"].apply(lambda x: float())
df_merged = pd.concat([df_testing_dataset, df_metadata_select], axis=1)
real_age = clock.get_y()
real_age_reshape = np.array(real_age).reshape(-1, 1)
observed_transcriptomic_age = clock.get_predicted_y()
df_merged["RNAage"] = observed_transcriptomic_age
df_merged["TAA"] = df_merged["RNAage"] - df_merged["Sample_Age"]

corr, pval = spearmanr(df_merged["Sample_Age"], df_merged["RNAage"])
mae = np.mean(abs(df_merged["RNAage"] - df_merged["Sample_Age"]))
r2 = r2_score(df_merged["Sample_Age"], df_merged["RNAage"])

fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=True)

scatter_kwargs = {"data": df_merged, "x": "Sample_Age", "y": "RNAage", "color": "w", 
                  "edgecolor": "k", "linewidth": 0.5, "s": 20, "alpha": 0.7, "zorder": 2}

sns.scatterplot(ax=axes[0], **scatter_kwargs)
axes[0].plot([min(real_age), max(real_age)], [min(real_age), max(real_age)], 
             color="firebrick", linestyle="dotted", linewidth=1, label="Ground truth", zorder=5)
axes[0].set_ylim(0, 200)
axes[0].set_ylabel("Predicted Age (years)", fontsize=10)
axes[0].set_yticks(np.arange(0, 201, 20))

sns.scatterplot(ax=axes[1], **scatter_kwargs)
axes[1].plot([min(real_age), max(real_age)], [min(real_age), max(real_age)], 
             color="firebrick", linestyle="dotted", linewidth=1, label="Ground truth", zorder=5)
axes[1].set_ylim(0, 4000)
axes[1].set_ylabel("Predicted Age (years)", fontsize=10)
axes[1].set_yticks(np.arange(0, 4001, 500))

for i, ax in enumerate(axes):
    ax.set_xlabel("Chronological Age (years)", fontsize=10)
    ax.grid(axis="both", alpha=0.3)
    ax.set_xticks(np.arange(0, 81, 10))
    ax.tick_params(axis='both', which='major', labelsize=9)
    ax.text(2, ax.get_ylim()[1] * 0.9, f"Pearson's r:{round(corr, 2)}\nMAE:{round(mae, 2)}\n$R^2$:{round(r2, 2)}",
            fontsize=9)    
    ax.text(-0.25, 1.05, chr(65 + i), transform=ax.transAxes, fontsize=12, 
            fontweight='bold', va='top', ha='left')


plt.tight_layout()

plt.savefig(f"../Figures/Supplementary_Figure_11.tiff", dpi=300, bbox_inches="tight")
plt.savefig(f"../Figures/Supplementary_Figure_11.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# %%
