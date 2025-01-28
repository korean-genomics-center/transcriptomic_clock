# %%
import os

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

import draw_enrichment_plot as enrich
from transcriptomic_clock import Clock

# %%
dict_name_conversion = {"GTEx": "GTEx (Whole Blood)"}
WORKDIR = ".."
cm = 1/2.54
path = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 9
width = 10*cm
height = 10*cm
figsize = (width, height)
plot_linewidth = 3
tick_spacing = 10

fig = plt.figure(figsize = figsize)
col = 100
row = 100
gsfig = gridspec.GridSpec(
                        row, col,
                        left = 0, 
                        right = 1, 
                        bottom = 0,
                        top = 1,
                        wspace = 1, 
                        hspace = 1)

gs1 = gsfig[0:100, 0:100]
ax1 = fig.add_subplot(gs1)
WORKDIR = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean"
dir_trained_model = f"{WORKDIR}/LASSO_INFO/healthy_cohort_train_prev_median_filter/corr_0.35_above_pval_0.05_below"
file_testing_data = f"{dir_trained_model}/{list(dict_name_conversion.keys())[0]}_std_scaled.txt"
os.makedirs(dir_trained_model, exist_ok=True)
clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
clock.make_dataframe_error_statistics_per_sample()
df_testing_dataset = clock.get_testing_dataset()
df_metadata = pd.read_csv("/BiO/Research/GeroExpressome/Resources/Data/GTEx/GTEx_Analysis_v8_Annotations_Metadata.txt", sep="\t")
df_metadata_select = df_metadata[["RACE", "DTHHRDY", "Sample_Trait"]]
df_merged = pd.concat([df_testing_dataset, df_metadata_select], axis=1)
real_age = clock.get_y()
real_age_reshape = np.array(real_age).reshape(-1, 1)
observed_transcriptomic_age = clock.get_predicted_y()
df_merged["RNAage"] = observed_transcriptomic_age
df_merged["TAA"] = df_merged["RNAage"] - df_merged["Sample_Age"]
hardy = 1
df_merged = df_merged[df_merged["DTHHRDY"].isin([hardy])]
corr, pval = spearmanr(df_merged["Sample_Age"], df_merged["RNAage"])
mae = np.mean(abs(df_merged["RNAage"] - df_merged["Sample_Age"]))
sns.scatterplot(data=df_merged, x="Sample_Age", y="RNAage", color="white", edgecolor="k", linewidth=1, s=20, alpha=0.7, zorder=2)
ax1.plot([min(real_age), max(real_age)], [min(real_age), max(real_age)], color="firebrick", linestyle="dotted", linewidth=1, label="Ground truth", zorder=5)
ax1.tick_params(axis='both', which='major', labelsize=plt.rcParams["font.size"]+1)
ax1.legend(loc='center', bbox_to_anchor=(1.2, 0.5), fontsize=plt.rcParams["font.size"]+1, frameon=False).set_zorder(30)
ax1.set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+2)
ax1.set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+2)
title = f"{dict_name_conversion[list(dict_name_conversion.keys())[0]]}, r={round(corr, 2)}, mae={round(mae, 2)}"
ax1.set_title(title, fontsize=plt.rcParams["font.size"]+1, weight="bold", pad=1.0)
ax1.set_xlim(0, 91)
ax1.set_ylim(0, )
ax1.grid(axis="both")
ax1.legend().set_visible(False)

# %%
# df_merged_nona = df_merged[df_merged["Sample_Trait"].astype(str) != "nan"]
# df_merged_cardiac_arrest = df_merged_nona[df_merged_nona["Sample_Trait"].str.lower().str.contains("myocardial")]
# df_merged_trauma = df_merged_nona[df_merged_nona[""].str.lower().str.contains("blunt")]