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
from scipy.stats import pearsonr

import draw_enrichment_plot as enrich
from transcriptomic_clock import Clock

# %%
WORKDIR = ".."
cm = 1/2.54
path = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 8
width = 10*cm
height = 25*cm
figsize = (width, height)
plot_linewidth = 3
tick_spacing = 10
dict_name_conversion = {"GSE273149": "COVID-19 ARDS Cohort (GSE273149)",
                        "GSE119117": "HCV Cohort (GSE119117)"}

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

list_test_group = list(dict_name_conversion.keys())

gs1 = gsfig[0:45, 0:100]
ax1 = fig.add_subplot(gs1)
idx = 0
hue = "time"
style = "disease status"
dict_conversion = {"time": "Stages", "disease status": "Group"}
dir_trained_model = f"{WORKDIR}/LASSO_INFO/healthy_cohort_train_prev_median_filter/corr_0.35_above_pval_0.05_below"
file_testing_data = f"{dir_trained_model}/{list_test_group[idx]}_std_scaled.txt"
os.makedirs(dir_trained_model, exist_ok=True)
clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
df_testing_dataset = clock.get_testing_dataset()
path_metadata = f"{WORKDIR}/Metadata/{list_test_group[idx]}.txt"
df_metadata = pd.read_csv(path_metadata, sep="\t")
df_metadata_selected = df_metadata[["Project-ID", hue, style]]
df_merged = pd.concat([df_testing_dataset, df_metadata_selected], axis=1)
list_samples_ards_con = list(df_merged[df_merged["disease status"] == "ARDS_CON"]["Sample_Age"].unique())
df_merged = df_merged[df_merged["disease status"] != "ARDS_CON"]

real_age = clock.get_y()[16:]
real_age_reshape = np.array(real_age).reshape(-1, 1)
observed_transcriptomic_age = clock.get_predicted_y()[16:]
df_merged["RNAage"] = observed_transcriptomic_age
df_merged["TAA"] = df_merged["RNAage"] - df_merged["Sample_Age"]
df_merged = df_merged.rename(columns=dict_conversion)
corr, pval = pearsonr(df_merged["Sample_Age"], df_merged["RNAage"])
mae = np.mean(abs(df_merged["RNAage"] - df_merged["Sample_Age"]))
sns.scatterplot(data=df_merged, x="Sample_Age", y="RNAage", color = "white", edgecolor="k", s=70, alpha=0.5, zorder=2, hue=dict_conversion.get(hue), style=dict_conversion.get(style), palette="RdBu")
ax1.text(-0.2, 1, s="A", transform=ax1.transAxes, fontsize=plt.rcParams["font.size"]+10, weight="bold")
ax1.plot([9, 91], [9, 91], color="firebrick", linestyle="dotted", linewidth=1, label="Ground truth", zorder=5)
ax1.tick_params(axis='both', which='major', labelsize=plt.rcParams["font.size"]+3)
ax1.legend(loc='center', bbox_to_anchor=(1.3, 0.5), fontsize=plt.rcParams["font.size"]+3, frameon=False).set_zorder(30)
ax1.set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+5)
ax1.set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+5)
title = f"{dict_name_conversion[list_test_group[idx]]}"
ax1.set_title(title, fontsize=plt.rcParams["font.size"]+5, weight="bold", pad=1.0)
ax1.grid(axis="both")
ax1.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
ax1.set_xlim(9, 71)
ax1.set_ylim(9, 171)

gs2 = gsfig[55:100, 0:100]
ax2 = fig.add_subplot(gs2)
idx = 1
hue = "time point"
style = "hcvgroup"
dict_conversion = {"time point": "Stages", "hcvgroup": "Group"}
dir_trained_model = f"{WORKDIR}/LASSO_INFO/healthy_cohort_train_prev_median_filter/corr_0.35_above_pval_0.05_below"
file_testing_data = f"{dir_trained_model}/{list_test_group[idx]}_std_scaled.txt"
os.makedirs(dir_trained_model, exist_ok=True)
clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
df_testing_dataset = clock.get_testing_dataset()
path_metadata = f"{WORKDIR}/Metadata/{list_test_group[idx]}.txt"
df_metadata = pd.read_csv(path_metadata, sep="\t")
df_metadata_selected = df_metadata[["Project-ID", hue, style]]
df_merged = pd.concat([df_testing_dataset, df_metadata_selected], axis=1)

real_age = clock.get_y()
real_age_reshape = np.array(real_age).reshape(-1, 1)
observed_transcriptomic_age = clock.get_predicted_y()
df_merged["RNAage"] = observed_transcriptomic_age
df_merged["TAA"] = df_merged["RNAage"] - df_merged["Sample_Age"]
df_merged = df_merged.rename(columns=dict_conversion)
corr, pval = pearsonr(df_merged["Sample_Age"], df_merged["RNAage"])
mae = np.mean(abs(df_merged["RNAage"] - df_merged["Sample_Age"]))
sns.scatterplot(data=df_merged, x="Sample_Age", y="RNAage", color = "white", edgecolor="k", s=70, alpha=0.5, zorder=2, hue=dict_conversion.get(hue), style=dict_conversion.get(style), palette="RdBu")
ax2.text(-0.2, 1, s="B", transform=ax2.transAxes, fontsize=plt.rcParams["font.size"]+10, weight="bold")
ax2.plot([-10, 91], [-10, 91], color="firebrick", linestyle="dotted", linewidth=1, label="Ground truth", zorder=5)
ax2.tick_params(axis='both', which='major', labelsize=plt.rcParams["font.size"]+3)
ax2.legend(loc='center', bbox_to_anchor=(1.25, 0.5), fontsize=plt.rcParams["font.size"]+3, frameon=False).set_zorder(30)
ax2.set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+5)
ax2.set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+5)
title = f"{dict_name_conversion[list_test_group[idx]]}"
ax2.set_title(title, fontsize=plt.rcParams["font.size"]+5, weight="bold", pad=1.0)
ax2.grid(axis="both")
ax2.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
ax2.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

ax2.set_xlim(19, 51)
ax2.set_ylim(-11, 111)

plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_10.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_10.tiff", dpi=300, bbox_inches="tight")
plt.show()
plt.close()