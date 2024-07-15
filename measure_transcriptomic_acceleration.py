# %%
import glob
import os
from collections import Counter
from itertools import chain

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib.cm import get_cmap
from scipy.stats import f_oneway, sem, ttest_1samp
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import fdrcorrection

# %%
WORKDIR = "./LASSO_INFO/standardized/healthy_illumina"
dict_name_conversion = {"healthy": "Healthy Cohort1", 
                        "bgi": "Healthy Cohort2", 
                        "viral_v1_long": "Acute Phase",
                        "viral_v2_long": "Mid Phase",
                        "viral_v3_long": "Late Phase",
                        "viral_recov": "Convalescence",
                        "anxiety": "Anxiety", 
                        "attempt": "Suicide Attempt",
                        "depression": "Major Depression"}

list_files = glob.glob(f"{WORKDIR}/corr_0.3_above_pval_0.05_below/sample_statistics*", recursive=True)
list_files_filtered = list(filter(lambda x: os.path.basename(x).replace("sample_statistics_", "").replace(".tsv", "") in dict_name_conversion.keys(), list_files)) 

dict_groupname_residual = dict()
for file in list_files_filtered:
    groupname = os.path.basename(file).split("sample_statistics")[-1][1:].replace(".tsv", "")
    groupname_new = dict_name_conversion[groupname]
    df_residual = pd.read_csv(file, sep="\t")
    df_residual_copy = df_residual.copy()
    array_error = df_residual_copy["error"].to_numpy()
    dict_groupname_residual[groupname_new] = array_error

dict_groupname_residual = dict(sorted(dict_groupname_residual.items(), key=lambda x: list(dict_name_conversion.values()).index(x[0])))

# %%

# Extracting necessary values from the dictionary and ordering y-axis
list_cohorts = list(dict_groupname_residual.keys())
list_errors = list(dict_groupname_residual.values())
list_mean_error = list(map(np.mean, list_errors))
list_sem_error = list(map(sem, list_errors))
list_me_error = list(map(lambda x: 1.96 * x, list_sem_error))
list_lower_ci = list(map(lambda x: x[0] - x[1], zip(list_mean_error, list_me_error)))
list_upper_ci = list(map(lambda x: x[0] + x[1], zip(list_mean_error, list_me_error)))
list_sample_sizes = list(map(len, list_errors))
list_pval = [ttest_1samp(errors, 0)[1] for errors in list_errors]
_, list_pval_corrected = fdrcorrection(list_pval)

# Define groups (not changing the grouping here, assuming it's already defined as needed)
dict_groups = {
    "Healthy": ['Healthy Cohort1', 'Healthy Cohort2'],
    "Viral\nInfection": ['Acute Phase', 'Mid Phase', 'Late Phase', 'Convalescence'],
    "Psychiatric\nDisorders": ['Anxiety', 'Suicide Attempt', 'Major Depression']

}

# %%
import compare_severity_infection_per_age_group as crf
from scipy.stats import mannwhitneyu

path_sev_info = "/BiO/Access/kyungwhan1998/transcriptomic_aging_paper/Clinical/infectomics_CRF_20230410_edit.xlsx"
path_acc_v1 = "./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below/sample_statistics_viral_v1_long.tsv"
path_acc_v2 = "./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below/sample_statistics_viral_v2_long.tsv"
path_acc_v3 = "./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below/sample_statistics_viral_v3_long.tsv"

df_excel = pd.read_excel(path_sev_info, engine="openpyxl", skiprows=1)
df_crf = df_excel.rename(columns={"Subject NO.(고유번호)": "SUBJID"})
df_crf["Individual"] = df_crf["SUBJID"].apply(lambda x: "-".join(x.split("-")[:-1]))
df_crf["sample"] = df_crf["SUBJID"]
df_crf["CRP"] = df_crf["CRP"].apply(lambda x: float(str(x).replace("0..31", "0.31")))
df_crf_v1 = df_crf[df_crf["VISIT_no."]==1]
dict_v1_crp = dict(zip(df_crf_v1["Individual"], df_crf_v1["CRP"]))
dict_crp_group = crf.group_patients_by_crp_level(dict_v1_crp)
list_crp_group = crf.get_column_crp_group(df_crf, dict_crp_group)
df_crf["CRP_Group"] = list_crp_group
df_viral_v1 = pd.read_csv(path_acc_v1, sep="\t")
df_viral_v2 = pd.read_csv(path_acc_v2, sep="\t")
df_viral_v3 = pd.read_csv(path_acc_v3, sep="\t")
df_merge_v1 = pd.merge(df_crf, df_viral_v1, on="SUBJID", how="inner")
df_merge_v2 = pd.merge(df_crf, df_viral_v2, on="SUBJID", how="inner")
df_merge_v3 = pd.merge(df_crf, df_viral_v3, on="SUBJID", how="inner")
df_merge_all = pd.concat([df_merge_v1, df_merge_v2, df_merge_v3], axis=0
)

# %%
list_crp_groups = df_merge_v1["CRP_Group"].to_list()
dict_cnt = Counter(list_crp_groups)
from scipy.stats import ttest_ind

error_crp_groups = df_merge_v1.groupby("CRP_Group")["error"].apply(list)
list_high_crp_error = error_crp_groups[0]
list_low_crp_error = error_crp_groups[1]
# stat, pval = ttest_ind(list_high_crp_error, list_low_crp_error, equal_var=False, alternative="greater")
stat, pval = mannwhitneyu(list_high_crp_error, list_low_crp_error, alternative='greater')
# %%
cm = 1/2.54
path = "./Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
width = 10*cm
height = 10*cm
figsize = (width, height)
plot_linewidth = 1

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

gs1 = gsfig[0:20, 70:100]
ax1 = fig.add_subplot(gs1)
ax1.text(-0.3, 1.1, "b", transform=ax1.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')

gs2 = gsfig[30:50, 70:100]
ax2 = fig.add_subplot(gs2)
ax2.text(-0.3, 1.1, "c", transform=ax2.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')

gs3 = gsfig[60:80, 70:100] 
ax3 = fig.add_subplot(gs3)
ax3.text(-0.3, 1.1, "d", transform=ax3.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')

gs4 = gsfig[0:80, 0:30]
ax4 = fig.add_subplot(gs4)
ax4.text(-0.75, 1.02, "a", transform=ax4.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')

plt.rcParams["font.size"] = 5
sns.barplot(df_merge_v1, 
               x="CRP_Group", 
               y = "error", 
               hue="CRP_Group",
               err_kws = {"linewidth": 1},
               ax=ax1)
ax1.text(0.5, 140, s=f"P={round(float(pval), 5)}", ha="center", weight="bold")
ax1.set_yscale("log")
ax1.set_ylim(0, 310)
ax1.hlines(120, xmin=0, xmax=1, color="k", linewidth=0.5)
ax1.set_xlabel("CRP Levels", fontsize=plt.rcParams["font.size"]+1)
ax1.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+1)
ax1.set_xticklabels([f"{group}\n(N={count})"for group, count in dict_cnt.items()])
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# Custom legend
legend_elements = [
    plt.Line2D([0], [0], 
               marker='o',
               markeredgecolor='tab:orange' if 'Low' in group else 'tab:blue',
               lw=0,
               label=f'{group} (N={count})',
               markerfacecolor="None", 
               markersize=plt.rcParams["font.size"])
    for group, count in dict_cnt.items()]

crf.hue_regplot(df_merge_v1, 
           x="Neutrophil", 
           y="error",
           hue="CRP_Group", 
           scatter_kws={"s": plot_linewidth+11,
                        "marker":'o',
                        "facecolor":"None",
                        "linewidths":0.6},
           line_kws={"linewidth": 0},
           ci=0,
           ax=ax2)
ax2.set_xlabel("Neutrophil Count (%)", fontsize=plt.rcParams["font.size"]+1)
ax2.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+1)
ax2.set_xlim(20, 100)
ax2.set_ylim(-50, 250)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.legend(handles=legend_elements, 
           loc='upper left', 
           title='CRP Levels', 
           title_fontsize=plt.rcParams["font.size"]+1, 
           fontsize=plt.rcParams["font.size"],
           frameon=False)

crf.hue_regplot(df_merge_v1,  
           x="Lymphocytes", 
           y="error", 
           hue="CRP_Group", 
           scatter_kws={"s": plot_linewidth+11,
                        "marker":'o',
                        "facecolor":"None",
                        "linewidths":0.5},
           line_kws={"linewidth": 0},
           ci=0,
           ax=ax3)
ax3.set_xlabel("Lymphocyte Count (%)", fontsize=plt.rcParams["font.size"]+1)
ax3.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+1)
ax3.set_xlim(0, 80)
ax3.set_ylim(-50, 250)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.legend(handles=legend_elements, 
           loc='upper right', 
           title='CRP Levels', 
           title_fontsize=plt.rcParams["font.size"]+1, 
           fontsize=plt.rcParams["font.size"],
           frameon=False)

# Add vertical line at summarized result for Group 1
ax4.axvline(x=0, ymin=0, ymax=0.95, color='black', linestyle='solid', linewidth=plot_linewidth-0.4, zorder=-999)
# Set labels for y-axis with cohort names and sample sizes
ax4.set_yticks(range(len(list_cohorts)))
yticklabels = list(map(lambda x: f"{x[0]} ({x[1]})", zip(list_cohorts, list_sample_sizes)))
ax4.set_yticklabels(yticklabels)  # Ordering y-axis according to dict_groupname_residual
# Invert y-axis to have the first study at the top
ax4.invert_yaxis()
ax4.spines['right'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax4.spines['top'].set_visible(False)

# Add group labels on the leftmost side
for group, cohorts in dict_groups.items():
    first_cohort_index = list_cohorts.index(cohorts[0])
    group_label_y = (first_cohort_index + list_cohorts.index(cohorts[-1])) / 2 # Position label in the middle of the group
    ax4.text(-55, group_label_y, group, va='center', ha='center', color='black', fontsize=plt.rcParams["font.size"], weight="bold", rotation=90)
   
line_adjust = 0.2
ax4.vlines(x=-50, ymin=0-line_adjust, ymax=1+line_adjust, color='black', linewidth=plot_linewidth-0.5, clip_on = False)
ax4.vlines(x=-50, ymin=2-line_adjust, ymax=5+line_adjust, color='black', linewidth=plot_linewidth-0.5, clip_on = False)
ax4.vlines(x=-50, ymin=6-line_adjust, ymax=8+line_adjust, color='black', linewidth=plot_linewidth-0.5, clip_on = False)

for i, cohort in enumerate(list_cohorts):
    if (list_pval_corrected[i] < 0.05):
        markercolor = 'tab:red'
        textcolor = 'tab:red'
        weight = 'bold'
    else:
        markercolor = 'grey'
        textcolor = "k"
        weight = 'normal'
    # Plot the mean effects
    ax4.plot(list_mean_error[i], i, 'd', color=markercolor, markersize=plot_linewidth+2.5, zorder=999, markeredgewidth=plot_linewidth-0.7, markeredgecolor=markercolor)
    ax4.hlines(i, list_lower_ci[i], list_upper_ci[i], color=markercolor, linewidth=plot_linewidth)
    ax4.text(1.2, i, f'{list_mean_error[i]:.2f} [{list_lower_ci[i]:.2f}, {list_upper_ci[i]:.2f}]', va='center', ha='center', color=textcolor, weight=weight, transform=ax4.get_yaxis_transform())
    ax4.text(1.7, i, f'{list_pval_corrected[i]:.2f}', va='center', ha='center', color=textcolor, weight=weight, transform=ax4.get_yaxis_transform())

ax4.text(-0.3, 1.0, "Cohorts (N)", ha='center', va='center', transform=ax4.transAxes, fontsize=plt.rcParams["font.size"], weight="bold")
ax4.text(1.2, 1.0, "Mean Acc. [95% CI]", ha='center', va='center', transform=ax4.transAxes, fontsize=plt.rcParams["font.size"], weight="bold")
ax4.text(1.7, 1.0, "FDR", ha='center', va='center', transform=ax4.transAxes, fontsize=plt.rcParams["font.size"], weight="bold")

# Add labels for x-axis and title
ax4.set_xlabel('TAA (years)', fontsize=plt.rcParams["font.size"]+1)
ax4.set_xlim(-11, 46)
ax4.tick_params(axis='y', which='both',length=0)

plt.subplots_adjust(left=0.8)

# Display the plot
plt.tight_layout()
plt.savefig("Figures/Fig_2.png", dpi=300, bbox_inches="tight")
plt.savefig("Figures/Fig_2.pdf", dpi=300, bbox_inches="tight")
plt.savefig("Figures/Fig_2.jpeg", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

 # %% [Prediction Error]

cm = 1/2.54
path = "./Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
width = 12*cm
height = 10*cm

# Determine the number of cohorts and subplots
num_cohorts = len(list_cohorts)
N = int(np.ceil(np.sqrt(num_cohorts)))  # Adjust N dynamically based on num_cohorts

# Create a figure and subplots
fig, axs = plt.subplots(N, N, figsize=(width, height), constrained_layout=True)

# Iterate through cohorts and plot each histogram
for i, cohort in enumerate(list_cohorts):
    error = dict_groupname_residual[cohort]
    mean_error = round(np.mean(error), 2)
    
    # Determine subplot indices
    ax = axs[i // N, i % N]
    
    # Plot histogram
    ax.hist(error, color="grey", alpha=0.5)
    ax.axvline(x=mean_error, color="tab:red", linestyle="dotted")
    ax.text(mean_error+1, 1, f"Mean Error\n={mean_error}", color="tab:red", weight="bold", fontsize=plt.rcParams["font.size"]+1)
    ax.set_xlabel("Prediction Error (years)", fontsize=plt.rcParams["font.size"])
    ax.set_ylabel("Proportion (%)", fontsize=plt.rcParams["font.size"])
    ax.set_title(cohort, fontsize=plt.rcParams["font.size"]+1, weight="bold", pad=0.1)
    ax.set_xlim(int(min(error))-2, int(max(error))+2)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    # Label subplot with letters a-z
    subplot_letter = chr(ord('a') + i)  # Convert integer index to character
    ax.text(-0.5, 1.2, subplot_letter, transform=ax.transAxes,
            fontsize=plt.rcParams["font.size"]+1, fontweight='bold', va='top', ha='right')

# Hide unused subplots if any
for j in range(num_cohorts, N * N):
    axs.flatten()[j].axis('off')

# Adjust layout and show plot
plt.savefig("Figures/Extended_Data_Fig_7.png", dpi=300, bbox_inches="tight")
plt.savefig("Figures/Extended_Data_Fig_7.pdf", dpi=300, bbox_inches="tight")
plt.savefig("Figures/Extended_Data_Fig_7.jpeg", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# %%
# Calculate summarized result for different groups
list_groups = list(dict_groups.keys())

def get_overall_mean_ci_group_cohorts(groupname):
    if groupname == "Viral\nInfection":
        group_indices = [list_cohorts.index(cohort) for cohort in dict_groups[groupname]][:-1]
    else:
        group_indices = [list_cohorts.index(cohort) for cohort in dict_groups[groupname]]
    group_mean = np.mean([list_mean_error[idx] for idx in group_indices])
    group_lower_ci = np.min([list_lower_ci[idx] for idx in group_indices])
    group_upper_ci = np.max([list_upper_ci[idx] for idx in group_indices])
    group_pval = ttest_1samp(list(chain(*[list_errors[idx] for idx in group_indices])), 0)[1] 

    return [group_mean, group_lower_ci, group_upper_ci, group_pval]

dict_group_stats = dict()
for group in list_groups:
    list_stats = get_overall_mean_ci_group_cohorts(group)
    group_name_fix = group.replace("\n", "")
    dict_group_stats[group_name_fix] = list_stats

df_group_stats = pd.DataFrame.from_dict(dict_group_stats, orient="index").reset_index(drop=False)
df_group_stats.columns = ["group", "mean", "lower_ci", "upper_ci", "pval"]
list_pval = df_group_stats["pval"].to_list()
_, list_pval_corrected = fdrcorrection(list_pval)
df_group_stats["padj"] = list_pval_corrected
df_group_stats.to_csv("./Supplementary_Tables/Table_S4.txt", sep="\t", index=False)

import math

# %%
from scipy.stats import pearsonr
from statsmodels.stats.multitest import fdrcorrection

cm = 1/2.54
path = "./Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 10
width = 20*cm
height = 4*cm
figsize = (width, height)
plot_linewidth = 1

# Clean and identify numerical columns (excluding 'error' and 'VISIT_no.')
df_merge_all.columns = list(map(lambda x: x.replace("진단시", "").replace("값", "").replace("  ", " ").replace("\n", "").replace("ORF", " ORF").replace("L D H", "LDH").replace(" CT", "CT").replace("error", "TAA (years)"), df_merge_all.columns))
numerical_columns = df_merge_all.select_dtypes(include=['float64', 'int64']).columns
numerical_columns = list(numerical_columns.drop(['error', 'VISIT_no.'], errors='ignore'))[20: -35]
numerical_columns += ["CRP", "TAA (years)"]
numerical_columns.remove("CPK")
numerical_columns.remove("FDP")
numerical_columns.remove("Troponin I")
numerical_columns.remove("D-dimer")
numerical_columns.remove("Monocytes")
numerical_columns = sorted(numerical_columns)
df_merge_all_filt = df_merge_all[numerical_columns]

list_corrs = []
list_pvals = []
for column in numerical_columns:
    df_merge_all_na_filt = df_merge_all_filt[df_merge_all_filt[column].notna()]
    corr, pval = pearsonr(df_merge_all_na_filt["TAA (years)"], df_merge_all_na_filt[column])
    list_corrs.append(corr)
    list_pvals.append(pval)

_, list_padjs = fdrcorrection(list_pvals)

n_rows = 7
n_cols = 5
fig, axes = plt.subplots(n_rows, n_cols, figsize=(width, height * n_rows), sharey=True)
axes = axes.flatten()

for ax, column, corr, padj in zip(axes, numerical_columns, list_corrs, list_padjs):
    color = "tab:red" if padj < 0.05 else "grey"
    sns.regplot(
        df_merge_all_filt,
        x=column,
        y="TAA (years)",
        color=color,
        scatter_kws={"s": 10, "marker": 'o', "facecolor": "None", "linewidths": 0.5},
        line_kws={"linewidth": 1},
        ci=95,
        ax=ax
    )
    weight = "bold" if padj < 0.05 else "regular"
    ax.set_xlabel(column, fontsize=plt.rcParams["font.size"] + 1, color=color, weight=weight)
    ax.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"] + 1)
    # Position the text in the upper left corner within the plot area
    ax.text(0.05, 0.95, f"r = {round(corr, 3)}\nFDR = {round(padj, 3)}", transform=ax.transAxes,
            verticalalignment='top', fontsize=plt.rcParams["font.size"], color=color)

# Hide unused subplots if any
for ax in axes[len(numerical_columns):]:
    ax.axis('off')

plt.tight_layout()
plt.savefig("./Figures/Extended_Data_Fig_9.pdf", dpi=300, bbox_inches="tight")
plt.savefig("./Figures/Extended_Data_Fig_9.png", dpi=300, bbox_inches="tight")
plt.savefig("./Figures/Extended_Data_Fig_9.jpeg", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# %%
# cm = 1/2.54
# path = "./Arial.ttf"
# prop = fm.FontProperties(fname=path)
# plt.rcParams['font.family'] = prop.get_name()
# plt.rcParams["font.size"] = 5
# width = 5*cm
# height = 12*cm
# figsize = (width, height)
# plot_linewidth = 1
# list_visit_groups = df_merge_all["VISIT_no."].to_list()
# dict_cnt = Counter(list_visit_groups)
# phase_map = {1: "Acute", 2: "Mid", 3: "Late"}
# color_map = {1: "tab:blue", 2: "tab:orange", 3: "tab:green"}

# legend_elements = [
#     plt.Line2D([0], [0], 
#                marker='o',
#                markeredgecolor=color_map[group],
#                lw=0,
#                label=f'{phase_map[group]} (N={count})',
#                markerfacecolor="None", 
#                markersize=plt.rcParams["font.size"])
#     for group, count in dict_cnt.items()]

# fig, axes = plt.subplots(3, 1, figsize=figsize)
# axes = axes.flatten()

# crf.hue_regplot(df_merge_all, 
#         x="CRP", 
#         y="TAA (years)",
#         hue="VISIT_no.", 
#         scatter_kws={"s": 10,
#                         "marker":'o',
#                         "facecolor":"None",
#                         "linewidths":0.5},
#         line_kws={"linewidth": 1},
#         ci=95,
#         ax=axes[0])
# axes[0].set_xlabel("CRP (mg/dL)", fontsize=plt.rcParams["font.size"]+1)
# axes[0].set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+1)
# axes[0].text(-0.3, 1.1, "a", transform=axes[0].transAxes,
#         fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')
# axes[0].legend(handles=legend_elements, 
#            loc='upper left', 
#            title='Infect. Phase', 
#            title_fontsize=plt.rcParams["font.size"]+1, 
#            fontsize=plt.rcParams["font.size"],
#            frameon=False)


# crf.hue_regplot(df_merge_all, 
#         x="Neutrophil", 
#         y="TAA (years)",
#         hue="VISIT_no.", 
#         scatter_kws={"s": 10,
#                         "marker":'o',
#                         "facecolor":"None",
#                         "linewidths":0.5},
#         line_kws={"linewidth": 1},
#         ci=95,
#         ax=axes[1])
# axes[1].set_xlabel("Neutrophil Count (%)", fontsize=plt.rcParams["font.size"]+1)
# axes[1].set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+1)
# axes[1].text(-0.3, 1.1, "b", transform=axes[1].transAxes,
#         fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')
# axes[1].legend(handles=legend_elements, 
#            loc='upper left', 
#            title='Infect. Phase',  
#            title_fontsize=plt.rcParams["font.size"]+1, 
#            fontsize=plt.rcParams["font.size"],
#            frameon=False)


# crf.hue_regplot(df_merge_all, 
#         x="Lymphocytes", 
#         y="TAA (years)",
#         hue="VISIT_no.", 
#         scatter_kws={"s": 10,
#                         "marker":'o',
#                         "facecolor":"None",
#                         "linewidths":0.5},
#         line_kws={"linewidth": 1},
#         ci=95,
#         ax=axes[2])
# axes[2].set_xlabel("Lymphocyte Count (%)", fontsize=plt.rcParams["font.size"]+1)
# axes[2].set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+1)
# axes[2].text(-0.3, 1.1, "c", transform=axes[2].transAxes,
#         fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')
# axes[2].legend(handles=legend_elements, 
#            loc='upper right', 
#            title='Infect. Phase', 
#            title_fontsize=plt.rcParams["font.size"]+1, 
#            fontsize=plt.rcParams["font.size"],
#            frameon=False)


# plt.tight_layout()
# plt.savefig("./Figures/Extended_Data_Fig_11.pdf", dpi=300, bbox_inches="tight")
# plt.savefig("./Figures/Extended_Data_Fig_11.png", dpi=300, bbox_inches="tight")
# plt.show()
# plt.close()

# %%
