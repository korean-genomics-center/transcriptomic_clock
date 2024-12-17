# %%
import glob
import json
import os
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import sem, ttest_1samp, ttest_ind, ranksums
from transcriptomic_clock import Clock
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns
from matplotlib import gridspec

# %%
def get_dictionary_groupname_residual(dir_trained_model, dict_name_conversion, tag_file_sample="sample_statistics"):

    list_files = glob.glob(os.path.join(dir_trained_model, f"{tag_file_sample}*"), recursive=True)
    list_files_filtered = list(filter(lambda x: os.path.basename(x).replace(f"{tag_file_sample}_", "").replace("_std_scaled.txt", "") in dict_name_conversion.keys(), list_files)) 
    dict_groupname_residual_ards = dict()
    for file in list_files_filtered:
        groupname = os.path.basename(file).split(tag_file_sample)[-1][1:].replace("_std_scaled.txt", "")
        groupname_new = dict_name_conversion[groupname]
        df_residual = pd.read_csv(file, sep="\t")
        df_residual_copy = df_residual.copy()
        array_error = df_residual_copy["error"].to_numpy()
        dict_groupname_residual_ards[groupname_new] = array_error

    dict_groupname_residual_ards = dict(sorted(dict_groupname_residual_ards.items(), key=lambda x: list(dict_name_conversion.values()).index(x[0])))

    return dict_groupname_residual_ards

def get_sample_statistics(dir_trained_model, dict_name_conversion):
    dict_groupname_residual = get_dictionary_groupname_residual(dir_trained_model, dict_name_conversion)
    list_cohorts = list(dict_groupname_residual.keys())
    list_errors = list(dict_groupname_residual.values())
    list_mean_error = list(map(np.mean, list_errors))
    list_sem_error = list(map(sem, list_errors))
    list_me_error_ = list(map(lambda x: 1.96 * x, list_sem_error))
    list_lower_ci = list(map(lambda x: x[0] - x[1], zip(list_mean_error, list_me_error_)))
    list_upper_ci = list(map(lambda x: x[0] + x[1], zip(list_mean_error, list_me_error_)))
    list_sample_size = list(map(len, list_errors))
    list_pval = [ttest_1samp(errors, 0)[1] for errors in list_errors]
    # _, list_pval_corrected = fdrcorrection(list_pval)

    return list_cohorts, list_mean_error, list_lower_ci, list_upper_ci, list_sample_size, list_pval

def get_dict_sample_cohort_sample_size(list_cohorts, list_sample_size):
    dict_cohort_sample_size = dict(zip(list_cohorts, list_sample_size))

    return dict_cohort_sample_size

# %%
WORKDIR = ".."
train_group = "healthy_cohort_train_prev_median_filter"
dir_trained_model = f"{WORKDIR}/LASSO_INFO/{train_group}/corr_0.35_above_pval_0.05_below"

# %%
path_dict_name_conversion_json = f"{WORKDIR}/Data/dictionary_name_conversion_hcv.json"
with open (path_dict_name_conversion_json, mode="rb") as fb1:
    dict_name_conversion = json.load(fb1)

path_dict_cohort_information_json = f"{WORKDIR}/Data/dictionary_cohort_information_hcv.json"
with open(path_dict_cohort_information_json, mode="rb") as fb2:
    dict_groups = json.load(fb2)

list_test_group = list(dict_name_conversion.keys())
for i, group in enumerate(list_test_group):
    file_testing_data = f"{dir_trained_model}/{group}_std_scaled.txt"
    os.makedirs(dir_trained_model, exist_ok=True)
    clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
    df_clock = clock.make_dataframe_error_statistics_per_sample()
    df_clock.to_csv(os.path.join(dir_trained_model, f"sample_statistics_{group}_std_scaled.txt"), sep="\t", index=False)
    print(os.path.join(dir_trained_model, f"sample_statistics_{group}_std_scaled.txt"))

path_dict_name_conversion_ards_json = f"{WORKDIR}/Data/dictionary_name_conversion_ards.json"
with open (path_dict_name_conversion_ards_json, mode="rb") as fb1:
    dict_name_conversion_ards = json.load(fb1)

path_dict_cohort_information_ards_json = f"{WORKDIR}/Data/dictionary_cohort_information_ards.json"
with open(path_dict_cohort_information_ards_json, mode="rb") as fb2:
    dict_groups_ards = json.load(fb2)

list_group_ards = list(dict_name_conversion_ards.keys())

for i, group in enumerate(list_group_ards):
    file_testing_data = f"{dir_trained_model}/{group}_std_scaled.txt"
    os.makedirs(dir_trained_model, exist_ok=True)
    clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
    df_clock = clock.make_dataframe_error_statistics_per_sample()
    df_clock.to_csv(os.path.join(dir_trained_model, f"sample_statistics_{group}_std_scaled.txt"), sep="\t", index=False)
    print(os.path.join(dir_trained_model, f"sample_statistics_{group}_std_scaled.txt"))

# %%
cm = 1/2.54
path_font = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path_font)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 7
width = 16*cm
height = 14*cm
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

list_cohorts, list_mean_error, list_lower_ci, list_upper_ci, list_sample_size, list_pval = get_sample_statistics(dir_trained_model, dict_name_conversion)

dict_cohort_sample_size = dict(zip(list_cohorts, list_sample_size))

gs3 = gsfig[60:100, 0:30]
ax3 = fig.add_subplot(gs3)
ax3.text(-1.05, 1.05, "C", transform=ax3.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')

ax3.axvline(x=0, ymin=0, ymax=0.95, color='black', linestyle='solid', linewidth=plot_linewidth-0.4, zorder=-999)
ax3.set_yticks(range(len(list_cohorts)))
yticklabels = list(map(lambda x: f"{x[0]} ({x[1]})", zip(list_cohorts, list_sample_size)))
ax3.set_yticklabels(yticklabels)
ax3.invert_yaxis()
ax3.spines['right'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.spines['top'].set_visible(False)

ax3.text(-1.02, 0.5, "GSE119117", va='center', ha='center', color='black', fontsize=plt.rcParams["font.size"]+1, weight="bold", rotation=90, transform=ax3.transAxes)
ax3.vlines(x=-0.97, ymin=0.02, ymax=0.98, color='black', linewidth=plot_linewidth-0.5, clip_on = False, transform=ax3.transAxes)

for i, cohort in enumerate(list_cohorts):
    if (list_pval[i] < 0.05):
        markercolor = 'tab:red'
        textcolor = 'tab:red'
        weight = 'bold'
    else:
        markercolor = 'grey'
        textcolor = 'k'
        weight = 'normal'

    ax3.plot(list_mean_error[i], i, 'd', color=markercolor, markersize=plot_linewidth+4, zorder=999, markeredgewidth=plot_linewidth-0.7, markeredgecolor=markercolor)
    ax3.hlines(i, list_lower_ci[i], list_upper_ci[i], color=markercolor, linewidth=plot_linewidth)
    ax3.text(1.3, i, f'{list_mean_error[i]:.2f} [{list_lower_ci[i]:.2f}, {list_upper_ci[i]:.2f}]', va='center', ha='center', color=textcolor, weight=weight, transform=ax3.get_yaxis_transform())
    ax3.text(1.7, i, f'{list_pval[i]:.3f}', va='center', ha='center', color=textcolor, weight=weight, transform=ax3.get_yaxis_transform())

ax3.text(-0.4, 1.02, "Cohorts (N)", ha='center', va='center', transform=ax3.transAxes, fontsize=plt.rcParams["font.size"]+1, weight="bold")
ax3.text(1.3, 1.02, "Mean TAA [95% CI]", ha='center', va='center', transform=ax3.transAxes, fontsize=plt.rcParams["font.size"]+1, weight="bold")
ax3.text(1.7, 1.02, "P", ha='center', va='center', transform=ax3.transAxes, fontsize=plt.rcParams["font.size"]+1, weight="bold")

ax3.set_xlabel('TAA (years)', fontsize=plt.rcParams["font.size"]+1)
ax3.set_xlim(-20, 30)
ax3.tick_params(axis='y', which='both', length=0)

gs4 = gsfig[60:100, 65:100]
ax4 = fig.add_subplot(gs4)
ax4.text(-0.25, 1.05, "D", transform=ax4.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')
ax4.tick_params(axis='y', which='both', length=0)

df_error_GSE119117 = pd.read_csv(f"{WORKDIR}/Data/GSE119117_prediction_result_KoreanBloodClock.txt", sep="\t", index_col=[0])
df_meta_GSE119117 = pd.read_csv(f"{WORKDIR}/Metadata/GSE119117.txt", sep="\t", index_col=[0])
df_GSE119117 = pd.concat([df_error_GSE119117, df_meta_GSE119117], axis=1).reset_index(drop=False)
df_GSE119117['Phase'] = pd.Categorical(
    df_GSE119117['Phase'],
    categories=["Pre-infection", "Early acute", "Late acute", "Follow up"],
    ordered=True
)

# Count samples with "High" and "Low" CRP levels across stages
df_counts = df_GSE119117.groupby(["Phase", "hcvgroup"]).size().unstack(fill_value=0)

# Statistical comparison between "High" and "Low" groups at each stage
stage_pvals = {}
for stage in df_counts.index:
    resolution = df_GSE119117[(df_GSE119117["Phase"] == stage) & (df_GSE119117["hcvgroup"] == "Resolution")]["AgeAccel"]
    chronic = df_GSE119117[(df_GSE119117["Phase"] == stage) & (df_GSE119117["hcvgroup"] == "Chronic")]["AgeAccel"]
    
    if len(resolution) > 1 and len(chronic) > 1:  # Ensure sufficient samples
        pval = ttest_ind(resolution, chronic, equal_var=False).pvalue
    else:
        pval = None  # Not enough samples for comparison
    stage_pvals[stage] = pval

df_GSE119117["Subject-ID"] = df_GSE119117["Project-ID"].str.split("-").str[0]
dict_hcv_group = dict(zip(df_GSE119117["Subject-ID"], df_GSE119117["hcvgroup"]))

list_subject = df_GSE119117["Subject-ID"].to_list()
dict_subj_hcv_color = dict()
for subj in list_subject:
    ards_group = dict_hcv_group.get(subj, "Unknown")
    if ards_group == "Resolution":
        color="tab:green"
    elif ards_group == "Chronic":
        color="tab:orange"
    else:
        color="tab:blue"
    
    dict_subj_hcv_color[subj] = color

sns.pointplot(
    data=df_GSE119117,
    x="Phase",
    y="AgeAccel",
    marker="o",
    hue="Group",
    palette={"Resolution":"tab:green", "Chronic":"tab:orange"},
    markersize=5,
    linewidth=2,
    capsize=0.05,
    errorbar="se",
    zorder=3,
    ax=ax4
)

list_subject = df_GSE119117["Subject-ID"].to_list()
for subj in list_subject:
    df_GSE119117_copy = df_GSE119117.copy()
    subj_data = df_GSE119117_copy[df_GSE119117_copy["Subject-ID"] == subj]
    trend_color = dict_subj_hcv_color.get(subj, "tab:blue")
    sns.lineplot(
        data=subj_data,
        x="Phase",
        y="AgeAccel",
        marker="o",
        markersize=1,
        markerfacecolor=trend_color,
        markeredgecolor=None,
        color=trend_color,
        linewidth=1,
        alpha=0.1,
        ax=ax4
    )

ax4.set_ylim(-30, 80)
ax4.set_xlabel("Stages", fontsize=plt.rcParams["font.size"]+2, labelpad=20)
ax4.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+2)
ax4.margins(0.1)
ax4.legend(frameon=False)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
list_pval = list()
for idx, stage in enumerate(df_counts.index):
    resolution_cnt = df_counts.loc[stage, "Resolution"]
    chronic_cnt = df_counts.loc[stage, "Chronic"]
    ax4.text(idx, -0.13, f"Res. ({resolution_cnt})", ha="center", va="bottom", color="black", transform=ax4.get_xaxis_transform())
    ax4.text(idx, -0.18, f"Chron. ({chronic_cnt})", ha="center", va="bottom", color="black", transform=ax4.get_xaxis_transform())
    
    # Add p-value annotation
    pval = stage_pvals[stage]
    if pval is not None:
        list_pval.append(pval)

for idx, (stage, pval) in enumerate(zip(df_counts.index, list_pval)):
    if pval is not None and pval <= 0.05:
        ax4.text(idx, 0, f"{pval:.2f}", ha="center", va="bottom", fontsize=10, transform=ax4.get_xaxis_transform(), color="red", weight="bold")
    elif pval is not None and pval > 0.05:
        ax4.text(idx, 0.01, f"{pval:.2f}", ha="center", va="bottom", fontsize=8, transform=ax4.get_xaxis_transform(), color="black", weight="bold")
    else:
        continue

dict_groupname_residual_ards = get_dictionary_groupname_residual(dir_trained_model, dict_name_conversion_ards)
list_cohorts_ards = list(dict_groupname_residual_ards.keys())
list_errors_ards = list(dict_groupname_residual_ards.values())
list_mean_error_ards = list(map(np.mean, list_errors_ards))
list_sem_error_ards = list(map(sem, list_errors_ards))
list_me_error_ards = list(map(lambda x: 1.96 * x, list_sem_error_ards))
list_lower_ci_ards = list(map(lambda x: x[0] - x[1], zip(list_mean_error_ards, list_me_error_ards)))
list_upper_ci_ards = list(map(lambda x: x[0] + x[1], zip(list_mean_error_ards, list_me_error_ards)))
list_sample_size_ards = list(map(len, list_errors_ards))
list_pval_ards = [ttest_1samp(errors, 0)[1] for errors in list_errors_ards]
# _, list_pval_corrected_ards = fdrcorrection(list_pval_ards)

dict_cohort_sample_size_ards = dict(zip(list_cohorts_ards, list_sample_size_ards))

gs1 = gsfig[0:45, 0:30]
ax1 = fig.add_subplot(gs1)
ax1.text(-1.05, 1.05, "A", transform=ax1.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')

ax1.axvline(x=0, ymin=0, ymax=0.95, color='black', linestyle='solid', linewidth=plot_linewidth-0.4, zorder=-999)
ax1.set_yticks(range(len(list_cohorts_ards)))
yticklabels = list(map(lambda x: f"{x[0]} ({x[1]})", zip(list_cohorts_ards, list_sample_size_ards)))
ax1.set_yticklabels(yticklabels)
ax1.invert_yaxis()
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax1.text(-1.02, 0.5, "GSE273149", va='center', ha='center', color='black', fontsize=plt.rcParams["font.size"]+1, weight="bold", rotation=90, transform=ax1.transAxes)
ax1.vlines(x=-0.97, ymin=0.02, ymax=0.98, color='black', linewidth=plot_linewidth-0.5, clip_on = False, transform=ax1.transAxes)

for i, cohort in enumerate(list_cohorts):
    if (list_pval_ards[i] < 0.05):
        markercolor = 'tab:red'
        textcolor = 'tab:red'
        weight = 'bold'
    else:
        markercolor = 'grey'
        textcolor = 'k'
        weight = 'normal'
    ax1.plot(list_mean_error_ards[i], i, 'd', color=markercolor, markersize=plot_linewidth+4, zorder=999, markeredgewidth=plot_linewidth-0.7, markeredgecolor=markercolor)
    ax1.hlines(i, list_lower_ci_ards[i], list_upper_ci_ards[i], color=markercolor, linewidth=plot_linewidth)
    ax1.text(1.3, i, f'{list_mean_error_ards[i]:.2f} [{list_lower_ci_ards[i]:.2f}, {list_upper_ci_ards[i]:.2f}]', va='center', ha='center', color=textcolor, weight=weight, transform=ax1.get_yaxis_transform())
    ax1.text(1.7, i, f'{list_pval_ards[i]:.3f}', va='center', ha='center', color=textcolor, weight=weight, transform=ax1.get_yaxis_transform())

ax1.text(-0.4, 1.02, "Cohorts (N)", ha='center', va='center', transform=ax1.transAxes, fontsize=plt.rcParams["font.size"]+1, weight="bold")
ax1.text(1.3, 1.02, "Mean TAA [95% CI]", ha='center', va='center', transform=ax1.transAxes, fontsize=plt.rcParams["font.size"]+1, weight="bold")
ax1.text(1.7, 1.02, "P", ha='center', va='center', transform=ax1.transAxes, fontsize=plt.rcParams["font.size"]+1, weight="bold")

ax1.set_xlabel('TAA (years)', fontsize=plt.rcParams["font.size"]+1)
ax1.set_xticks(range(-20, 121, 20))
ax1.set_xlim(-20, 120)
ax1.tick_params(axis='y', which='both', length=0)

gs2 = gsfig[0:45, 65:100]
ax2 = fig.add_subplot(gs2)
ax2.text(-0.25, 1.05, "B", transform=ax2.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')
ax2.tick_params(axis='y', which='both', length=0)

df_GSE273149 = pd.read_csv(f"{WORKDIR}/Data/GSE273149_prediction_result_KoreanBloodClock.txt", sep="\t")

df_error_GSE273149= pd.read_csv(f"{WORKDIR}/Data/GSE273149_prediction_result_KoreanBloodClock.txt", sep="\t", index_col=[0])
df_meta_GSE273149 = pd.read_csv(f"{WORKDIR}/Metadata/GSE273149.txt", sep="\t", index_col=[0])
df_GSE273149= pd.concat([df_error_GSE273149, df_meta_GSE273149], axis=1).reset_index(drop=False)
df_GSE273149['time'] = pd.Categorical(
    df_GSE273149['time'],
    categories=["Day 1", "Day 3", "Day 7", "Day 10"],
    ordered=True
)

# Count samples with "High" and "Low" CRP levels across stages
df_counts = df_GSE273149.groupby(["time", "disease status"]).size().unstack(fill_value=0)

# Statistical comparison between "High" and "Low" groups at each stage
stage_pvals = {}
for stage in df_counts.index:
    survivor = df_GSE273149[(df_GSE273149["time"] == stage) & (df_GSE273149["disease status"] == "COVID_survivor")]["AgeAccel"]
    non_survivor = df_GSE273149[(df_GSE273149["time"] == stage) & (df_GSE273149["disease status"] == "COVID_nonsurvivor")]["AgeAccel"]
    
    if len(survivor) > 1 and len(non_survivor) > 1:  # Ensure sufficient samples
        pval = ttest_ind(survivor, non_survivor, equal_var=False).pvalue
    else:
        pval = None  # Not enough samples for comparison
    stage_pvals[stage] = pval

df_GSE273149 = df_GSE273149[df_GSE273149["disease status"].str.contains("COVID")]
df_GSE273149["Subject-ID"] = df_GSE273149["Sample_Age"]
dict_ards_group = dict(zip(df_GSE273149["Subject-ID"], df_GSE273149["disease status"]))

list_subject = df_GSE273149["Subject-ID"].to_list()
dict_subj_ards_color = dict()
for subj in list_subject:
    ards_group = dict_ards_group.get(subj, "Unknown")
    if ards_group == "COVID_survivor":
        color="tab:green"
    elif ards_group == "COVID_nonsurvivor":
        color="tab:orange"
    else:
        color="tab:blue"
    
    dict_subj_ards_color[subj] = color

sns.pointplot(
    data=df_GSE273149,
    x="time",
    y="AgeAccel",
    marker="o",
    hue="disease status",
    palette={"COVID_survivor":"tab:green", "COVID_nonsurvivor":"tab:orange"},
    markersize=5,
    linewidth=2,
    capsize=0.05,
    errorbar="se",
    zorder=3,
    ax=ax2
)

list_subject = df_GSE273149["Subject-ID"].to_list()
for subj in list_subject:
    df_GSE273149_copy = df_GSE273149.copy()
    subj_data = df_GSE273149_copy[df_GSE273149_copy["Subject-ID"] == subj]
    trend_color = dict_subj_ards_color.get(subj, "tab:blue")
    sns.lineplot(
        data=subj_data,
        x="time",
        y="AgeAccel",
        marker="o",
        markersize=1,
        markerfacecolor=trend_color,
        markeredgecolor=None,
        color=trend_color,
        linewidth=1,
        alpha=0.1,
        zorder=-1,
        ax=ax2
    )

ax2.set_ylim(0, 130)
# xticklabels = list_stages_shortened
# ax4.set_xticklabels(xticklabels)
ax2.set_xlabel("Stages", fontsize=plt.rcParams["font.size"]+2, labelpad=20)
ax2.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+2)
ax2.margins(0.1)
ax2.legend(frameon=False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
list_pval = list()
for idx, stage in enumerate(df_counts.index):
    survivor_cnt = df_counts.loc[stage, "COVID_survivor"]
    nonsurvivor_cnt = df_counts.loc[stage, "COVID_nonsurvivor"]
    ax2.text(idx, -0.12, f"Surv. ({survivor_cnt})", ha="center", va="bottom", color="black", transform=ax2.get_xaxis_transform())
    ax2.text(idx, -0.16, f"Non-Surv. ({nonsurvivor_cnt})", ha="center", va="bottom", color="black", transform=ax2.get_xaxis_transform())
    
    # Add p-value annotation
    pval = stage_pvals[stage]
    if pval is not None:
        list_pval.append(pval)

for idx, (stage, pval) in enumerate(zip(df_counts.index, list_pval)):
    if pval is not None and pval <= 0.05:
        ax2.text(idx, 0, f"{pval:.2f}", ha="center", va="bottom", fontsize=10, transform=ax2.get_xaxis_transform(), color="red", weight="bold")
    elif pval is not None and pval > 0.05:
        ax2.text(idx, 0.01, f"{pval:.2f}", ha="center", va="bottom", fontsize=8, transform=ax2.get_xaxis_transform(), color="black", weight="bold")
    else:
        continue

fig.subplots_adjust(left=0.1, hspace=1, wspace=1)
plt.savefig(f"{WORKDIR}/Figures/Figure_3.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{WORKDIR}/Figures/Figure_3.tiff", dpi=300, bbox_inches="tight")
plt.show()
plt.close()
# %%
