# %%
import glob
import json
import os
import seaborn as sns
import compare_clinical_values as crf
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from scipy.stats import sem, ttest_1samp, ttest_ind
from transcriptomic_clock import Clock


# %%
def get_dictionary_groupname_residual(dir_trained_model, dict_name_conversion, tag_file_sample="sample_statistics"):

    list_files = glob.glob(os.path.join(dir_trained_model, f"{tag_file_sample}*"), recursive=True)
    list_files_filtered = list(filter(lambda x: os.path.basename(x).replace(f"{tag_file_sample}_", "").replace("_std_scaled.txt", "") in dict_name_conversion.keys(), list_files)) 
    dict_groupname_residual = dict()
    for file in list_files_filtered:
        groupname = os.path.basename(file).split(tag_file_sample)[-1][1:].replace("_std_scaled.txt", "")
        groupname_new = dict_name_conversion[groupname]
        df_residual = pd.read_csv(file, sep="\t")
        df_residual_copy = df_residual.copy()
        array_error = df_residual_copy["error"].to_numpy()
        dict_groupname_residual[groupname_new] = array_error

    dict_groupname_residual = dict(sorted(dict_groupname_residual.items(), key=lambda x: list(dict_name_conversion.values()).index(x[0])))

    return dict_groupname_residual

def get_error_statistics_per_sample(dir_trained_model, dict_name_conversion):
    list_test_group = list(dict_name_conversion.keys())
    for i, testname in enumerate(list_test_group):
        file_testing_data = f"{dir_trained_model}/{testname}_std_scaled.txt"
        os.makedirs(dir_trained_model, exist_ok=True)
        clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
        df_clock = clock.make_dataframe_error_statistics_per_sample()
        df_clock.to_csv(os.path.join(dir_trained_model, f"sample_statistics_{testname}_std_scaled.txt"), sep="\t", index=False)
        print(os.path.join(dir_trained_model, f"sample_statistics_{testname}_std_scaled.txt"))

# %%
WORKDIR = ".."
train_group = "healthy_cohort_train_prev_median_filter"
list_test_group = [train_group, "healthy_cohort_valid_prev", "healthy_cohort2"]
dir_trained_model = f"{WORKDIR}/LASSO_INFO/{train_group}/corr_0.35_above_pval_0.05_below"

path_dict_name_conversion_json = f"{WORKDIR}/Data/dictionary_name_conversion.json"
with open (path_dict_name_conversion_json, mode="rb") as fb1:
    dict_name_conversion = json.load(fb1)

path_dict_cohort_information_json = f"{WORKDIR}/Data/dictionary_cohort_information.json"
with open(path_dict_cohort_information_json, mode="rb") as fb2:
    dict_groups = json.load(fb2)

get_error_statistics_per_sample(dir_trained_model, dict_name_conversion)

path_sev_info = f"{WORKDIR}/Clinical/infectomics_CRF_20230410_edit.xlsx"
path_acc_v1 = f"{dir_trained_model}/sample_statistics_viral_v1_long_std_scaled.txt"
path_acc_v2 = f"{dir_trained_model}/sample_statistics_viral_v2_long_std_scaled.txt"
path_acc_v3 = f"{dir_trained_model}/sample_statistics_viral_v3_long_std_scaled.txt"
path_acc_recov = f"{dir_trained_model}/sample_statistics_viral_recov_prev_std_scaled.txt"

df_excel = pd.read_excel(path_sev_info, engine="openpyxl", skiprows=1)
df_crf = df_excel.rename(columns={"Subject NO.(고유번호)": "Project-ID"})
df_crf["Subject-ID"] = df_crf["Project-ID"].apply(lambda x: "-".join(x.split("-")[:-1]))
df_crf["CRP"] = df_crf["CRP"].apply(lambda x: float(str(x).replace("0..31", "0.31")))
df_crf_v1 = df_crf[df_crf["VISIT_no."] == 1]
dict_v1_crp = dict(zip(df_crf_v1["Subject-ID"], df_crf_v1["CRP"]))
dict_crp_group = crf.group_patients_by_crp_level(dict_v1_crp)
list_crp_group = crf.get_column_crp_group(df_crf, dict_crp_group, colsample="Subject-ID")
from collections import Counter
dict_crp_cnt = dict(Counter(list_crp_group))
df_crf["CRP_Group"] = list_crp_group
cond_filt = np.logical_and(df_crf["Project-ID"].str.split("-").str[1].str.startswith("R"), df_crf["VISIT_no."].astype(str)=="1")
df_crf.loc[cond_filt, "VISIT_no."] = 5
df_viral_v1 = pd.read_csv(path_acc_v1, sep="\t")
df_merge_v1 = pd.merge(df_crf, df_viral_v1, on="Project-ID", how="inner")
df_viral_v2 = pd.read_csv(path_acc_v2, sep="\t")
df_merge_v2 = pd.merge(df_crf, df_viral_v2, on="Project-ID", how="inner")
df_viral_v3 = pd.read_csv(path_acc_v3, sep="\t")
df_merge_v3 = pd.merge(df_crf, df_viral_v3, on="Project-ID", how="inner")
df_viral_recov = pd.read_csv(path_acc_recov, sep="\t")
df_merge_recov = pd.merge(df_crf, df_viral_recov, on="Project-ID", how="inner")
df_merge_all = pd.concat([df_merge_v1, df_merge_v2, df_merge_v3, df_merge_recov], axis=0)
dict_stage_info = {1: "Acute Phase", 2: "Mid Phase", 3: "Late Phase", 5: "Convalescence"}
df_merge_all["Subject-ID"] = df_merge_all["Project-ID"].apply(lambda x: "-".join(x.split("-")[:-1]))
df_merge_all["Stage"] = df_merge_all["VISIT_no."].apply(lambda x: dict_stage_info.get(x, x))
df_merge_all['Stage'] = pd.Categorical(
    df_merge_all['Stage'], 
    categories=['Acute Phase', 'Mid Phase', 'Late Phase', 'Convalescence'],
    ordered=True
)
df_merge_all = df_merge_all.sort_values(['Subject-ID', 'Stage'])

# %%
dict_groupname_residual = get_dictionary_groupname_residual(dir_trained_model, dict_name_conversion)
list_cohorts = list(dict_groupname_residual.keys())
list_errors = list(dict_groupname_residual.values())
list_mean_error = list(map(np.mean, list_errors))
list_sem_error = list(map(sem, list_errors))
list_me_error = list(map(lambda x: 1.96 * x, list_sem_error))
list_lower_ci = list(map(lambda x: x[0] - x[1], zip(list_mean_error, list_me_error)))
list_upper_ci = list(map(lambda x: x[0] + x[1], zip(list_mean_error, list_me_error)))
list_sample_size = list(map(len, list_errors))
list_pval = [ttest_1samp(errors, 0)[1] for errors in list_errors]
dict_cohort_sample_size = dict(zip(list_cohorts, list_sample_size))

 # %%
cm = 1/2.54
path_font = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path_font)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 7
width = 15*cm
height = 12*cm
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

gs1 = gsfig[0:100, 0:30]
ax1 = fig.add_subplot(gs1)
ax1.text(-0.9, 1.05, "A", transform=ax1.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')

ax1.axvline(x=0, ymin=0, ymax=0.95, color='black', linestyle='solid', linewidth=plot_linewidth-0.4, zorder=-999)
ax1.set_yticks(range(len(list_cohorts)))
yticklabels = list(map(lambda x: f"{x[0]} ({x[1]})", zip(list_cohorts, list_sample_size)))
ax1.set_yticklabels(yticklabels)
ax1.invert_yaxis()
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['top'].set_visible(False)

for group, cohorts in dict_groups.items():
    first_cohort_index = list_cohorts.index(cohorts[0])
    group_label_y = (first_cohort_index + list_cohorts.index(cohorts[-1])) / 2
    ax1.text(-64, group_label_y, group, va='center', ha='center', color='black', fontsize=plt.rcParams["font.size"]+1, weight="bold", rotation=90)
   
line_adjust = 0.2
ax1.vlines(x=-60, ymin=0-line_adjust, ymax=2+line_adjust, color='black', linewidth=plot_linewidth-0.5, clip_on = False)
ax1.vlines(x=-60, ymin=3-line_adjust, ymax=6+line_adjust, color='black', linewidth=plot_linewidth-0.5, clip_on = False)
ax1.vlines(x=-60, ymin=7-line_adjust, ymax=9+line_adjust, color='black', linewidth=plot_linewidth-0.5, clip_on = False)

for i, cohort in enumerate(list_cohorts):
    if (list_pval[i] < 0.05):
        markercolor = 'tab:red'
        textcolor = 'tab:red'
        weight = 'bold'
    else:
        markercolor = 'grey'
        textcolor = 'k'
        weight = 'normal'

    ax1.plot(list_mean_error[i], i, 'd', color=markercolor, markersize=plot_linewidth+4, zorder=999, markeredgewidth=plot_linewidth-0.7, markeredgecolor=markercolor)
    ax1.hlines(i, list_lower_ci[i], list_upper_ci[i], color=markercolor, linewidth=plot_linewidth)
    ax1.text(1.2, i, f'{list_mean_error[i]:.2f} [{list_lower_ci[i]:.2f}, {list_upper_ci[i]:.2f}]', va='center', ha='center', color=textcolor, weight=weight, transform=ax1.get_yaxis_transform())
    ax1.text(1.7, i, f'{list_pval[i]:.3f}', va='center', ha='center', color=textcolor, weight=weight, transform=ax1.get_yaxis_transform())

ax1.text(-0.4, 1.0, "Cohorts (N)", ha='center', va='center', transform=ax1.transAxes, fontsize=plt.rcParams["font.size"], weight="bold")
ax1.text(1.2, 1.0, "Mean TAA [95% CI]", ha='center', va='center', transform=ax1.transAxes, fontsize=plt.rcParams["font.size"], weight="bold")
ax1.text(1.7, 1.0, "P", ha='center', va='center', transform=ax1.transAxes, fontsize=plt.rcParams["font.size"], weight="bold")

ax1.set_xlabel('TAA (years)', fontsize=plt.rcParams["font.size"]+1)
ax1.set_xlim(-11, 46)
ax1.tick_params(axis='y', which='both', length=0)

gs2 = gsfig[0:100, 67:100]
ax2 = fig.add_subplot(gs2)
ax2.text(-0.3, 1.05, "B", transform=ax2.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')
ax2.tick_params(axis='y', which='both', length=0)

df_merge_all_excl_unk = df_merge_all.copy()
df_merge_all_excl_unk = df_merge_all_excl_unk[df_merge_all_excl_unk["CRP_Group"] != "Unknown"]

# Count samples with "High" and "Low" CRP levels across stages
df_counts = df_merge_all_excl_unk.groupby(["Stage", "CRP_Group"]).size().unstack(fill_value=0)

# Statistical comparison between "High" and "Low" groups at each stage
stage_pvals = {}
for stage in df_counts.index:
    high_values = df_merge_all_excl_unk[(df_merge_all_excl_unk["Stage"] == stage) & (df_merge_all_excl_unk["CRP_Group"] == "High")]["error"]
    low_values = df_merge_all_excl_unk[(df_merge_all_excl_unk["Stage"] == stage) & (df_merge_all_excl_unk["CRP_Group"] == "Low")]["error"]
    
    if len(high_values) > 1 and len(low_values) > 1:  # Ensure sufficient samples
        pval = ttest_ind(high_values, low_values, equal_var=False).pvalue
    else:
        pval = None  # Not enough samples for comparison
    stage_pvals[stage] = pval

sns.pointplot(
    data=df_merge_all_excl_unk,
    x="Stage",
    y="error",
    hue="CRP_Group",
    hue_order=["High", "Low"],
    marker="o",
    palette={"High":"tab:orange", "Low":"tab:green"},
    markersize=5,
    linewidth=2,
    capsize=0.05,
    errorbar="se",
    zorder=3,
    ax=ax2
)

list_pval_line = list()
for idx, stage in enumerate(df_counts.index):
    high_count = df_counts.loc[stage, "High"]
    low_count = df_counts.loc[stage, "Low"]
    ax2.text(idx, -0.065, f"High ({high_count})", ha="center", va="bottom", color="black", transform=ax2.get_xaxis_transform())
    ax2.text(idx, -0.09, f"Low ({low_count})", ha="center", va="bottom", color="black", transform=ax2.get_xaxis_transform())
    
    # Add p-value annotation
    pval = stage_pvals[stage]
    if pval is not None:
        list_pval_line.append(pval)

for idx, (stage, pval) in enumerate(zip(df_counts.index, list_pval_line)):
    if pval is not None and pval <= 0.05:
        ax2.text(idx, 0.01, f"{pval:.2f}", ha="center", va="bottom", fontsize=8, transform=ax2.get_xaxis_transform(), color="red", weight="bold")
    elif pval is not None and pval > 0.05:
        ax2.text(idx, 0.01, f"{pval:.2f}", ha="center", va="bottom", fontsize=8, transform=ax2.get_xaxis_transform(), color="black", weight="bold")
    else:
        continue

ax2.set_ylim(-0.2, 1.2)

list_subject = df_merge_all["Subject-ID"].to_list()
dict_subj_crp_color = dict()
for subj in list_subject:
    crp_group = dict_crp_group.get(subj, "Unknown")
    if crp_group == "High":
        color="tab:orange"
    elif crp_group == "Low":
        color="tab:green"
    else:
        color="k"
    
    dict_subj_crp_color[subj] = color

for subj in list_subject:
    df_merge_all_copy = df_merge_all.copy()
    subj_data = df_merge_all_copy[df_merge_all_copy["Subject-ID"] == subj]
    trend_color = dict_subj_crp_color.get(subj, "k")
    sns.lineplot(
        data=subj_data,
        x="Stage",
        y="error",
        marker="o",
        markersize=2,
        markerfacecolor=trend_color,
        markeredgecolor=None,
        color=trend_color,
        linewidth=1,
        alpha=0.1,
        zorder=-1,
        ax=ax2
    )

ax2.set_ylim(-30, 250)
list_stages = list(dict_stage_info.values())
ax2.set_xticks(range(len(list_stages)))
list_stages_shortened = list(map(lambda x: x[:6]+"." if len(x) > 12 else x, list_stages))
xticklabels = list_stages_shortened
ax2.set_xticklabels(xticklabels)
ax2.set_xlabel("Stages", fontsize=plt.rcParams["font.size"]+2, labelpad=25)
ax2.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+2)
ax2.margins()
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.margins(0.1)

legend_elements = [
    plt.Line2D(
        [0], [0],
        linestyle='-',
        color='tab:green' if group == 'Low' else 'tab:orange' if group == 'High' else 'k',
        lw=1,
        label=f"{group}",
        marker='o',
        markerfacecolor='tab:green' if group == 'Low' else 'tab:orange' if group == 'High' else 'k',
        markersize=5
    )
    for group in ["Low", "High", "Unknown"]
]

# Add legend to the plot
ax2.legend(
    handles=legend_elements,
    loc='upper left',
    title='CRP Levels',
    title_fontsize=plt.rcParams["font.size"] + 1,
    fontsize=plt.rcParams["font.size"],
    bbox_to_anchor=(0.5, 1),
    frameon=False
)

plt.subplots_adjust(left=0.1)
plt.tight_layout()
plt.savefig(f"{WORKDIR}/Figures/Figure_2.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{WORKDIR}/Figures/Figure_2.tiff", dpi=300, bbox_inches="tight")
plt.show()
plt.close()
# %%
