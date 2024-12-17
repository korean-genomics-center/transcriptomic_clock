# %%
import glob
import json
import os
from collections import Counter
import statsmodels.api as sm

import compare_clinical_values as crf
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from scipy.stats import mannwhitneyu, sem, ttest_1samp
from statsmodels.stats.multitest import fdrcorrection
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

list_test_group = list(dict_name_conversion.keys())
for i, testname in enumerate(list_test_group):
    file_testing_data = f"{dir_trained_model}/{testname}_std_scaled.txt"
    os.makedirs(dir_trained_model, exist_ok=True)
    clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
    df_clock = clock.make_dataframe_error_statistics_per_sample()
    df_clock.to_csv(os.path.join(dir_trained_model, f"sample_statistics_{testname}_std_scaled.txt"), sep="\t", index=False)
    print(os.path.join(dir_trained_model, f"sample_statistics_{testname}_std_scaled.txt"))

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
_, list_pval_corrected = fdrcorrection(list_pval)
dict_cohort_sample_size = dict(zip(list_cohorts, list_sample_size))

# %%
path_sev_info = f"{WORKDIR}/Clinical/infectomics_CRF_20230410_edit.xlsx"
path_acc_v1 = f"{dir_trained_model}/sample_statistics_viral_v1_prev_std_scaled.txt"
# path_acc_v2 = f"{dir_trained_model}/sample_statistics_viral_v2_prev_std_scaled.txt"
# path_acc_v3 = f"{dir_trained_model}/sample_statistics_viral_v3_prev_std_scaled.txt"
df_excel = pd.read_excel(path_sev_info, engine="openpyxl", skiprows=1)
df_crf = df_excel.rename(columns={"Subject NO.(고유번호)": "Project-ID"})
df_crf["Individual"] = df_crf["Project-ID"].apply(lambda x: "-".join(x.split("-")[:-1]))
df_crf["sample"] = df_crf["Project-ID"]
df_crf["CRP"] = df_crf["CRP"].apply(lambda x: float(str(x).replace("0..31", "0.31")))
df_crf_v1 = df_crf[df_crf["VISIT_no."]==1]
dict_v1_crp = dict(zip(df_crf_v1["Individual"], df_crf_v1["CRP"]))
dict_crp_group = crf.group_patients_by_crp_level(dict_v1_crp)
list_crp_group = crf.get_column_crp_group(df_crf, dict_crp_group, colsample="Project-ID")
df_crf["CRP_Group"] = list_crp_group
df_viral_v1 = pd.read_csv(path_acc_v1, sep="\t")
# df_viral_v2 = pd.read_csv(path_acc_v2, sep="\t")
# df_viral_v3 = pd.read_csv(path_acc_v3, sep="\t")
df_merge_v1 = pd.merge(df_crf, df_viral_v1, on="Project-ID", how="inner")
# df_merge_v2 = pd.merge(df_crf, df_viral_v2, on="Project-ID", how="inner")
# df_merge_v3 = pd.merge(df_crf, df_viral_v3, on="Project-ID", how="inner")
# df_merge_all = pd.concat([df_merge_v1, df_merge_v2, df_merge_v3], axis=0)

list_clinical_crp, list_error_crp = crf.get_clinical_and_error(df_merge_v1, clinical="CRP", error="error")
list_clinical_neut, list_error_neut = crf.get_clinical_and_error(df_merge_v1, clinical="Neutrophil", error="error")
list_clinical_lymp, list_error_lymp = crf.get_clinical_and_error(df_merge_v1, clinical="Lymphocytes", error="error")

list_clinical_neut_hig_crp, list_clinical_neut_low_crp, list_error_neut_high_crp, list_error_neut_low_crp = crf.get_clinical_and_error_divided_by_CRP_group(df_merge_v1, group="CRP_Group", clinical="Neutrophil", error="error")
stat, pval = mannwhitneyu(list_error_neut_high_crp, list_error_neut_low_crp, alternative='greater')
dict_cnt = dict(zip(["High", "Low"], [len(list_clinical_neut_hig_crp), len(list_clinical_neut_low_crp)]))
list_clinical_lymp_hig_crp, list_clinical_lymp_low_crp, list_error_lymp_high_crp, list_error_lymp_low_crp = crf.get_clinical_and_error_divided_by_CRP_group(df_merge_v1, group="CRP_Group", clinical="Lymphocytes", error="error")

r_crp, p_crp = crf.calculate_correlation(list_clinical_crp, list_error_crp)
r_lymp, p_lymp = crf.calculate_correlation(list_clinical_lymp, list_error_lymp)
r_neut, p_neut = crf.calculate_correlation(list_clinical_neut, list_error_neut)

r_neut_high_crp, p_neut_high_crp = crf.calculate_correlation(list_clinical_neut_hig_crp, list_error_neut_high_crp)
r_neut_low_crp, p_neut_low_crp = crf.calculate_correlation(list_clinical_neut_low_crp, list_error_neut_low_crp)
r_lymp_high_crp, p_lymp_high_crp = crf.calculate_correlation(list_clinical_lymp_hig_crp, list_error_lymp_high_crp)
r_lymp_low_crp, p_lymp_low_crp = crf.calculate_correlation(list_clinical_lymp_low_crp, list_error_lymp_low_crp)

# model_neut = crf.compare_regression_slopes_by_CRP_group(df_merge_v1, var_x = "Neutrophil", var_y = "error", group="CRP_Group")
# effect_size_interaction_neut = model_neut.params['interaction']
# pval_interaction_neut = model_neut.pvalues['interaction']

# model_lymp = crf.compare_regression_slopes_by_CRP_group(df_merge_v1, var_x = "Lymphocytes", var_y = "error", group="CRP_Group")
# effect_size_interaction_lymp = model_lymp.params['interaction']
# pval_interaction_lymp = model_lymp.pvalues['interaction']

# %%
cm = 1/2.54
path_font = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path_font)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 7
width = 13*cm
height = 13*cm
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

gs1 = gsfig[0:50, 0:20]
ax1 = fig.add_subplot(gs1)
ax1.text(-0.3, 1.1, "A", transform=ax1.transAxes,
        fontsize=plt.rcParams["font.size"]+1, fontweight='bold', va='top', ha='right')

gs2 = gsfig[0:50, 30:60]
ax2 = fig.add_subplot(gs2)
ax2.text(-0.3, 1.1, "B", transform=ax2.transAxes,
        fontsize=plt.rcParams["font.size"]+1, fontweight='bold', va='top', ha='right')

gs3 = gsfig[0:50, 70:100] 
ax3 = fig.add_subplot(gs3)
ax3.text(-0.3, 1.1, "C", transform=ax3.transAxes,
        fontsize=plt.rcParams["font.size"]+1, fontweight='bold', va='top', ha='right')


sns.barplot(df_merge_v1, 
               x="CRP_Group", 
               y = "error", 
               hue="CRP_Group",
               err_kws = {"linewidth": 1},
               ax=ax1)
ax1.text(0.5, 102, s=f"P={round(float(pval), 4)}", ha="center", weight="bold", fontsize=plt.rcParams["font.size"]+1)
ax1.set_ylim(0, 100)
ax1.hlines(100, xmin=0, xmax=1, color="k", linewidth=0.5)
ax1.set_xlabel("CRP Levels", fontsize=plt.rcParams["font.size"]+1)
ax1.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+1)
ax1.set_xticklabels([f"{group}\n(N={count})"for group, count in dict_cnt.items()])
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

legend_elements = [
    plt.Line2D([0], [0], 
               linestyle='-',
               color='tab:orange' if 'Low' in group else 'tab:blue',
               lw=1,
               label=f'{group} (N={count})',
               markerfacecolor="None", 
               markersize=plt.rcParams["font.size"])
    for group, count in dict_cnt.items()]

crf.draw_regplot(df_merge_v1,  
           x="Lymphocytes", 
           y="error", 
           hue="CRP_Group", 
           scatter_kws={"s": plot_linewidth+11,
                        "marker":'o',
                        "facecolor":"None",
                        "linewidths":0.5},
           line_kws={"linewidth": 1,},
           ci=95,
           ax=ax2)

ax2.set_xlabel("Lymphocyte Count (%)", fontsize=plt.rcParams["font.size"]+1)
ax2.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+1)
ax2.set_xlim(0, 80)
ax2.set_ylim(-50, 250)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.text(35, 200, f"R={round(r_lymp_high_crp,1)}, P={round(p_lymp_high_crp, 3)}", color="tab:blue", weight="bold", fontsize=plt.rcParams["font.size"]+1)
ax2.text(35, 180, f"R={round(r_lymp_low_crp,1)}, P={round(p_lymp_low_crp, 3)}", color="tab:orange", weight="bold", fontsize=plt.rcParams["font.size"]+1)


crf.draw_regplot(df_merge_v1, 
           x="Neutrophil", 
           y="error",
           hue="CRP_Group", 
           scatter_kws={"s": plot_linewidth+11,
                        "marker":'o',
                        "facecolor":"None",
                        "linewidths":0.6},
           line_kws={"linewidth": 1,},
           ci=95,
           ax=ax3)

ax3.set_xlabel("Neutrophil Count (%)", fontsize=plt.rcParams["font.size"]+1)
ax3.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"]+1)
ax3.set_xlim(20, 100)
ax3.set_ylim(-50, 250)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.legend(handles=legend_elements, 
           loc='upper left', 
           title='CRP Levels', 
           title_fontsize=plt.rcParams["font.size"]+1, 
           fontsize=plt.rcParams["font.size"], 
           bbox_to_anchor=(1., .6),  
           frameon=False)
ax3.text(25, 200, f"R={round(r_neut_high_crp,1)}, P={round(p_neut_high_crp, 3)}", color="tab:blue", weight="bold", fontsize=plt.rcParams["font.size"]+1)
ax3.text(25, 180, f"R={round(r_neut_low_crp,1)}, P={round(p_neut_low_crp, 3)}", color="tab:orange", weight="bold", fontsize=plt.rcParams["font.size"]+1)

# %%
