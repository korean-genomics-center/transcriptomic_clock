# %%
import json
import glob
import os
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import compare_clinical_values as crf
from transcriptomic_clock import Clock
import seaborn as sns

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
df_merge_all = df_merge_all.reset_index(drop=True)

#%% 
cm = 1/2.54
path = "/BiO/Share/Fonts/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 7
width = 8*cm
height = 20*cm
figsize = (width, height)
plot_linewidth = 1
list_visit_groups = df_merge_all["VISIT_no."].to_list()
dict_cnt = Counter(list_visit_groups)
phase_map = dict_stage_info
color_map = {1: "tab:blue", 2: "tab:orange", 3: "tab:green", 5: "tab:red"}

legend_elements = [
    plt.Line2D([0], [0], 
               marker='o',
               markeredgecolor=color_map[group],
               lw=0,
               label=f'{phase_map[group]} (N={count})',
               markerfacecolor="None", 
               markersize=plt.rcParams["font.size"])
    for group, count in dict_cnt.items()]

fig, axes = plt.subplots(3, 1, figsize=figsize)
axes = axes.flatten()
list_clinicals = ["CRP", "Neutrophil", "Lymphocytes"]
for i, clinical in enumerate(list_clinicals):
    sns.scatterplot(
        data=df_merge_all,
        x="CRP",
        y="error",
        hue="Stage",
        s=10,  # Marker size
        marker='o',
        facecolor="None",
        linewidth=0.5,
        ax=axes[i]
    )
    axes[i].set_xlabel(clinical, fontsize=plt.rcParams["font.size"] + 1)
    axes[i].set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"] + 1)
    axes[i].text(
        -0.3, 1.1, "a",
        transform=axes[i].transAxes,
        fontsize=plt.rcParams["font.size"] + 2,
        fontweight='bold',
        va='top', 
        ha='right'
    )
    axes[i].legend(
        handles=legend_elements,
        loc='upper left',
        title='Infect. Phase',
        title_fontsize=plt.rcParams["font.size"] + 1,
        fontsize=plt.rcParams["font.size"],
        frameon=False
    )

plt.tight_layout()
# plt.savefig("./Figures/Supplementary_Figure_11.pdf", dpi=330, bbox_inches="tight")
# plt.savefig("./Figures/Supplementary_Figure_11.png", dpi=330, bbox_inches="tight")
plt.show()
plt.close()

# %%
