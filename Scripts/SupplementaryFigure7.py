# %%
import glob
import json
import os
import seaborn as sns
import compare_clinical_values as crf
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.font_manager as fm
import statsmodels.api as sm
import numpy as np
import matplotlib.patches as mpatches
from scipy.stats import sem, ttest_1samp
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

# Prepare data excluding the "Convalescence" phase
df_merge_all_filtered = df_merge_all[~df_merge_all["Stage"].str.contains("Convalescence", na=False)]

# Redefine stages to exclude "Convalescence"
list_stages = [stage for stage in dict_stage_info.values() if stage != "Convalescence"]

# Update xticklabels excluding "Convalescence"
xticklabels = [
    f"{stage}\n(N={dict_cohort_sample_size.get(stage, 'NA')})" for stage in list_stages
]

# Redefine subject list and trend coloring
list_sample = df_merge_all_filtered["Subject-ID"].unique()
dict_subj_bool_dec = {}

for subj in list_sample:
    subj_data = df_merge_all_filtered[df_merge_all_filtered["Subject-ID"] == subj].sort_values("Stage")
    if subj_data["error"].isnull().any() or len(subj_data) < 2:
        dict_subj_bool_dec[subj] = "transparent"
    else:
        trend = "Decreasing" if subj_data["error"].iloc[0] > subj_data["error"].iloc[-1] else "Increasing"
        dict_subj_bool_dec[subj] = trend

df_merge_all_filtered["Color_Trend_Desc"] = df_merge_all_filtered["Subject-ID"].map(dict_subj_bool_dec)

# Plotting settings
cm = 1 / 2.54
path_font = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path_font)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 10
figsize = (8 * cm, 16 * cm)

fig, (ax_line) = plt.subplots(
    nrows=1,
    figsize=figsize,
)

# Lineplot for individual subjects
for subj in list_sample:
    subj_data = df_merge_all_filtered[df_merge_all_filtered["Subject-ID"] == subj]
    trend_color = dict_subj_bool_dec[subj]
    if trend_color != "transparent":
        sns.lineplot(
            data=subj_data,
            x="Stage",
            y="error",
            marker="o",
            color="tab:orange" if trend_color == "Decreasing" else "tab:blue",  # Keep the color scheme
            linewidth=1,
            ax=ax_line,
        )

# Add the regression line for the overall trend
# To compute the regression line, we need to create a numeric representation of the stages
df_merge_all_filtered["Stage_numeric"] = pd.Categorical(df_merge_all_filtered["Stage"]).codes

# Fit the regression model
X = sm.add_constant(df_merge_all_filtered["Stage_numeric"])  # Add constant for the intercept
y = df_merge_all_filtered["error"]

model = sm.OLS(y, X).fit()
slope = model.params[1]  # Getting the slope (change in error with respect to stage)
p_value = model.pvalues[1]  # Getting the p-value for the slope

# Custom legend handles to fix the shape
legend_handles = [
    mpatches.Patch(hatch="o", color="tab:orange", label="Decreasing"),
    mpatches.Patch(hatch="o", color="tab:blue", label="Increasing")
]

# Update legend with fixed shapes and labels
ax_line.legend(
    handles=legend_handles,
    title="Trend",
    title_fontsize=plt.rcParams["font.size"]-1,
    fontsize=plt.rcParams["font.size"]-2,
    frameon=False
)

# Regression plot with 95% CI interval and thickened line
sns.regplot(
    data=df_merge_all_filtered,
    x="Stage_numeric",
    y="error",
    scatter=False,
    ax=ax_line,
    line_kws={"color": "black", "linewidth": 3},
    ci=95
)

# Display slope and p-value on the right side of the line plot
ax_line.text(0.95, 0.75, f'$R$={slope:.2f}, $P$={p_value:.2f}', 
             ha='right', va='top', transform=ax_line.transAxes, fontsize=10)

# Configure line plot without xticklabels
ax_line.set_xticks(range(len(list_stages)))
ax_line.set_xlabel("Stage", fontsize=plt.rcParams["font.size"] + 2)
ax_line.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"] + 1)
ax_line.set_ylim(-20, 250)
ax_line.spines["top"].set_visible(False)
ax_line.spines["right"].set_visible(False)
ax_line.tick_params(axis="x", which="both", length=0)
ax_line.margins(0.1)

# # Boxplots for TAA distributions by trend at each stage
# sns.boxplot(
#     data=df_merge_all_filtered,
#     x="Stage",
#     y="error",
#     hue="Color_Trend_Desc",
#     hue_order=["Decreasing", "Increasing"],  # Updated hue order
#     ax=ax_box,
#     palette={"Decreasing": "green", "Increasing": "red"},  # Color palette for Decreasing and Increasing
#     showmeans=True,
#     meanline=True,
#     meanprops={"linestyle": "--", "color": "black"},
# )

# ax_box.legend().set_visible(False)
# ax_box.set_xticks(range(len(list_stages)))
# ax_box.set_xticklabels(xticklabels)
# ax_box.set_ylim(-20, 300)
# ax_box.set_xlabel("Stages", fontsize=plt.rcParams["font.size"] + 1)
# ax_box.set_ylabel("TAA (years)", fontsize=plt.rcParams["font.size"] + 1)

# # Perform and display the ranksum test for each stage
# for i, stage in enumerate(list_stages):
#     # Data for each trend (Decreasing vs Increasing)
#     data_decreasing = df_merge_all_filtered[(df_merge_all_filtered["Stage"] == stage) & (df_merge_all_filtered["Color_Trend_Desc"] == "Decreasing")]["error"]
#     data_increasing = df_merge_all_filtered[(df_merge_all_filtered["Stage"] == stage) & (df_merge_all_filtered["Color_Trend_Desc"] == "Increasing")]["error"]
    
#     # Perform the rank-sum test
#     if len(data_decreasing) > 0 and len(data_increasing) > 0:
#         stat, p_val = ranksums(data_decreasing, data_increasing)
#         # Display p-value on the plot for each stage
#         if float(p_val) < 0.05:
#             ax_box.text(i, 250, f'{p_val:.3f}', ha='center', fontsize=10, color='red')
#         else:
#             ax_box.text(i, 250, f'{p_val:.3f}', ha='center', fontsize=10, color='black')

# ax_box.spines["top"].set_visible(False)
# ax_box.spines["right"].set_visible(False)
# ax_box.tick_params(axis="both", which="both", length=0)

plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_7.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_7.tiff", dpi=300, bbox_inches="tight")
plt.tight_layout()
plt.show()
 
# %%
# import json
# import matplotlib.font_manager as fm
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.ticker import MaxNLocator
# import glob
# import os
# import numpy as np
# import pandas as pd
# from scipy.stats import sem, ttest_1samp
# from statsmodels.stats.multitest import fdrcorrection
# from transcriptomic_clock import Clock

# # %%
# def get_dictionary_groupname_residual(dir_trained_model, dict_name_conversion, tag_file_sample="sample_statistics"):

#     list_files = glob.glob(os.path.join(dir_trained_model, f"{tag_file_sample}*"), recursive=True)
#     list_files_filtered = list(filter(lambda x: os.path.basename(x).replace(f"{tag_file_sample}_", "").replace("_std_scaled.txt", "") in dict_name_conversion.keys(), list_files)) 
#     dict_groupname_residual = dict()
#     for file in list_files_filtered:
#         groupname = os.path.basename(file).split(tag_file_sample)[-1][1:].replace("_std_scaled.txt", "")
#         groupname_new = dict_name_conversion[groupname]
#         df_residual = pd.read_csv(file, sep="\t")
#         df_residual_copy = df_residual.copy()
#         array_error = df_residual_copy["error"].to_numpy()
#         dict_groupname_residual[groupname_new] = array_error

#     dict_groupname_residual = dict(sorted(dict_groupname_residual.items(), key=lambda x: list(dict_name_conversion.values()).index(x[0])))

#     return dict_groupname_residual

# # %%
# WORKDIR = ".."
# DIR_LASSO_INFO = f"{WORKDIR}/LASSO_INFO"
# train_group = "healthy_cohort_train_prev_median_filter"
# dir_trained_model = f"{DIR_LASSO_INFO}/{train_group}/corr_0.35_above_pval_0.05_below"
# list_test_group = [train_group, "healthy_cohort_valid_prev", "healthy_cohort2"]

# path_dict_name_conversion_json = f"{WORKDIR}/Data/dictionary_name_conversion.json"
# with open (path_dict_name_conversion_json, mode="rb") as fb1:
#     dict_name_conversion = json.load(fb1)

# path_dict_cohort_information_json = f"{WORKDIR}/Data/dictionary_cohort_information.json"
# with open(path_dict_cohort_information_json, mode="rb") as fb2:
#     dict_groups = json.load(fb2)

# list_test_group = list(dict_name_conversion.keys())
# for i, testname in enumerate(list_test_group):
#     file_testing_data = f"{dir_trained_model}/{testname}_std_scaled.txt"
#     os.makedirs(dir_trained_model, exist_ok=True)
#     clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
#     df_clock = clock.make_dataframe_error_statistics_per_sample()
#     df_clock.to_csv(os.path.join(dir_trained_model, f"sample_statistics_{testname}_std_scaled.txt"), sep="\t", index=False)
#     print(os.path.join(dir_trained_model, f"sample_statistics_{testname}_std_scaled.txt"))

# dict_groupname_residual = get_dictionary_groupname_residual(dir_trained_model, dict_name_conversion)
# list_cohorts = list(dict_groupname_residual.keys())
# list_errors = list(dict_groupname_residual.values())
# list_mean_error = list(map(np.mean, list_errors))
# list_sem_error = list(map(sem, list_errors))
# list_me_error = list(map(lambda x: 1.96 * x, list_sem_error))
# list_lower_ci = list(map(lambda x: x[0] - x[1], zip(list_mean_error, list_me_error)))
# list_upper_ci = list(map(lambda x: x[0] + x[1], zip(list_mean_error, list_me_error)))
# list_sample_size = list(map(len, list_errors))
# list_pval = [ttest_1samp(errors, 0)[1] for errors in list_errors]
# _, list_pval_corrected = fdrcorrection(list_pval)
# dict_cohort_sample_size = dict(zip(list_cohorts, list_sample_size))


#  # %% 
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.ticker import MaxNLocator
# from matplotlib import font_manager as fm

# cm = 1 / 2.54
# path = f"{WORKDIR}/Arial.ttf"
# prop = fm.FontProperties(fname=path)
# plt.rcParams['font.family'] = prop.get_name()
# plt.rcParams["font.size"] = 5

# width = 15 * cm
# height = 12 * cm

# # Determine the number of cohorts and subplots
# num_cohorts = len(list_cohorts)
# num_rows = int(np.ceil(np.sqrt(num_cohorts))) - 1
# num_cols = int(np.ceil(np.sqrt(num_cohorts)))

# # Create a figure and subplots
# fig, axes = plt.subplots(num_rows, num_cols, figsize=(width, height), constrained_layout=True)

# # Ensure axes is iterable
# if num_rows * num_cols == 1:
#     axes = [axes]  # Wrap single Axes in a list
# else:
#     axes = axes.flatten()  # Flatten array of Axes for easy iteration

# # Iterate through cohorts and plot each histogram
# for i, (ax, cohort) in enumerate(zip(axes, list_cohorts)):
#     error = dict_groupname_residual[cohort]
#     mean_error = round(np.mean(error), 2)
    
#     # Plot histogram
#     ax.hist(error, color="grey", alpha=0.5)
#     ax.axvline(x=mean_error, color="tab:blue", linestyle="dotted")
#     ax.text(mean_error + 1, 1, f"Mean Error\n={mean_error}", 
#             color="tab:blue", weight="bold", fontsize=plt.rcParams["font.size"] + 1)
#     ax.set_xlabel("Prediction Error (years)", fontsize=plt.rcParams["font.size"])
#     ax.set_ylabel("Proportion (%)", fontsize=plt.rcParams["font.size"])
#     ax.set_title(cohort, fontsize=plt.rcParams["font.size"] + 1, weight="bold", pad=0.1)
#     ax.set_xlim(int(min(error)) - 2, int(max(error)) + 2)
#     ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    
#     # Label subplot with letters a-z
#     subplot_letter = chr(ord('A') + i)  # Convert integer index to character
#     ax.text(-0.25, 1.1, subplot_letter, transform=ax.transAxes,
#             fontsize=plt.rcParams["font.size"] + 1, fontweight='bold', va='top', ha='right')

# # Hide unused subplots if any
# for j in range(len(list_cohorts), len(axes)):
#     axes[j].axis('off')

# # Adjust layout and show plot
# plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_7.png", dpi=300, bbox_inches="tight")
# plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_7.tiff", dpi=300, bbox_inches="tight")
# plt.show()
# plt.close()
# # %%

# %%
