# %%
import glob
import json
import os
import numpy as np
import pandas as pd
from scipy.stats import sem, ttest_1samp
from statsmodels.stats.multitest import fdrcorrection
from transcriptomic_clock import Clock
from itertools import chain

# %%
WORKDIR = ".."
train_group = "healthy_cohort_train_prev_median_filter"
dir_trained_model = f"{WORKDIR}/LASSO_INFO/{train_group}/corr_0.35_above_pval_0.05_below"
path_dict_name_conversion_json = f"{WORKDIR}/Data/dictionary_name_conversion.json"
with open (path_dict_name_conversion_json, mode="rb") as fb1:
    dict_name_conversion = json.load(fb1)

list_test_group = list(dict_name_conversion.keys())

path_dict_cohort_information_json = f"{WORKDIR}/Data/dictionary_cohort_information.json"
with open (path_dict_cohort_information_json, mode="rb") as fb2:
    dict_groups = json.load(fb2)

for i, testname in enumerate(list_test_group):
    file_testing_data = f"{dir_trained_model}/{testname}_std_scaled.txt"
    os.makedirs(dir_trained_model, exist_ok=True)
    clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
    df_clock = clock.make_dataframe_error_statistics_per_sample()
    df_clock.to_csv(os.path.join(dir_trained_model, f"sample_statistics_{testname}_std_scaled.txt"), sep="\t", index=False)
    print(os.path.join(dir_trained_model, f"sample_statistics_{testname}_std_scaled.txt"))

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

dict_cohort_stats = {"cohort": list_cohorts, "mean_error": list_mean_error, "lower_ci":list_lower_ci, "upper_ci":list_upper_ci, "pval":list_pval}
# %%
# Calculate summarized result for different groups
list_groups = list(dict_groups.keys())

def get_overall_mean_ci_group_cohorts(groupname):
    if groupname == "COVID-19":
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
df_group_stats.to_csv(f"{WORKDIR}/Tables/Table_S3.txt", sep="\t", index=False)