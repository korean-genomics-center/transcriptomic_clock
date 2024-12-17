# %%
import os
import json
from transcriptomic_clock import Clock

# %% 
WORKDIR = ".."
DIR_DATA = os.path.join(WORKDIR, "Data")
train_group = "healthy_cohort_train_prev_median_filter"
list_test_group = ["healthy_cohort_valid_prev", "healthy_cohort2", "GSE134080"]
dir_trained_model = f"{WORKDIR}/LASSO_INFO/{train_group}/corr_0.35_above_pval_0.05_below"
path_dict_name_conversion_json = f"{WORKDIR}/Data/dictionary_name_conversion.json"
with open (path_dict_name_conversion_json, mode="rb") as fb1:
    dict_name_conversion = json.load(fb1)

# %%
for i in range(len(list_test_group)):
    file_testing_data = f"{dir_trained_model}/{list_test_group[i]}_std_scaled.txt"
    os.makedirs(dir_trained_model, exist_ok=True)
    clock = Clock(dir_trained_model, file_testing_data)
    df_stat = clock.make_dataframe_error_statistics_per_sample()
    df_stat.columns = ["Project-ID", "ChronAge", "RNAAge", "AgeAccel", "AbsAgeAccel"]
    file_stat = os.path.join(DIR_DATA, f"{list_test_group[i]}_prediction_result_KoreanBloodClock.txt")
    df_stat.to_csv(file_stat, sep="\t", index=False)