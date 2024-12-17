# %%
import os
import numpy as np
import pandas as pd
from utility import Utility

# %%
WORKDIR = ".."
DIR_META = os.path.join(WORKDIR, "Metadata")
DIR_FEAT = os.path.join(WORKDIR, "FeatureTable")
colsample = "Project-ID"
col_y = "Sample_Age"

# %%
file_original_data = os.path.join(DIR_FEAT, "healthy_cohort1.txt")
file_preprocess_train_data = os.path.join(DIR_FEAT, "healthy_cohort_train_prev.txt")
file_train_data = os.path.join(DIR_FEAT, "healthy_cohort_train_prev_added_back_10_80_90s.txt")
file_sample_10_80_90s = os.path.join(DIR_META, "list_samples_removed_10_80_90.txt")

# %%
# file_original_data = os.path.join(DIR_FEAT, "healthy_cohort1.txt")
# file_preprocess_train_data = os.path.join(DIR_FEAT, "healthy_cohort1_removed_20s_selected_out_10_80_90s_random_split_80.txt")
# file_train_data = os.path.join(DIR_FEAT, "healthy_cohort1_removed_20s_selected_out_10_80_90s_random_split_80_added_back_10_80_90s.txt")
# file_sample_10_80_90s = os.path.join(DIR_META, "list_samples_removed_10_80_90.txt")

# file_preprocess_test_data = os.path.join(DIR_FEAT, "healthy_cohort1_removed_20s_selected_out_10_80_90s_random_split_20.txt")
# file_test_data = os.path.join(DIR_FEAT, "healthy_cohort1_removed_20s_selected_out_10_80_90s_random_split_20_added_back_20s.txt")
# file_sample_20s = os.path.join(DIR_META, "list_samples_removed_20.txt")

# %%
util = Utility()
df_original_exp, df_original_meta = util.get_feature_table(file_original_data)
df_original_exp_, df_sample_age = util.filter_input_data(df_original_exp, df_original_meta, col_y)
df_original_data_colage = pd.concat([df_original_exp_, df_sample_age], axis=1)
df_original_data_indexed = df_original_data_colage.set_index(colsample)

# %%
with open(file_sample_10_80_90s, mode="r") as fr:
    record = fr.readlines()
    list_sample_10_80_90s = list(map(lambda x: x.rstrip("\n"), record))
df_sample_10_80_90s = df_original_data_indexed.loc[list_sample_10_80_90s].reset_index(drop=False)
df_preprocess_train_data = util.read_expmtx(file_preprocess_train_data)
if np.shape(df_sample_10_80_90s)[1] != np.shape(df_preprocess_train_data)[1]:
    X, ys = util.get_feature_table(file_preprocess_train_data)
    X_, y = util.filter_input_data(X, ys, col_y)
    df_preprocess_train_data = pd.concat([X_, y], axis=1)
df_train_data = pd.concat([df_preprocess_train_data, df_sample_10_80_90s], axis=0)
util.save_expmtx(df_train_data, file_train_data)

# %%
# with open(file_sample_20s, mode="r") as fr:
#     record = fr.readlines()
#     list_sample_20s = list(map(lambda x: x.rstrip("\n"), record))
# df_sample_20s = df_original_data_indexed.loc[list_sample_20s].reset_index(drop=False)
# df_preprocess_test_data = util.read_expmtx(file_preprocess_test_data)
# if np.shape(df_sample_20s)[1] != np.shape(df_preprocess_test_data)[1]:
#     X, ys = util.get_feature_table(file_preprocess_test_data)
#     X_, y = util.filter_input_data(X, ys, col_y)
#     df_preprocess_test_data = pd.concat([X_, y], axis=1)
# df_test_data = pd.concat([df_preprocess_test_data, df_sample_20s], axis=0)
# util.save_expmtx(df_test_data, file_test_data)
