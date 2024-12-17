# %%
import os
import pandas as pd
import numpy as np

# %%
def write_sample_list_file(file_samples, list_samples):
    with open(file_samples, mode="w") as fw:
        for sample in list_samples:
            fw.write(sample + "\n")

# %%
path_healthy_cohort1 = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/healthy_cohort1.txt"
df_healthy_cohort1 = pd.read_csv(path_healthy_cohort1, sep="\t", index_col=0)
list_samples_healthy_cohort1 = list(df_healthy_cohort1.columns)
path_stress = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/stress.txt"
df_stress = pd.read_csv(path_stress, sep="\t", index_col=0)
list_samples_stress = list(df_stress.columns)
path_viral = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/viral.txt"
df_viral = pd.read_csv(path_viral, sep="\t", index_col=0)
list_samples_viral = list(df_viral.columns)

# %%
path_metadata_healthy = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/healthy.txt"
df_meta_healthy = pd.read_csv(path_metadata_healthy, sep="\t", index_col=0)
path_metadata_stress = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/stress.txt"
df_meta_stress = pd.read_csv(path_metadata_stress, sep="\t", index_col=0)
path_metadata_viral = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/viral.txt"
df_meta_viral = pd.read_csv(path_metadata_viral, sep="\t", index_col=0)

# %%
path_samples_qc_fail = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/list_qc_fail.txt"
with open(path_samples_qc_fail, mode="r") as fr:
    list_samples_qc_fail = list(map(lambda x: x.rstrip("\n"), fr.readlines()))

dir_sample_list = "/BiO/Research/GeroExpressome/Resources/Miscellaneous"
dir_sample_list_intsc = "/BiO/Access/kyungwhan1998/transcriptomic_clock/Data"

path_samples_healthy_train = f"{dir_sample_list}/sample_list_healthy_train_prev.txt"
path_samples_healthy_train_intsc  = f"{dir_sample_list_intsc}/sample_list_healthy_train_prev.txt"
with open(path_samples_healthy_train, mode="r") as fr:
    list_samples_healthy_train = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
list_intsc_train = list(set(list_samples_healthy_cohort1).intersection(set(list_samples_healthy_train)))
write_sample_list_file(path_samples_healthy_train_intsc, list_intsc_train)
list_intsc_train_qc_fail = list(set(list_samples_qc_fail).intersection(set(list_samples_healthy_train)))

path_samples_healthy_validation = f"{dir_sample_list}/sample_list_healthy_valid_prev.txt"
path_samples_healthy_validation_intsc  = f"{dir_sample_list_intsc}/sample_list_healthy_valid_prev.txt"
with open(path_samples_healthy_validation, mode="r") as fr:
    list_samples_healthy_valid = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
list_intsc_valid = list(set(list_samples_healthy_cohort1).intersection(set(list_samples_healthy_valid)))
write_sample_list_file(path_samples_healthy_validation_intsc, list_intsc_valid)
list_intsc_valid_qc_fail = list(set(list_samples_qc_fail).intersection(set(list_samples_healthy_valid)))

path_samples_v1 = f"{dir_sample_list}/sample_list_viral_v1_long_prev.txt"
path_samples_v1_intsc = f"{dir_sample_list_intsc}/sample_list_viral_v1_long_prev.txt"
with open(path_samples_v1, mode="r") as fr:
    list_samples_v1 = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
list_intsc_v1 = list(set(list_samples_v1).intersection(set(list_samples_viral)))
write_sample_list_file(path_samples_v1_intsc, list_intsc_v1)

path_samples_v2 = f"{dir_sample_list}/sample_list_viral_v2_long_prev.txt"
path_samples_v2_intsc = f"{dir_sample_list_intsc}/sample_list_viral_v2_long_prev.txt"
with open(path_samples_v2, mode="r") as fr:
    list_samples_v2 = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
list_intsc_v2 = list(set(list_samples_v2).intersection(set(list_samples_viral)))
write_sample_list_file(path_samples_v2_intsc, list_intsc_v2)

path_samples_v3 = f"{dir_sample_list}/sample_list_viral_v3_long_prev.txt"
path_samples_v3_intsc = f"{dir_sample_list_intsc}/sample_list_viral_v3_long_prev.txt"
with open(path_samples_v3, mode="r") as fr:
    list_samples_v3 = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
list_intsc_v3 = list(set(list_samples_v3).intersection(set(list_samples_viral)))
write_sample_list_file(path_samples_v3_intsc, list_intsc_v3)

path_samples_viral_recov = f"{dir_sample_list}/sample_list_viral_recov_prev.txt"
path_samples_viral_recov_intsc = f"{dir_sample_list_intsc}/sample_list_viral_recov_prev.txt"
with open(path_samples_viral_recov, mode="r") as fr:
    list_samples_recov = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
list_intsc_recov = list(set(list_samples_recov).intersection(set(list_samples_viral)))
write_sample_list_file(path_samples_viral_recov_intsc, list_intsc_recov)

path_samples_anxiety = f"{dir_sample_list}/sample_list_anxiety_prev.txt"
path_samples_anxiety_intsc = f"{dir_sample_list_intsc}/sample_list_anxiety_prev.txt"
with open(path_samples_anxiety, mode="r") as fr:
    list_samples_anxiety = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
list_intsc_anxiety = list(set(list_samples_anxiety).intersection(set(list_samples_stress)))
write_sample_list_file(path_samples_anxiety_intsc, list_intsc_anxiety)
list_intsc_anxiety_qc_fail = list(set(list_samples_qc_fail).intersection(set(list_samples_anxiety)))

path_samples_suicide = f"{dir_sample_list}/sample_list_attempt_prev.txt"
path_samples_suicide_intsc = f"{dir_sample_list_intsc}/sample_list_attempt_prev.txt"
with open(path_samples_suicide, mode="r") as fr:
    list_samples_suicide = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
list_intsc_suicide = list(set(list_samples_suicide).intersection(set(list_samples_stress)))
write_sample_list_file(path_samples_suicide_intsc, list_intsc_suicide)
list_intsc_suicide_qc_fail = list(set(list_samples_qc_fail).intersection(set(list_samples_suicide)))

path_samples_depression = f"{dir_sample_list}/sample_list_depression_prev.txt"
path_samples_depression_intsc = f"{dir_sample_list_intsc}/sample_list_depression_prev.txt"
with open(path_samples_depression, mode="r") as fr:
    list_samples_depression = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
list_intsc_depression = list(set(list_samples_depression).intersection(set(list_samples_stress)))
write_sample_list_file(path_samples_depression_intsc, list_intsc_depression)
list_intsc_depression_qc_fail = list(set(list_samples_qc_fail).intersection(set(list_samples_depression)))


# %%
path_healthy_cohort1_train_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/healthy_cohort_train_prev.txt"
df_healthy_cohort1_train = df_healthy_cohort1[list_intsc_train]
df_healthy_cohort1_train.to_csv(path_healthy_cohort1_train_prev, sep="\t")

path_healthy_cohort1_valid_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/healthy_cohort_valid_prev.txt"
df_healthy_cohort1_valid = df_healthy_cohort1[list_intsc_valid]
df_healthy_cohort1_valid.to_csv(path_healthy_cohort1_valid_prev, sep="\t")

path_viral_v1_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/viral_v1_prev.txt"
df_viral_v1 = df_viral[list_intsc_v1]
df_viral_v1.to_csv(path_viral_v1_prev, sep="\t")

path_viral_v2_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/viral_v2_prev.txt"
df_viral_v2 = df_viral[list_intsc_v2]
df_viral_v2.to_csv(path_viral_v2_prev, sep="\t")

path_viral_v3_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/viral_v3_prev.txt"
df_viral_v3 = df_viral[list_intsc_v3]
df_viral_v3.to_csv(path_viral_v3_prev, sep="\t")

path_viral_recov_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/viral_recov_prev.txt"
df_viral_recov = df_viral[list_intsc_recov]
df_viral_recov.to_csv(path_viral_recov_prev, sep="\t")

path_anxiety_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/anxiety_prev.txt"
df_anxiety = df_stress[list_intsc_anxiety]
df_anxiety.to_csv(path_anxiety_prev, sep="\t")

path_suicide_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/suicide_attempt_prev.txt"
df_suicide = df_stress[list_intsc_suicide]
df_suicide.to_csv(path_suicide_prev, sep="\t")

path_depression_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression/depression_prev.txt"
df_depression = df_stress[list_intsc_depression]
df_depression.to_csv(path_depression_prev, sep="\t")

# %%
path_meta_healthy_cohort1_train_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/healthy_cohort_train_prev.txt"
df_meta_healthy_cohort1_train = df_meta_healthy.loc[list_intsc_train]
df_meta_healthy_cohort1_train.to_csv(path_meta_healthy_cohort1_train_prev, sep="\t")

path_meta_healthy_cohort1_valid_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/healthy_cohort_valid_prev.txt"
df_meta_healthy_cohort1_valid = df_meta_healthy.loc[list_intsc_valid]
df_meta_healthy_cohort1_valid.to_csv(path_meta_healthy_cohort1_valid_prev, sep="\t")

path_meta_viral_v1_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/viral_v1_prev.txt"
df_meta_viral_v1 = df_meta_viral.loc[list_intsc_v1]
df_meta_viral_v1.to_csv(path_meta_viral_v1_prev, sep="\t")

path_meta_viral_v2_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/viral_v2_prev.txt"
df_meta_viral_v2 = df_meta_viral.loc[list_intsc_v2]
df_meta_viral_v2.to_csv(path_meta_viral_v2_prev, sep="\t")

path_meta_viral_v3_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/viral_v3_prev.txt"
df_meta_viral_v3 = df_meta_viral.loc[list_intsc_v3]
df_meta_viral_v3.to_csv(path_meta_viral_v3_prev, sep="\t")

path_meta_viral_recov_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/viral_recov_prev.txt"
df_meta_viral_recov = df_meta_viral.loc[list_intsc_recov]
df_meta_viral_recov.to_csv(path_meta_viral_recov_prev, sep="\t")

path_meta_anxiety_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/anxiety_prev.txt"
df_meta_anxiety = df_meta_stress.loc[list_intsc_anxiety]
df_meta_anxiety.to_csv(path_meta_anxiety_prev, sep="\t")

path_meta_suicide_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/suicide_attempt_prev.txt"
df_meta_suicide = df_meta_stress.loc[list_intsc_suicide]
df_meta_suicide.to_csv(path_meta_suicide_prev, sep="\t")

path_meta_depression_prev = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata/depression_prev.txt"
df_meta_depression = df_meta_stress.loc[list_intsc_depression]
df_meta_depression.to_csv(path_meta_depression_prev, sep="\t")