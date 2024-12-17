# %%
import os
import numpy as np
import pandas as pd

# %%
dir_scaled_exp = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/LASSO_INFO/healthy_cohort_train_prev_median_filter/corr_0.35_above_pval_0.05_below"
dir_metadata = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata"

# %%
GEO = "GSE107993"
colage = "Sample_Age"
colsample = "Project-ID"
coltrait = "group"

path_scaled_exp = f"{dir_scaled_exp}/{GEO}_std_scaled.txt"

path_metadata = f"{dir_metadata}/{GEO}.txt"

df_scaled_exp = pd.read_csv(path_scaled_exp, sep="\t")
list_columns = list(df_scaled_exp.columns)
df_scaled_exp = df_scaled_exp.iloc[:, :-1]
df_meta = pd.read_csv(path_metadata, sep="\t")
df_meta_set_idx = df_meta.set_index(colsample)
df_merged = pd.merge(df_scaled_exp, df_meta, how="inner", on=colsample)
list_traits = list(set(df_merged[coltrait].to_list()))

for trait in list_traits:
    path_new_scaled_exp = f"{dir_scaled_exp}/{GEO}_{trait}_std_scaled.txt"
    df_merged_selec_trait = df_merged[df_merged[coltrait] == trait]
    df_merged_selec_trait = df_merged_selec_trait[list_columns]
    df_merged_selec_trait.to_csv(path_new_scaled_exp, sep="\t", index=False)

    path_new_meta = f"{dir_metadata}/{GEO}_{trait}.txt"
    list_samples_selec_trait = df_merged_selec_trait[colsample].to_list()
    df_meta_filt_samples = df_meta_set_idx.loc[list_samples_selec_trait]
    df_meta_filt_samples.to_csv(path_new_meta, sep="\t", index=True)
# %%
GEO = "GSE107994"
colage = "Sample_Age"
colsample = "Project-ID"
coltrait = "group"

path_scaled_exp = f"{dir_scaled_exp}/{GEO}_std_scaled.txt"

path_metadata = f"{dir_metadata}/{GEO}.txt"

df_scaled_exp = pd.read_csv(path_scaled_exp, sep="\t")
list_columns = list(df_scaled_exp.columns)
df_scaled_exp = df_scaled_exp.iloc[:, :-1]
df_meta = pd.read_csv(path_metadata, sep="\t")
df_meta_set_idx = df_meta.set_index(colsample)
df_merged = pd.merge(df_scaled_exp, df_meta, how="inner", on=colsample)
list_traits = list(set(df_merged[coltrait].to_list()))

for trait in list_traits:
    path_new_scaled_exp = f"{dir_scaled_exp}/{GEO}_{trait}_std_scaled.txt"
    df_merged_selec_trait = df_merged[df_merged[coltrait] == trait]
    df_merged_selec_trait = df_merged_selec_trait[list_columns]
    df_merged_selec_trait.to_csv(path_new_scaled_exp, sep="\t", index=False)

    path_new_meta = f"{dir_metadata}/{GEO}_{trait}.txt"
    list_samples_selec_trait = df_merged_selec_trait[colsample].to_list()
    df_meta_filt_samples = df_meta_set_idx.loc[list_samples_selec_trait]
    df_meta_filt_samples.to_csv(path_new_meta, sep="\t", index=True)
# %%
GEO = "GSE119117"
colage = "Sample_Age"
colsample = "Project-ID"
coltrait = "Phase"
colgroup = "Group"

path_scaled_exp = f"{dir_scaled_exp}/{GEO}_std_scaled.txt"

path_metadata = f"{dir_metadata}/{GEO}.txt"

df_scaled_exp = pd.read_csv(path_scaled_exp, sep="\t")
list_columns = list(df_scaled_exp.columns)
df_scaled_exp = df_scaled_exp.iloc[:, :-1]
df_meta = pd.read_csv(path_metadata, sep="\t")
df_meta_set_idx = df_meta.set_index(colsample)
df_merged = pd.merge(df_scaled_exp, df_meta, how="inner", on=colsample)
list_traits = list(set(df_merged[coltrait].to_list()))
list_groups = list(set(df_merged[colgroup].to_list()))

for trait in list_traits:
    for group in list_groups:
        path_new_scaled_exp = f"{dir_scaled_exp}/{GEO}_{'_'.join(trait.split())}_{group}_std_scaled.txt"
        cond_filt = np.logical_and(df_merged[coltrait] == trait, df_merged[colgroup] == group)
        df_merged_selec_trait = df_merged[cond_filt]
        df_merged_selec_trait = df_merged_selec_trait[list_columns]
        df_merged_selec_trait.to_csv(path_new_scaled_exp, sep="\t", index=False)

        path_new_meta = f"{dir_metadata}/{GEO}_{'_'.join(trait.split())}_{group}.txt"
        list_samples_selec_trait = df_merged_selec_trait[colsample].to_list()
        df_meta_filt_samples = df_meta_set_idx.loc[list_samples_selec_trait]
        df_meta_filt_samples.to_csv(path_new_meta, sep="\t", index=True)
# %%
GEO = "GSE160299"
colage = "Sample_Age"
colsample = "Project-ID"
coltrait = "Sample_characteristics_ch1"

path_scaled_exp = f"{dir_scaled_exp}/{GEO}_std_scaled.txt"

path_metadata = f"{dir_metadata}/{GEO}.txt"

df_scaled_exp = pd.read_csv(path_scaled_exp, sep="\t")
list_columns = list(df_scaled_exp.columns)
df_scaled_exp = df_scaled_exp.iloc[:, :-1]
df_meta = pd.read_csv(path_metadata, sep="\t")
df_meta_set_idx = df_meta.set_index(colsample)
df_merged = pd.merge(df_scaled_exp, df_meta, how="inner", on=colsample)
list_traits = list(set(df_merged[coltrait].to_list()))

for trait in list_traits:
    path_new_scaled_exp = f"{dir_scaled_exp}/{GEO}_{trait.split(':')[-1].lstrip()}_std_scaled.txt"
    df_merged_selec_trait = df_merged[df_merged[coltrait] == trait]
    df_merged_selec_trait = df_merged_selec_trait[list_columns]
    df_merged_selec_trait.to_csv(path_new_scaled_exp, sep="\t", index=False)
    
    path_new_meta = f"{dir_metadata}/{GEO}_{trait.split(':')[-1].lstrip()}.txt"
    list_samples_selec_trait = df_merged_selec_trait[colsample].to_list()
    df_meta_filt_samples = df_meta_set_idx.loc[list_samples_selec_trait]
    df_meta_filt_samples.to_csv(path_new_meta, sep="\t", index=True)
# %%
GEO = "GSE273149"
colage = "Sample_Age"
colsample = "Project-ID"
coltrait = "disease status"
coltime = "time"

path_scaled_exp = f"{dir_scaled_exp}/{GEO}_std_scaled.txt"

path_metadata = f"{dir_metadata}/{GEO}.txt"

df_scaled_exp = pd.read_csv(path_scaled_exp, sep="\t")
list_columns = list(df_scaled_exp.columns)
df_scaled_exp = df_scaled_exp.iloc[:, :-1]
df_meta = pd.read_csv(path_metadata, sep="\t")
df_meta_set_idx = df_meta.set_index(colsample)
df_merged = pd.merge(df_scaled_exp, df_meta, how="inner", on=colsample)
list_traits = list(set(df_merged[coltrait].to_list()))
list_time = list(set(df_merged[coltime].to_list()))

for trait in list_traits:
    for time in list_time:
        path_new_scaled_exp = f"{dir_scaled_exp}/{GEO}_{trait}_{'_'.join(time.split())}_std_scaled.txt"
        cond_filt = np.logical_and(df_merged[coltrait] == trait, df_merged[coltime] == time)
        df_merged_selec_trait = df_merged[cond_filt]
        df_merged_selec_trait = df_merged_selec_trait[list_columns]
        df_merged_selec_trait.to_csv(path_new_scaled_exp, sep="\t", index=False)

        path_new_meta = f"{dir_metadata}/{GEO}_{trait}_{'_'.join(time.split())}.txt"
        list_samples_selec_trait = df_merged_selec_trait[colsample].to_list()
        df_meta_filt_samples = df_meta_set_idx.loc[list_samples_selec_trait]
        df_meta_filt_samples.to_csv(path_new_meta, sep="\t", index=True)
# %%
