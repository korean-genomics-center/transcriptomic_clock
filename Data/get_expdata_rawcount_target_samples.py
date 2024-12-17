#%%
import os
import pandas as pd

# %%
path_exp = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/KOGIC_STAR_RSEM_Freeze_V1_RawCount_For_Norm.gct"
df_exp = pd.read_csv(path_exp, sep="\t", index_col=[0])
list_genes = list(map(lambda x: "_".join(x.split("_")[1:]), list(df_exp.index)))
df_exp.index = list_genes
list_columns = list(df_exp.columns)

# %%
dir_sample_list = "/BiO/Access/kyungwhan1998/transcriptomic_clock/Data"
dir_exp = dir_sample_list
# list_cohorts = ["healthy_valid_prev", "healthy_cohort2"]
# list_filenames = ["healthy_cohort_valid_prev", "healthy_cohort2"]
list_cohorts = ["healthy_cohort2"]
list_filenames = ["healthy_cohort2"]
for cohort, filename in zip(list_cohorts, list_filenames):
    path_sample_list = os.path.join(dir_sample_list, f"sample_list_{cohort}.txt")
    with open(path_sample_list, mode="r") as fr:
        list_samples = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
        list_intsc_columns = list(set(list_columns).intersection(set(list_samples)))
        df_exp_target = df_exp[list_intsc_columns]
        filename_exp = f"{filename}_rawcount.csv"
        df_exp_target.to_csv(os.path.join(dir_exp, filename_exp), sep=",", index=True)