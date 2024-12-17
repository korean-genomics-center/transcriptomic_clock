# %%
import os
import pandas as pd

# %%
WORKDIR = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean"
DIR_LASSO = os.path.join(WORKDIR, "LASSO_INFO/healthy_cohort_train_prev_median_filter/corr_0.35_above_pval_0.05_below")
DIR_META = os.path.join(WORKDIR, "Metadata")
DIR_EXP = os.path.join(WORKDIR, "Expression")
path_viral_v1 = f"{DIR_LASSO}/viral_v1_prev_std_scaled.txt"
path_viral_v2 = f"{DIR_LASSO}/viral_v2_prev_std_scaled.txt"
path_viral_v3 = f"{DIR_LASSO}/viral_v3_prev_std_scaled.txt"
colsample = "Project-ID"
colgene = "Gene_ID"

# %%
def get_target_samples(path, list_samples_excl=[]):
    with open(path, mode="r") as fr:
        list_samples = list()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            index = record[0]
            list_samples.append(index)
    
    list_samples = list_samples[1:]
    list_target_samples = list(set(list_samples).difference(set(list_samples_excl)))

    return list_target_samples

# %%
list_sample_viral_v1 = get_target_samples(path_viral_v1)
list_sample_viral_v2 = get_target_samples(path_viral_v2)
list_sample_viral_v3 = get_target_samples(path_viral_v3)
set_subj_v1 = set(list(map(lambda x: "-".join(x.split("-")[:-1]), list_sample_viral_v1)))
list_sample_viral_v2_long = list(filter(lambda x: "-".join(str(x).split("-")[:-1]) in set_subj_v1, list_sample_viral_v2))
list_sample_viral_v3_long = list(filter(lambda x: "-".join(str(x).split("-")[:-1]) in set_subj_v1, list_sample_viral_v3))

# %%
path_meta = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/KOGIC_Metadata_Freeze_V1.txt"
df_meta = pd.read_csv(path_meta, sep="\t")
df_meta_set_idx = df_meta.set_index(colsample)

df_meta_v1 = df_meta_set_idx.loc[list_sample_viral_v1, :]
df_meta_v1.to_csv(f"{DIR_META}/viral_v1_long.txt", sep="\t")

df_meta_v2 = df_meta_set_idx.loc[list_sample_viral_v2_long, :]
df_meta_v2.to_csv(f"{DIR_META}/viral_v2_long.txt", sep="\t")

df_meta_v3 = df_meta_set_idx.loc[list_sample_viral_v3_long, :]
df_meta_v3.to_csv(f"{DIR_META}/viral_v3_long.txt", sep="\t")

# %%
path_exp = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/KOGIC_STAR_RSEM_Freeze_V1_DESeq2_NormCount.gct"
df_exp = pd.read_csv(path_exp, sep="\t")
df_exp_reset_idx = df_exp.reset_index(drop=False).rename(columns={"index": colgene})
df_exp_v1 = df_exp_reset_idx[[colgene, *list_sample_viral_v1]]
df_exp_v1.to_csv(f"{DIR_EXP}/viral_v1_long.txt", sep="\t", index=False)
df_exp_v2 = df_exp_reset_idx[[colgene, *list_sample_viral_v2_long]]
df_exp_v2.to_csv(f"{DIR_EXP}/viral_v2_long.txt", sep="\t", index=False)
df_exp_v3 = df_exp_reset_idx[[colgene, *list_sample_viral_v3_long]]
df_exp_v3.to_csv(f"{DIR_EXP}/viral_v3_long.txt", sep="\t", index=False)