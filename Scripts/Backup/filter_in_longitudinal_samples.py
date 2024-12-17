# %%
import os
import pandas as pd

# %%
WORKDIR = "./LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below"
path_healthy_illumina = f"{WORKDIR}/training_dataset.tsv"
path_viral_v1 = f"{WORKDIR}/testing_dataset_selected_viral_v1.tsv"
path_viral_v2 = f"{WORKDIR}/testing_dataset_selected_viral_v2.tsv"
path_viral_v3 = f"{WORKDIR}/testing_dataset_selected_viral_v3.tsv"

# %%
def get_target_samples(path):
    with open (path, mode="r") as fr:
        list_samples = list()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            index = record[0]
            list_samples.append(index)
    
    list_samples = list_samples[1:]

    return list_samples

# %%
list_sample_healthy = get_target_samples(path_healthy_illumina)
list_sample_viral_v1 = get_target_samples(path_viral_v1)
list_sample_viral_v2 = get_target_samples(path_viral_v2)
list_sample_viral_v3 = get_target_samples(path_viral_v3)
set_subj_v1 = set(list(map(lambda x: "-".join(x.split("-")[:-1]), list_sample_viral_v1)))
list_sample_viral_v2_long = list(filter(lambda x: "-".join(str(x).split("-")[:-1]) in set_subj_v1, list_sample_viral_v2))
list_sample_viral_v3_long = list(filter(lambda x: "-".join(str(x).split("-")[:-1]) in set_subj_v1, list_sample_viral_v3))

# %%
path_meta = "./Metadata/Metadata_RNA.tsv"
df_meta = pd.read_csv(path_meta, sep="\t")
df_meta_set_idx = df_meta.set_index("Project-ID")
df_meta_healthy = df_meta_set_idx.loc[list_sample_healthy, :]
df_meta_healthy.to_csv("./Metadata/healthy_train.txt", sep="\t")

df_meta_v1 = df_meta_set_idx.loc[list_sample_viral_v1, :]
df_meta_v1.to_csv("./Metadata/viral_v1_long.txt", sep="\t")

df_meta_v2 = df_meta_set_idx.loc[list_sample_viral_v2_long, :]
df_meta_v2.to_csv("./Metadata/viral_v2_long.txt", sep="\t")

df_meta_v3 = df_meta_set_idx.loc[list_sample_viral_v3_long, :]
df_meta_v3.to_csv("./Metadata/viral_v3_long.txt", sep="\t")


# %%
path_exp = "./Expression/deseq2_normalized_count.tsv"
df_exp = pd.read_csv(path_exp, sep="\t")
df_exp_reset_idx = df_exp.reset_index(drop=False).rename(columns={"index": "ID"})
df_exp_healthy = df_exp_reset_idx[["ID", *list_sample_healthy]]
df_exp_v1 = df_exp_reset_idx[["ID", *list_sample_viral_v1]]
df_exp_v1.to_csv("./Expression/viral_v1_long.txt", sep="\t", index=False)
df_exp_v2 = df_exp_reset_idx[["ID", *list_sample_viral_v2_long]]
df_exp_v2.to_csv("./Expression/viral_v2_long.txt", sep="\t", index=False)
df_exp_v3 = df_exp_reset_idx[["ID", *list_sample_viral_v3_long]]
df_exp_v3.to_csv("./Expression/viral_v3_long.txt", sep="\t", index=False)
# %%
