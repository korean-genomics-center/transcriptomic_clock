# %%
import os
import numpy as np
import pandas as pd

# %%
def read_metadata(metadata):
    df_meta = pd.read_csv(metadata, sep="\t")

    return df_meta

def select_samples_meta(df_meta, list_select_samples, colsample):
    df_meta_select = df_meta[df_meta[colsample].isin(list_select_samples)]

    return df_meta_select  

def drop_samples_meta(df_meta, list_drop_samples, colsample):
    df_meta_drop = df_meta[~df_meta[colsample].isin(list_drop_samples)]

    return df_meta_drop

def read_expmtx(expmtx):
    df_exp = pd.read_csv(expmtx, sep="\t")
    
    return df_exp

def drop_samples_exp(df_exp, list_drop_samples):
    df_exp_drop = df_exp.drop(columns=list_drop_samples)

    return df_exp_drop

def select_samples_df(df_exp, list_select_samples, colsample="ID"):
    df_exp_indexed = df_exp.set_index(colsample)
    df_exp_select = df_exp_indexed[list_select_samples]
    df_exp_select_reidx = df_exp_select.reset_index(drop=False)

    return df_exp_select_reidx

def save_expmtx(df_save, path_output):
    df_save.to_csv(path_output, sep="\t", index=False)

    return None

def get_list_drop_samples(file_drop_list):
    with open(file_drop_list, mode="r") as fr:
        list_drop_samples = fr.readlines()
    list_drop_samples = list(map(lambda x: x.rstrip("\n"), list_drop_samples))
    
    return list_drop_samples

# %% [initialize]
path_exp_mtx = "/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240621/Expression/deseq2_normalized_count.tsv"
path_metadata = "/BiO/Research/RNASeqReference/Workspace/KU10K/Resources/Data/Metadata/KOGIC_Metadata_RNA_20240621.tsv"
dir_exp_out = "/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240621/Expression"
dir_meta_out = "/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240621/Metadata"
file_rare_disease_samples = "/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240621/Script/rare_disease_sample_list.txt"
file_bgi_samples = "/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240621/Script/bgi500_2019_sample_list.txt"
os.makedirs(dir_exp_out, exist_ok=True)
os.makedirs(dir_meta_out, exist_ok=True)
list_drop_samples = get_list_drop_samples(file_rare_disease_samples)
list_samples_bgi = get_list_drop_samples(file_bgi_samples)

# %% [all]
df_meta = read_metadata(path_metadata)
list_select_samples = df_meta["Project-ID"].to_list()
print(len(list_select_samples))
# %% [healthy]
df_meta_healthy = df_meta[np.logical_and(df_meta["Sample_Trait"]=="Healthy", df_meta["Sample_Trait_Verbose"] != "Vaccination")]
df_meta_healthy_incl = drop_samples_meta(df_meta_healthy, list_drop_samples, colsample="Project-ID")
df_meta_healthy_incl.to_csv(os.path.join(dir_meta_out, "healthy.txt"), sep="\t", index=False)
list_healthy_illumina = df_meta_healthy_incl["Project-ID"].to_list()
print(len(list_healthy_illumina))
# %% [healthy illumina]
df_meta_healthy = df_meta[np.logical_and(df_meta["Sample_Trait"]=="Healthy", df_meta["Sample_Trait_Verbose"] != "Vaccination")]
df_meta_healthy_incl = drop_samples_meta(df_meta_healthy, list_drop_samples, colsample="Project-ID")
df_meta_healthy_illumina = drop_samples_meta(df_meta_healthy_incl, list_samples_bgi, colsample="Project-ID")
df_meta_healthy_illumina.to_csv(os.path.join(dir_meta_out, "healthy_illumina.txt"), sep="\t", index=False)
list_healthy_illumina = df_meta_healthy_illumina["Project-ID"].to_list()
print(len(list_healthy_illumina))
# %% [healthy bgi]
df_meta_healthy = df_meta[np.logical_and(df_meta["Sample_Trait"]=="Healthy", df_meta["Sample_Trait_Verbose"] != "Vaccination")]
df_meta_healthy_incl = drop_samples_meta(df_meta_healthy, list_drop_samples, colsample="Project-ID")
df_meta_healthy_bgi = select_samples_meta(df_meta_healthy_incl, list_samples_bgi, colsample="Project-ID")
df_meta_healthy_bgi.to_csv(os.path.join(dir_meta_out, "healthy_bgi.txt"), sep="\t", index=False)
list_healthy_bgi = df_meta_healthy_bgi["Project-ID"].to_list()
print(len(list_healthy_bgi))
# %% [stress]
df_meta_stress = df_meta[df_meta["Sample_Trait"] == "MentalDisorder"]
df_meta_stress.to_csv(os.path.join(dir_meta_out, "stress.txt"), sep="\t", index=False)
list_mental_samples = df_meta_stress["Project-ID"].to_list()
print(len(list_mental_samples))
# %% [suicide]suicide
df_meta_suicide = df_meta[df_meta["Sample_Trait_Verbose"] == "SuicideAttempter"]
df_meta_suicide.to_csv(os.path.join(dir_meta_out, "suicide_attempt.txt"), sep="\t", index=False)
list_suicide_samples = df_meta_suicide["Project-ID"].to_list()
print(len(list_suicide_samples))
# %% [anxiety]
df_meta_anxiety = df_meta[df_meta["Sample_Trait_Verbose"] == "Anxiety"]
df_meta_anxiety.to_csv(os.path.join(dir_meta_out, "anxiety.txt"), sep="\t", index=False)
list_anxiety_samples = df_meta_anxiety["Project-ID"].to_list()
print(len(list_anxiety_samples))
# %% [depression]
df_meta_depress = df_meta[df_meta["Sample_Trait_Verbose"] == "Depression"]
df_meta_depress.to_csv(os.path.join(dir_meta_out, "depression.txt"), sep="\t", index=False)
list_depress_samples = df_meta_depress["Project-ID"].to_list()
print(len(list_depress_samples))
# %% [cardio]
df_meta_cardio = df_meta[df_meta["Sample_Trait"] == "CardioVascularDisease"]
df_meta_cardio.to_csv(os.path.join(dir_meta_out, "cardio.txt"), sep="\t", index=False)
list_cardio_samples = df_meta_cardio["Project-ID"].to_list()
print(len(list_cardio_samples))
# %% [covid19]
df_meta_viral = df_meta[df_meta["Project_Name"]=="COVID19"]
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral.txt"), sep="\t", index=False)
list_viral_samples = df_meta_viral["Project-ID"].to_list()
print(len(list_viral_samples))
# %% [covid19 v1]
df_meta_viral = df_meta[np.logical_and(df_meta["Project-ID"].str.contains("-V1"), df_meta["Project-ID"].str.contains("C19-C"))]
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_v1.txt"), sep="\t", index=False)
list_viral_v1_samples = df_meta_viral["Project-ID"].to_list()
print(len(list_viral_v1_samples))
# %% [covid19 v2]
df_meta_viral = df_meta[np.logical_and(df_meta["Project-ID"].str.contains("-V2"), df_meta["Project-ID"].str.contains("C19-C"))]
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_v2.txt"), sep="\t", index=False)
list_viral_v2_samples = df_meta_viral["Project-ID"].to_list()
print(len(list_viral_v2_samples))
# %% [covid19 v3]
df_meta_viral = df_meta[np.logical_and(df_meta["Project-ID"].str.contains("-V3"), df_meta["Project-ID"].str.contains("C19-C"))]
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_v3.txt"), sep="\t", index=False)
list_viral_v3_samples = df_meta_viral["Project-ID"].to_list()
print(len(list_viral_v3_samples))
# %% [covid19 v4]
df_meta_viral = df_meta[np.logical_and(df_meta["Project-ID"].str.contains("-V4"), df_meta["Project-ID"].str.contains("C19-C"))]
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_v4.txt"), sep="\t", index=False)
list_viral_v4_samples = df_meta_viral["Project-ID"].to_list()
print(len(list_viral_v4_samples))
# %% [covid19 recover]
df_meta_viral = df_meta[np.logical_and(df_meta["Project_Name"]=="COVID19", df_meta["Project-ID"].str.contains("C19-R"))]
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_recov.txt"), sep="\t", index=False)
list_viral_recov_samples = df_meta_viral["Project-ID"].to_list()
print(len(list_viral_recov_samples))
# %%
df_exp = read_expmtx(path_exp_mtx)
df_exp_reidx = df_exp.reset_index(drop=False).rename(columns={"index": "ID"})
df_exp_drop = drop_samples_exp(df_exp_reidx, list_drop_samples)
df_exp_healthy = select_samples_df(df_exp_drop, list_healthy_illumina)
df_exp_healthy_illumina = select_samples_df(df_exp_drop, list_healthy_illumina)
df_exp_healthy_bgi = select_samples_df(df_exp_drop, list_healthy_bgi)
df_exp_stress = select_samples_df(df_exp_drop, list_mental_samples)
df_exp_suicide = select_samples_df(df_exp_drop, list_suicide_samples)
df_exp_anxiety = select_samples_df(df_exp_drop, list_anxiety_samples)
df_exp_depress = select_samples_df(df_exp_drop, list_depress_samples)
df_exp_cardio = select_samples_df(df_exp_drop, list_cardio_samples)
# df_exp_rare = select_samples_df(df_exp_drop, list_rare_samples)
df_exp_viral = select_samples_df(df_exp_drop, list_viral_samples)
df_exp_viral_v1 = select_samples_df(df_exp_drop, list_viral_v1_samples)
df_exp_viral_v2 = select_samples_df(df_exp_drop, list_viral_v2_samples)
df_exp_viral_v3 = select_samples_df(df_exp_drop, list_viral_v3_samples)
df_exp_viral_v4 = select_samples_df(df_exp_drop, list_viral_v4_samples)
df_exp_viral_recov = select_samples_df(df_exp_drop, list_viral_recov_samples)
# %%
save_expmtx(df_exp_suicide, os.path.join(dir_exp_out, "suicide_attempt.txt"))
save_expmtx(df_exp_anxiety, os.path.join(dir_exp_out, "anxiety.txt"))
save_expmtx(df_exp_depress, os.path.join(dir_exp_out, "depression.txt"))
save_expmtx(df_exp_healthy, os.path.join(dir_exp_out, "healthy.txt"))
save_expmtx(df_exp_healthy_illumina, os.path.join(dir_exp_out, "healthy_illumina.txt"))
save_expmtx(df_exp_healthy_bgi, os.path.join(dir_exp_out, "healthy_bgi.txt"))
save_expmtx(df_exp_stress, os.path.join(dir_exp_out, "stress.txt"))
save_expmtx(df_exp_cardio, os.path.join(dir_exp_out, "cardio.txt"))
save_expmtx(df_exp_viral, os.path.join(dir_exp_out, "viral.txt"))
save_expmtx(df_exp_viral_v1, os.path.join(dir_exp_out, "viral_v1.txt"))
save_expmtx(df_exp_viral_v2, os.path.join(dir_exp_out, "viral_v2.txt"))
save_expmtx(df_exp_viral_v3, os.path.join(dir_exp_out, "viral_v3.txt"))
save_expmtx(df_exp_viral_v4, os.path.join(dir_exp_out, "viral_v4.txt"))
save_expmtx(df_exp_viral_recov, os.path.join(dir_exp_out, "viral_recov.txt"))
