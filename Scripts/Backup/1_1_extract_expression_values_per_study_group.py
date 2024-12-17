# %%
import os
import numpy as np
from utility import Utility

# %% [initialize]
WORKDIR = ".."
util = Utility()
path_exp_mtx = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/KOGIC_STAR_RSEM_Freeze_V1_DESeq2_NormCount.gct"
path_metadata = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/KOGIC_Metadata_Freeze_V1.txt"
dir_exp_out = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Expression"
dir_meta_out = "/BiO/Research/GeroExpressome/Workspace/RNAClock/Revision_final_trained_with_only_korean/Metadata"
file_rare_disease_samples = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/rare_disease_sample_list.txt"
file_bgi_samples = "/BiO/Research/GeroExpressome/Resources/Data/KOGIC/bgi500_2019_sample_list.txt"
os.makedirs(dir_exp_out, exist_ok=True)
os.makedirs(dir_meta_out, exist_ok=True)
list_drop_samples = util.get_list_drop_samples(file_rare_disease_samples)
list_samples_bgi = util.get_list_drop_samples(file_bgi_samples)
list_samples_PC_outlier = ["KU10K-02988", 
                           "KU10K-02959", 
                           "KU10K-02765", 
                           "KU10K-02972", 
                           "KU10K-03628", 
                           "KU10K-03745"]
colsample = "Project-ID"
colgene = "Gene_ID"

# %% [all]
df_meta = util.read_metadata(path_metadata, colsample=colsample)
list_select_samples = df_meta[colsample].to_list()
print(len(list_select_samples))
# %% [healthy]
df_meta_healthy = df_meta[np.logical_and(df_meta["Sample_Trait"]=="Healthy", df_meta["Sample_Trait_Verbose"] != "Vaccination")]
df_meta_healthy_incl = util.drop_samples_meta(df_meta_healthy, list_drop_samples, colsample=colsample)
df_meta_healthy_incl.to_csv(os.path.join(dir_meta_out, "healthy.txt"), sep="\t", index=False)
list_healthy = df_meta_healthy_incl[colsample].to_list()
print(len(list_healthy))
# %% [healthy cohort1]
df_meta_healthy = df_meta[np.logical_and(df_meta["Sample_Trait"]=="Healthy", df_meta["Sample_Trait_Verbose"] != "Vaccination")]
df_meta_healthy_incl = util.drop_samples_meta(df_meta_healthy, list_drop_samples, colsample=colsample)
df_meta_healthy_cohort1 = util.drop_samples_meta(df_meta_healthy_incl, list_samples_bgi, colsample=colsample)
df_meta_healthy_cohort1 = util.drop_samples_meta(df_meta_healthy_cohort1, list_samples_PC_outlier, colsample=colsample)
df_meta_healthy_cohort1.to_csv(os.path.join(dir_meta_out, "healthy_cohort1.txt"), sep="\t", index=False)
list_healthy_cohort1 = df_meta_healthy_cohort1[colsample].to_list()
print(len(list_healthy_cohort1))
# %% [healthy cohort2]
df_meta_healthy = df_meta[np.logical_and(df_meta["Sample_Trait"]=="Healthy", df_meta["Sample_Trait_Verbose"] != "Vaccination")]
df_meta_healthy_incl = util.drop_samples_meta(df_meta_healthy, list_drop_samples, colsample=colsample)
df_meta_healthy_cohort2 = util.select_samples_meta(df_meta_healthy_incl, list_samples_bgi, colsample=colsample)
df_meta_healthy_cohort2 = util.drop_samples_meta(df_meta_healthy_cohort2, list_samples_PC_outlier, colsample=colsample)
df_meta_healthy_cohort2.to_csv(os.path.join(dir_meta_out, "healthy_cohort2.txt"), sep="\t", index=False)
list_healthy_cohort2 = df_meta_healthy_cohort2[colsample].to_list()
print(len(list_healthy_cohort2))
# %% [stress]
df_meta_stress = df_meta[df_meta["Sample_Trait"] == "MentalDisorder"]
df_meta_stress.to_csv(os.path.join(dir_meta_out, "stress.txt"), sep="\t", index=False)
list_mental_samples = df_meta_stress[colsample].to_list()
print(len(list_mental_samples))
# %% [suicide]suicide
df_meta_suicide = df_meta[df_meta["Sample_Trait_Verbose"] == "SuicideAttempter"]
df_meta_suicide.to_csv(os.path.join(dir_meta_out, "suicide_attempt.txt"), sep="\t", index=False)
list_suicide_samples = df_meta_suicide[colsample].to_list()
print(len(list_suicide_samples))
# %% [anxiety]
df_meta_anxiety = df_meta[df_meta["Sample_Trait_Verbose"] == "Anxiety"]
df_meta_anxiety.to_csv(os.path.join(dir_meta_out, "anxiety.txt"), sep="\t", index=False)
list_anxiety_samples = df_meta_anxiety[colsample].to_list()
print(len(list_anxiety_samples))
# %% [depression]
df_meta_depress = df_meta[df_meta["Sample_Trait_Verbose"] == "Depression"]
df_meta_depress.to_csv(os.path.join(dir_meta_out, "depression.txt"), sep="\t", index=False)
list_depress_samples = df_meta_depress[colsample].to_list()
print(len(list_depress_samples))
# %% [cardio]
df_meta_cardio = df_meta[df_meta["Sample_Trait"] == "CardioVascularDisease"]
df_meta_cardio.to_csv(os.path.join(dir_meta_out, "cardio.txt"), sep="\t", index=False)
list_cardio_samples = df_meta_cardio[colsample].to_list()
print(len(list_cardio_samples))
# %% [covid19]
df_meta_viral = df_meta[df_meta["Project_Name"]=="COVID19"]
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral.txt"), sep="\t", index=False)
list_viral_samples = df_meta_viral[colsample].to_list()
print(len(list_viral_samples))
# %% [covid19 v1]
df_meta_viral = df_meta[np.logical_and(df_meta[colsample].str.contains("-V1"), df_meta[colsample].str.contains("C19-C"))]
df_meta_viral = util.drop_samples_meta(df_meta_viral, list_samples_PC_outlier, colsample=colsample)
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_v1.txt"), sep="\t", index=False)
list_viral_v1_samples = df_meta_viral[colsample].to_list()
print(len(list_viral_v1_samples))
# %% [covid19 v2]
df_meta_viral = df_meta[np.logical_and(df_meta[colsample].str.contains("-V2"), df_meta[colsample].str.contains("C19-C"))]
df_meta_viral = util.drop_samples_meta(df_meta_viral, list_samples_PC_outlier, colsample=colsample)
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_v2.txt"), sep="\t", index=False)
list_viral_v2_samples = df_meta_viral[colsample].to_list()
print(len(list_viral_v2_samples))
# %% [covid19 v3]
df_meta_viral = df_meta[np.logical_and(df_meta[colsample].str.contains("-V3"), df_meta[colsample].str.contains("C19-C"))]
df_meta_viral = util.drop_samples_meta(df_meta_viral, list_samples_PC_outlier, colsample=colsample)
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_v3.txt"), sep="\t", index=False)
list_viral_v3_samples = df_meta_viral[colsample].to_list()
print(len(list_viral_v3_samples))
# %% [covid19 v4]
df_meta_viral = df_meta[np.logical_and(df_meta[colsample].str.contains("-V4"), df_meta[colsample].str.contains("C19-C"))]
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_v4.txt"), sep="\t", index=False)
list_viral_v4_samples = df_meta_viral[colsample].to_list()
print(len(list_viral_v4_samples))
# %% [covid19 recover]
df_meta_viral = df_meta[np.logical_and(df_meta["Project_Name"]=="COVID19", df_meta[colsample].str.contains("C19-R"))]
df_meta_viral = util.drop_samples_meta(df_meta_viral, list_samples_PC_outlier, colsample=colsample)
df_meta_viral.to_csv(os.path.join(dir_meta_out, "viral_recov.txt"), sep="\t", index=False)
list_viral_recov_samples = df_meta_viral[colsample].to_list()
print(len(list_viral_recov_samples))

# %%
df_exp = util.read_expmtx(path_exp_mtx)
if colgene not in list(df_exp.columns):
    df_exp = df_exp.reset_index(drop=False).rename(columns={"index": colgene})
df_exp_drop = util.drop_samples_exp(df_exp, list_drop_samples)
# df_exp_healthy = util.select_samples_df(df_exp_drop, list_healthy, colgene=colgene)
# df_exp_healthy_cohort1 = util.select_samples_df(df_exp_drop, list_healthy_cohort1, colgene=colgene)
# df_exp_healthy_cohort2 = util.select_samples_df(df_exp_drop, list_healthy_cohort2, colgene=colgene)
# df_exp_stress = util.select_samples_df(df_exp_drop, list_mental_samples, colgene=colgene)
# df_exp_suicide = util.select_samples_df(df_exp_drop, list_suicide_samples, colgene=colgene)
# df_exp_anxiety = util.select_samples_df(df_exp_drop, list_anxiety_samples, colgene=colgene)
# df_exp_depress = util.select_samples_df(df_exp_drop, list_depress_samples, colgene=colgene)
df_exp_cardio = util.select_samples_df(df_exp_drop, list_cardio_samples, colgene=colgene)
# # df_exp_rare = select_samples_df(df_exp_drop, list_rare_samples)
# df_exp_viral = util.select_samples_df(df_exp_drop, list_viral_samples, colgene=colgene)
# df_exp_viral_v1 = util.select_samples_df(df_exp_drop, list_viral_v1_samples, colgene=colgene)
# df_exp_viral_v2 = util.select_samples_df(df_exp_drop, list_viral_v2_samples, colgene=colgene)
# df_exp_viral_v3 = util.select_samples_df(df_exp_drop, list_viral_v3_samples, colgene=colgene)
# df_exp_viral_v4 = util.select_samples_df(df_exp_drop, list_viral_v4_samples, colgene=colgene)
# df_exp_viral_recov = util.select_samples_df(df_exp_drop, list_viral_recov_samples, colgene=colgene)

# %%
# util.save_expmtx(df_exp_suicide, os.path.join(dir_exp_out, "suicide_attempt.txt"))
# util.save_expmtx(df_exp_anxiety, os.path.join(dir_exp_out, "anxiety.txt"))
# util.save_expmtx(df_exp_depress, os.path.join(dir_exp_out, "depression.txt"))
# util.save_expmtx(df_exp_healthy, os.path.join(dir_exp_out, "healthy.txt"))
# util.save_expmtx(df_exp_healthy_cohort1, os.path.join(dir_exp_out, "healthy_cohort1.txt"))
# util.save_expmtx(df_exp_healthy_cohort2, os.path.join(dir_exp_out, "healthy_cohort2.txt"))
# util.save_expmtx(df_exp_stress, os.path.join(dir_exp_out, "stress.txt"))
util.save_expmtx(df_exp_cardio, os.path.join(dir_exp_out, "cardio.txt"))
# util.save_expmtx(df_exp_viral, os.path.join(dir_exp_out, "viral.txt"))
# util.save_expmtx(df_exp_viral_v1, os.path.join(dir_exp_out, "viral_v1.txt"))
# util.save_expmtx(df_exp_viral_v2, os.path.join(dir_exp_out, "viral_v2.txt"))
# util.save_expmtx(df_exp_viral_v3, os.path.join(dir_exp_out, "viral_v3.txt"))
# util.save_expmtx(df_exp_viral_v4, os.path.join(dir_exp_out, "viral_v4.txt"))
# util.save_expmtx(df_exp_viral_recov, os.path.join(dir_exp_out, "viral_recov.txt"))
# %%