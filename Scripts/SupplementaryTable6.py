# %%
import glob
import os
import pandas as pd

# %%
WORKDIR = ".."
DIR_DATA = os.path.join(WORKDIR, "Data")
DIR_QC = f"{WORKDIR}/QCCollection"
DIR_TABLE = f"{WORKDIR}/Tables"
fastp_stat_ku10k = f"{DIR_QC}/fastp_qc_ku10k.tsv"
fastp_stat_covid19 = f"{DIR_QC}/fastp_qc_covid19.tsv"
star_stat_ku10k = f"{DIR_QC}/star_qc_ku10k.tsv"
star_stat_covid19 = f"{DIR_QC}/star_qc_covid19.tsv"
list_file_samples = glob.glob(f"{DIR_DATA}/sample_list_*_prev.txt", recursive=True)

# %%
list_samples_all = list()
list_groups_all = list()
for file_samples in list_file_samples:
    with open(file_samples, mode="r") as fr:
        list_samples = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
        group = os.path.basename(file_samples).replace("sample_list_", "").replace("_prev.txt", "")
        list_groups = [group] * len(list_samples)
    list_samples_all.extend(list_samples)
    list_groups_all.extend(list_groups)

dict_sample_group_all = dict(zip(list_samples_all, list_groups_all))
df_sample_group = pd.DataFrame.from_dict(dict_sample_group_all, orient="index").reset_index(drop=False)
df_sample_group.columns = ["SampleID", "Group"]

df_fastp_stat_ku10k = pd.read_csv(fastp_stat_ku10k, sep="\t")
df_fastp_stat_covid19 = pd.read_csv(fastp_stat_covid19, sep="\t")
df_fastp_stat_all = pd.concat([df_fastp_stat_ku10k, df_fastp_stat_covid19], axis=0)
df_star_stat_ku10k = pd.read_csv(star_stat_ku10k, sep="\t")
df_star_stat_covid19 = pd.read_csv(star_stat_covid19, sep="\t")
df_star_stat_all = pd.concat([df_star_stat_ku10k, df_star_stat_covid19], axis=0)

df_fastp_stat_all_in_sample_list = df_fastp_stat_all[df_fastp_stat_all["SubjectID"].isin(list_samples_all)]
df_fastp_stat_all_in_sample_list = df_fastp_stat_all_in_sample_list.rename(columns={"SampleID": "LaneID", "SubjectID": "SampleID"})
df_fastp_stat_all_in_sample_list = df_fastp_stat_all_in_sample_list[["SampleID", "LaneID", "Q30_Rate_After_Filtering", "GC_Rate_After_Filtering"]]
df_star_stat_all_in_sample_list = df_star_stat_all[df_star_stat_all["SampleID"].isin(list_samples_all)]

df_qc_fastp = pd.merge(df_fastp_stat_all_in_sample_list, df_sample_group, on="SampleID", how="inner")
df_qc_fastp.to_csv(f"{DIR_TABLE}/Table_S6.txt", sep="\t", index=False)
df_qc_star = pd.merge(df_star_stat_all_in_sample_list, df_sample_group, on="SampleID", how="inner")
df_qc_star.to_csv(f"{DIR_TABLE}/Table_S7.txt", sep="\t", index=False)
