# %%
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

# %%
def get_list_file_qc(dir_qc, dict_qc, pattern="**/*.tsv"):
    list_all_file_qc = glob.glob(f"{dir_qc}/{pattern}", recursive = True)

    list_file_qc = list()
    for file_qc in list_all_file_qc:
        filename_qc = os.path.basename(file_qc)
        if filename_qc.split(".")[0] in dict_qc.keys():
            list_file_qc.append(file_qc)

    return list_file_qc

def get_id_samples(file_qc_stat, mode="fastp"):
    dict_id = dict()

    with open(file_qc_stat, mode="r") as fr:

        list_stat_id = list()
        id_header = fr.readline().rstrip("\n").split("\t")[0]
        for line in fr:
            record = line.rstrip("\n").split("\t")
            stat_id = record[0]
            list_stat_id.append(stat_id)

        list_lane_id = list()
        for stat_id, cnt in dict(Counter(list_stat_id)).items():
            for lane_num in range(1, int(cnt)+1, 1):
                lane_id = f"{stat_id}_L0{lane_num}"
                list_lane_id.append(lane_id)
        
        if mode == "fastp":
            dict_id[id_header] = list_lane_id

        elif mode == "star":
            dict_id[id_header] = list_stat_id

        else:
            print("unsupported mode")

    return dict_id

def remove_id_samples(df_stat_merge, list_samples_excl):
    df_stat_merge_sample_excl = df_stat_merge[~df_stat_merge["SampleID"].str.split("_").str[0].isin(list_samples_excl)]
    
    return df_stat_merge_sample_excl

def modify_id_samples(list_old_id, mode="fastp"):
    list_modif_id = list(map(lambda x: x.replace("U10K", "KU10K") if str(x).startswith("U10K") else x, list_old_id))
    if mode == "fastp":
        list_modif_id_1 = list(map(lambda x: x.replace("Ax-", "AX").replace("AX-", "AX").replace("_L01_L01", "_L01").replace("_L02_L01", "_L02").replace("_S218_L004_R1", "").replace("_S203_L004_R1", "").replace("-A1", ""), list_modif_id))

    elif mode == "star":
        list_modif_id_1 = list(map(lambda x: x.replace("Ax-", "AX").replace("AX-", "AX").replace("_L01", "").replace("_L02", "").replace("_S218_L004_R1", "").replace("_S203_L004_R1", "").replace("-A1", ""), list_modif_id))
    
    else:
        print("unsupported mode")

    return list_modif_id_1

def get_dict_id_conversion(conversion_table):
    dict_conversion = dict()
    with open(conversion_table, mode="r") as fr:
        fr.readline()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            new_id = record[0]
            old_id = record[1]
            dict_conversion[old_id] = new_id
    
    return dict_conversion

def get_list_converted_id(list_modif_id, dict_conversion, mode="fastp", verbose=True):

    if mode == "fastp":
        list_converted_id = list(map(lambda x: dict_conversion[x.split("_")[0]] + "_" + x.split("_")[1] if x.split("_")[0] in dict_conversion else x, list_modif_id))

    elif mode == "star":
        list_converted_id = list(map(lambda x: dict_conversion[x] if x in dict_conversion else x, list_modif_id))
    
    else:
        print("unsupported mode")

    if verbose:
        for modif, conv in zip(list_modif_id, list_converted_id):
            print(modif, conv)

    return list_converted_id

def get_dict_statistics(file_qc_stat, dict_qc):
    dict_stats = dict()
    with open(file_qc_stat, mode="r") as fr:
        stat_header = fr.readline().rstrip("\n").split("\t")[-1]
        stat_newheader = dict_qc[stat_header]
        
        list_stat_content = list()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            stat_content = (record[-1])
            list_stat_content.append(stat_content)
        
        dict_stats[stat_newheader] = list_stat_content

    return dict_stats

def get_list_samples_expression_matrix(path_exp):
    df_exp = pd.read_csv(path_exp, sep="\t")
    list_samples_exp = list(df_exp.columns)[1:]

    return list_samples_exp

def filter_qc_exp_samples(df_stat_merge_reordered, list_samples_exp, mode="fastp"):
    if mode == "fastp":
        df_stat_merge_reordered["tmp"] = df_stat_merge_reordered["SampleID"].apply(lambda x: x.split("_")[0])
        df_stat_merge_reordered_filtered = df_stat_merge_reordered[df_stat_merge_reordered["tmp"].isin(list_samples_exp)]
        df_stat_merge_reordered_filtered = df_stat_merge_reordered_filtered.drop(columns=["tmp"])

    if mode == "star":
        df_stat_merge_reordered_filtered = df_stat_merge_reordered[df_stat_merge_reordered["SampleID"].isin(list_samples_exp)]

    return df_stat_merge_reordered_filtered

def draw_histogram_stat(df_qc, stats):
    df_qc[stats]= df_qc[stats].astype(float)
    # cnt, bins = np.histogram(df_qc[stats])
    # start = math.floor(df_qc[stats].min() / bin_size) * bin_size
    # end = math.ceil(df_qc[stats].max() / bin_size) * bin_size
    plt.hist(df_qc[stats], zorder=1)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(stats, fontsize=16)
    plt.ylabel("Count", fontsize=16)
    plt.grid(axis="both", zorder=0)
    plt.tight_layout()
    # plt.savefig(outfile, dpi=600)
    plt.show()
    plt.close()

# %% [fastp qc - KU10K]
mode = "fastp"
project = "ku10k"
dir_qc = "/BiO/Store/KOGIC/RNASeq/QCCollection/KU10K/StatsSummary"
dict_fastp_qc = {
    "summary___before_filtering___total_bases": "Total_Number_of_Bases_Before_Filtering",
    "summary___after_filtering___total_bases": "Total_Number_of_Bases_After_Filtering",
    "summary___before_filtering___read1_mean_length": "Read1_Mean_Length_Before_Filtering",
    "summary___before_filtering___read2_mean_length": "Read2_Mean_Length_Before_Filtering",
    "summary___after_filtering___read1_mean_length": "Read1_Mean_Length_After_Filtering",
    "summary___after_filtering___read2_mean_length": "Read2_Mean_Length_After_Filtering",
    "summary___before_filtering___q30_rate": "Q30_Rate_Before_Filtering",
    "summary___after_filtering___q30_rate": "Q30_Rate_After_Filtering",
    "summary___before_filtering___gc_content": "GC_Rate_Before_Filtering",
    "summary___after_filtering___gc_content": "GC_Rate_After_Filtering",
}
conversion_table = "/BiO/Store/KOGIC/RNASeq/QCCollection/conversion_table.txt"
path_exp =  "/BiO/Store/KOGIC/RNASeq/ReleaseData/20240225/expression_matrix_genes.results_expected_count.tsv"
list_samples_excl = ["KU10K-00922-A1", "KU10K-01025-A1", "KU10K-01027-A1", "KU10K-02029-A1"]

list_file_qc = get_list_file_qc(dir_qc, dict_fastp_qc)
first_file = list_file_qc[0]
dict_id = get_id_samples(first_file)

dict_stat_merge = dict(dict_id)
for file_qc in list_file_qc:
    dict_stat = get_dict_statistics(file_qc, dict_fastp_qc)
    dict_stat_merge.update(dict_stat)
 
df_stat_merge = pd.DataFrame(dict_stat_merge)
df_stat_merge_sample_excl = remove_id_samples(df_stat_merge, list_samples_excl)
new_colheader = ["SampleID"] + list(dict_fastp_qc.values())
df_stat_merge_reordered = df_stat_merge_sample_excl[new_colheader]
list_old_id = df_stat_merge_reordered["SampleID"]
list_new_id = modify_id_samples(list_old_id)
dict_conversion = get_dict_id_conversion(conversion_table)
list_converted_id = get_list_converted_id(list_new_id, dict_conversion, verbose=True)
df_stat_merge_reordered["SampleID"] = list_converted_id
list_samples_exp = get_list_samples_expression_matrix(path_exp)
df_stat_merge_reordered_filtered = filter_qc_exp_samples(df_stat_merge_reordered, list_samples_exp, mode=mode)
df_stat_merge_reordered_filtered["SubjectID"] = df_stat_merge_reordered_filtered["SampleID"].str.split("_").str[0]
df_stat_merge_reordered_filtered.to_csv(f"/BiO/Store/KOGIC/RNASeq/QCCollection/{mode}_qc_{project}.tsv", sep="\t", index=False)

for stats in list(dict_fastp_qc.values()):
    draw_histogram_stat(df_stat_merge_reordered_filtered, stats=stats)

# %% [star qc - KU10K]
mode="star"
project="ku10k"
dir_qc = "/BiO/Store/KOGIC/RNASeq/QCCollection/KU10K/StatsSummary"
dict_star_qc = {
    "Total_mapped_reads_number": "Total_Mapping_Read_Number",
    "Total_mapped_reads_percentage": "Total_Mapping_Rate"
}
conversion_table = "/BiO/Store/KOGIC/RNASeq/QCCollection/conversion_table.txt"
path_exp =  "/BiO/Store/KOGIC/RNASeq/ReleaseData/20240225/expression_matrix_genes.results_expected_count.tsv"
list_samples_excl = ["KU10K-00922-A1", "KU10K-01025-A1", "KU10K-01027-A1", "KU10K-02029-A1"]

list_file_qc = get_list_file_qc(dir_qc, dict_star_qc)
first_file = list_file_qc[0]
dict_id = get_id_samples(first_file, mode=mode)

dict_stat_merge = dict(dict_id)
for file_qc in list_file_qc:
    dict_stat = get_dict_statistics(file_qc, dict_star_qc)
    dict_stat_merge.update(dict_stat)
 
df_stat_merge = pd.DataFrame(dict_stat_merge)
df_stat_merge_sample_excl = remove_id_samples(df_stat_merge, list_samples_excl)
new_colheader = ["SampleID"] + list(dict_star_qc.values())
df_stat_merge_reordered = df_stat_merge_sample_excl[new_colheader]
list_old_id = df_stat_merge_reordered["SampleID"]
list_new_id = modify_id_samples(list_old_id, mode=mode)
dict_conversion = get_dict_id_conversion(conversion_table)
list_converted_id = get_list_converted_id(list_new_id, dict_conversion, mode=mode, verbose=True)
df_stat_merge_reordered["SampleID"] = list_converted_id
list_samples_exp = get_list_samples_expression_matrix(path_exp)
df_stat_merge_reordered_filtered = filter_qc_exp_samples(df_stat_merge_reordered, list_samples_exp, mode=mode)
df_stat_merge_reordered_filtered.to_csv(f"/BiO/Store/KOGIC/RNASeq/QCCollection/{mode}_qc_{project}.tsv", sep="\t", index=False)

for stats in list(dict_star_qc.values()):
    draw_histogram_stat(df_stat_merge_reordered_filtered, stats=stats)

# %% [fastp qc - COVID19]
mode = "fastp"
dir_qc = "/BiO/Store/KOGIC/RNASeq/QCCollection/COVID19"
dict_fastp_qc = {
    "summary___before_filtering___total_bases": "Total_Number_of_Bases_Before_Filtering",
    "summary___after_filtering___total_bases": "Total_Number_of_Bases_After_Filtering",
    "summary___before_filtering___read1_mean_length": "Read1_Mean_Length_Before_Filtering",
    "summary___before_filtering___read2_mean_length": "Read2_Mean_Length_Before_Filtering",
    "summary___after_filtering___read1_mean_length": "Read1_Mean_Length_After_Filtering",
    "summary___after_filtering___read2_mean_length": "Read2_Mean_Length_After_Filtering",
    "summary___before_filtering___q30_rate": "Q30_Rate_Before_Filtering",
    "summary___after_filtering___q30_rate": "Q30_Rate_After_Filtering",
    "summary___before_filtering___gc_content": "GC_Rate_Before_Filtering",
    "summary___after_filtering___gc_content": "GC_Rate_After_Filtering",
}
list_samples_excl = []

list_file_qc_first = get_list_file_qc(dir_qc, dict_fastp_qc, pattern="1*/**/*.tsv")
list_file_qc_second = get_list_file_qc(dir_qc, dict_fastp_qc, pattern="2*/**/*.tsv")

first_file = list_file_qc_first[0]
second_file = list_file_qc_second[0]

dict_id_first = get_id_samples(first_file)
dict_id_second = get_id_samples(second_file)

dict_stat_merge_first = dict(dict_id_first)
for file_qc in list_file_qc_first:
    dict_stat_first = get_dict_statistics(file_qc, dict_fastp_qc)
    dict_stat_merge_first.update(dict_stat_first)

dict_stat_merge_second = dict(dict_id_second)
for file_qc in list_file_qc_second:
    dict_stat_second = get_dict_statistics(file_qc, dict_fastp_qc)
    dict_stat_merge_second.update(dict_stat_second)

df_stat_merge_first = pd.DataFrame(dict_stat_merge_first)
new_colheader = ["SampleID"] + list(dict_fastp_qc.values())
df_stat_merge_first_reordered = df_stat_merge_first[new_colheader]

df_stat_merge_second = pd.DataFrame(dict_stat_merge_second)
new_colheader = ["SampleID"] + list(dict_fastp_qc.values())
df_stat_merge_second_reordered = df_stat_merge_second[new_colheader]

df_fastp_qc = pd.concat([df_stat_merge_first_reordered, df_stat_merge_second_reordered], axis=0)
df_fastp_qc_sorted = df_fastp_qc.sort_values(by=["SampleID"], ascending=True).reset_index(drop=True)
df_fastp_qc_sorted["SubjectID"] = df_fastp_qc_sorted["SampleID"].str.split("_").str[0]
df_fastp_qc_sorted.to_csv(f"/BiO/Store/KOGIC/RNASeq/QCCollection/{mode}_qc_covid19.tsv", sep="\t", index=False)

for stats in list(dict_fastp_qc.values()):
    draw_histogram_stat(df_fastp_qc_sorted, stats=stats)

# %% [star qc - COVID19]
mode="star"
dir_qc = "/BiO/Store/KOGIC/RNASeq/QCCollection/COVID19"
dict_star_qc = {
    "Total_mapped_reads_number": "Total_Mapping_Read_Number",
    "Total_mapped_reads_percentage": "Total_Mapping_Rate"
}
list_samples_excl = []

list_file_qc_first = get_list_file_qc(dir_qc, dict_star_qc, pattern="1*/**/*.tsv")
list_file_qc_second = get_list_file_qc(dir_qc, dict_star_qc, pattern="2*/**/*.tsv")

first_file = list_file_qc_first[0]
second_file = list_file_qc_second[0]

dict_id_first = get_id_samples(first_file, mode=mode)
dict_id_second = get_id_samples(second_file, mode=mode)

dict_stat_merge_first = dict(dict_id_first)
for file_qc in list_file_qc_first:
    dict_stat_first = get_dict_statistics(file_qc, dict_star_qc)
    dict_stat_merge_first.update(dict_stat_first)

dict_stat_merge_second = dict(dict_id_second)
for file_qc in list_file_qc_second:
    dict_stat_second = get_dict_statistics(file_qc, dict_star_qc)
    dict_stat_merge_second.update(dict_stat_second)

df_stat_merge_first = pd.DataFrame(dict_stat_merge_first)
new_colheader = ["SampleID"] + list(dict_star_qc.values())
df_stat_merge_first_reordered = df_stat_merge_first[new_colheader]

df_stat_merge_second = pd.DataFrame(dict_stat_merge_second)
new_colheader = ["SampleID"] + list(dict_star_qc.values())
df_stat_merge_second_reordered = df_stat_merge_second[new_colheader]

df_star_qc = pd.concat([df_stat_merge_first_reordered, df_stat_merge_second_reordered], axis=0)
df_star_qc_sorted = df_star_qc.sort_values(by=["SampleID"], ascending=True).reset_index(drop=True)
df_star_qc_sorted.to_csv(f"/BiO/Store/KOGIC/RNASeq/QCCollection/{mode}_qc_covid19.tsv", sep="\t", index=False)

for stats in list(dict_star_qc.values()):
    draw_histogram_stat(df_star_qc_sorted, stats=stats)
# %%
