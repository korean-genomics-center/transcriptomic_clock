#%%
import os
import warnings
from collections import Counter

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import PercentFormatter

warnings.filterwarnings("ignore")

# %% [global variable]
# list_groups = ["viral", "viral_v1", "viral_v2", "viral_v3", "viral_v4", "viral_recov", "stress", "suicide_attempt", "anxiety", "depression"]
WORKDIR = "/BiO/Access/kyungwhan1998/transcriptomic_aging_paper"
list_groups = ["viral_v1_long", "viral_v2_long", "viral_v3_long"]
num_samples_selected = 0
DIR_EXP = f"{WORKDIR}/Expression"
DIR_META = f"{WORKDIR}/Metadata"
# %%
def draw_barplot(group, list_age_group, outdir, figname, width=7, color="grey"):
    dict_cnt_age_group = dict(sorted(dict(Counter(list_age_group)).items()))
    list_keys = list(dict_cnt_age_group.keys())
    list_vals = list(dict_cnt_age_group.values())
    sum_vals = sum(list_vals)
    list_props = list(map(lambda x: int(x)/int(sum_vals)*100, list_vals))
    list_props_rounded = list(map(lambda x: round(x, 1), list_props))
    _, ax = plt.subplots(figsize=(5,5),layout="constrained")
    ax.bar(list_keys, list_props_rounded, width=width, color=color, zorder=111)
    for i in ax.containers:
        ax.bar_label(i, size=13)
    plt.xlabel("Age (yrs)", fontsize=18)
    plt.ylabel("Sample Proportion (%)", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(visible=True, zorder=-999)
    plt.title(group)
    figpath = os.path.join(outdir, f"age_distribution_{figname}.png")
    plt.savefig(figpath)
    plt.show()
    plt.close()

def draw_histogram(list_age, color="grey"):
    plt.hist(list_age, color=color, weights=np.ones(len(list_age)) / len(list_age), zorder=111)
    plt.xlabel("Age (yrs)", fontsize=18)
    plt.ylabel("Sample Proportion (%)", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.grid(visible=True, zorder=-999)
    plt.show()
    plt.close()

def read_metadata(metadata):
    df_meta = pd.read_csv(metadata, sep="\t")

    return df_meta

def select_samples_meta(df_meta, list_select_samples, colsample="Project-ID"):
    df_meta_select = df_meta[df_meta[colsample].isin(list_select_samples)]

    return df_meta_select  

def drop_samples_meta(df_meta, list_drop_samples, colsample):
    df_meta_drop = df_meta[~df_meta[colsample].isin(list_drop_samples)]

    return df_meta_drop

def save_metadata(df_save, path_output):
    df_save.to_csv(path_output, sep="\t", index=False)

    return None

def read_expmtx(expmtx):
    df_exp = pd.read_csv(expmtx, sep="\t")
    
    return df_exp

def drop_samples_exp(df_exp, list_drop_samples):
    df_exp_drop = df_exp.drop(columns=list_drop_samples)

    return df_exp_drop

def select_samples_exp(df_exp, list_select_samples, colsample="ID"):
    df_exp_indexed = df_exp.set_index(colsample)
    df_exp_select = df_exp_indexed[list_select_samples]
    df_exp_select_reidx = df_exp_select.reset_index(drop=False)

    return df_exp_select_reidx

def save_expmtx(df_save, path_output):
    df_save.to_csv(path_output, sep="\t", index=False)

    return None
    
def add_df_sample_age_group(df_meta):
    df_meta["Sample_Age_Group"] = df_meta["Sample_Age"].map(lambda x: int(str(x)[0] + "0"))
    
    return df_meta

def remove_df_certain_age_group(df_meta, list_excl_age_group):
    df_meta_filtered_out = df_meta[df_meta["Sample_Age_Group"].isin(list_excl_age_group)]
    
    return df_meta_filtered_out

def get_sampleid_removed_age_group(df_meta_filtered_out):
    list_samples_excl = df_meta_filtered_out["Project-ID"].to_list()
    
    return list_samples_excl

def filter_out_df_excl_samples(df_meta_ori, list_samples_excl):
    df_meta_filtered_in = df_meta_ori[~df_meta_ori["Project-ID"].isin(list_samples_excl)]
    
    return df_meta_filtered_in

def get_list_column_val(df_meta_filtered_in, column="Sample_Age_Group"):
    list_column_val = df_meta_filtered_in[column].to_list()
    
    return list_column_val

def get_dictionary_sample_age_group(list_samples, list_age_group):
    dict_samples_age_group = dict(zip(list_samples, list_age_group))
    
    return dict_samples_age_group

def get_randomly_selected_samples_certain_age_group(dict_samples_age_group, age_group_selected=10, num_samples_selected=0, random_seed=1):
    dict_samples_age_group_selected = {sample: age_group for sample, age_group in dict_samples_age_group.items() if int(age_group)==age_group_selected}
    np.random.seed(random_seed)
    selected_samples = list(np.random.choice(list(dict_samples_age_group_selected.keys()), size=num_samples_selected, replace=False))
    unselected_samples = list(filter(lambda x: x not in selected_samples, list(dict_samples_age_group_selected)))
    
    return unselected_samples

# %%
# healthy, cardio
# list_excl_age_group = [10, 80, 90]
# viral
# list_excl_age_group = [10, 70, 80, 90]

list_excl_age_group = []
for group in list_groups:
    group = group
    age_group_selected = 0
    num_samples_selected = num_samples_selected
    path_exp_mtx = f"{DIR_EXP}/{group}.txt"
    path_metadata = f"{DIR_META}/{group}.txt"
    dir_exp_out = f"{DIR_EXP}/selected"
    dir_meta_out = f"{DIR_META}/selected"
    os.makedirs(dir_exp_out, exist_ok=True)
    os.makedirs(dir_meta_out, exist_ok=True)

    df_meta = read_metadata(path_metadata)
    df_meta_added_age_group = add_df_sample_age_group(df_meta)
    df_meta_filtered_out = remove_df_certain_age_group(df_meta_added_age_group, list_excl_age_group)
    list_sample_filtered_out = get_sampleid_removed_age_group(df_meta_filtered_out)
    df_meta_filtered_in = filter_out_df_excl_samples(df_meta, list_sample_filtered_out)
    list_age_group = get_list_column_val(df_meta_filtered_in, column="Sample_Age_Group")
    list_age = get_list_column_val(df_meta_filtered_in, column="Sample_Age")
    list_samples = get_list_column_val(df_meta_filtered_in, column="Project-ID")
    draw_barplot(group, list_age_group, dir_meta_out, group)

    dict_samples_age_group = get_dictionary_sample_age_group(list_samples, list_age_group)
    unselected_samples = get_randomly_selected_samples_certain_age_group(dict_samples_age_group, age_group_selected=age_group_selected, num_samples_selected=num_samples_selected)
    df_meta_dbl_filtered_in = filter_out_df_excl_samples(df_meta_filtered_in, unselected_samples)
    list_age_group_modified = get_list_column_val(df_meta_dbl_filtered_in, column="Sample_Age_Group")
    draw_barplot(group, list_age_group_modified, dir_meta_out, f"{group}_modified")

    list_sample_selected = get_list_column_val(df_meta_dbl_filtered_in, column="Project-ID")

    df_meta = read_metadata(path_metadata)
    df_meta_healthy_illumina = select_samples_meta(df_meta, list_sample_selected)
    path_out_meta = os.path.join(dir_meta_out, f"{group}.txt")
    save_metadata(df_meta_healthy_illumina, path_out_meta)


    df_exp = read_expmtx(path_exp_mtx)
    df_exp_healthy_illumina = select_samples_exp(df_exp, list_sample_selected)
    path_out_exp = os.path.join(dir_exp_out, f"{group}.txt")
    save_expmtx(df_exp_healthy_illumina, path_out_exp)

# %%
def get_list_sample_id(inexp: str, colsample="ID") -> str:
    list_sample_id = list()
    with open(inexp, mode="r") as fr:
        sampleid = fr.readline().rstrip("\n").split("\t")
        list_sample_id.extend(sampleid)
    list_sample_id.remove(colsample)

    return list_sample_id

def get_dictionary_sample_age(inmeta: str) -> str:
    dict_sample_age = dict()
    with open(inmeta, mode="r") as fr_meta:
        header = fr_meta.readline()
        idx_sample = header.index("ID")
        idx_age = header.index("AGE")
        for line in fr_meta:
            record = line.rstrip("\n").split("\t")
            sampleid = record[idx_sample]
            age = record[idx_age]

            if len(age.split("-")) > 0:
                age = (int(age.split("-")[0]) + int(age.split("-")[1]) + 1) / 2

            dict_sample_age[sampleid] = age
    
    return dict_sample_age

def get_dictionary_gene_expression(inexp: str) -> dict:
    dict_gene_exp = dict()
    with open(inexp, mode="r") as fr:
        header = fr.readline()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            gene = record[0]
            expval = record[1:]
            dict_gene_exp[gene] = expval
    
    return dict_gene_exp

def get_dictionary_age_expression(dict_gene_exp: dict, dict_sample_age: dict):
    dict_age_exp = dict()
    for sampleexp, exp in dict_gene_exp.items():
        for samplemeta, age in dict_sample_age.items():
            if sampleexp == samplemeta:
                dict_age_exp[age] = list()
                dict_age_exp[age].extend(exp)
    
    return dict_age_exp

def filter_dictionary_gene_name(dict_gene_exp: dict, remove_tag: list) -> dict:
    dict_gene_exp_name_filt = dict()
    for gene, list_exp in dict_gene_exp.items():
        for tag in remove_tag:
            if tag in gene:
                continue
            else:
                dict_gene_exp_name_filt[gene] = list_exp
    
    return dict_gene_exp_name_filt

def filter_dictionary_expression_cut_off(dict_gene_exp: dict, expr_cut_off: float) -> dict:
    dict_gene_exp_tpm_filt = dict()
    for gene, list_exp in dict_gene_exp.items():
        try:
            array_exp = np.array(list(map(float, list_exp)))
            median_exp = np.median(array_exp)
            if (median_exp > expr_cut_off):
                dict_gene_exp_tpm_filt[gene] = array_exp

        except:
            print(f"bad line skipped: {gene}")
            continue

    return dict_gene_exp_tpm_filt

def get_exp_thres_cutoff(dict_gene_exp_filt: dict):
    list_medianexp = list()
    for _, list_exp in dict_gene_exp_filt.items():
        list_exp_float = list(map(float, list_exp))
        medianexp = np.median(list_exp_float)
        list_medianexp.append(medianexp) 
        
    medexpall = np.median(list_medianexp)
    
    return medexpall

def make_feature_table(dict_gene_exp_filt: dict, list_sample_id: list, colsample="SUBJID") -> pd.DataFrame:
    df_feature = pd.DataFrame.from_records(dict_gene_exp_filt)
    df_feature.index = list_sample_id
    df_feature_reset_idx = df_feature.reset_index(drop=False)
    df_feature_reset_idx_rename = df_feature_reset_idx.rename(columns={"index":colsample})
    
    return df_feature_reset_idx_rename

def read_file_metadata(inmeta):
    df_meta = pd.read_csv(inmeta, sep="\t")
    
    return df_meta

def merge_feature_table_with_metadata(df_feature: pd.DataFrame, df_meta: pd.DataFrame, merge_how: str, merge_on: str) -> pd.DataFrame:
    df_feature_meta_merged = pd.merge(df_feature, df_meta, how=merge_how, on=merge_on)

    return df_feature_meta_merged

def save_merged_feature_table(df_save: pd.DataFrame, path_feature:str) -> pd.DataFrame:
    df_save.to_csv(path_feature, sep="\t", index=False)

# %%
flag_filter = False
remove_tag = ["PAR_Y", "CTC-338M12"]
flag_find_median_exp = False
dir_feature = f"{WORKDIR}/FeatureTable"
expr_cut_off = 0
os.makedirs(dir_feature, exist_ok=True)
for group in list_groups:
    featfilename = f"{group}.txt"
    path_feature = os.path.join(dir_feature, featfilename)
    list_sample_id = get_list_sample_id(path_out_exp)
    dict_gene_exp = get_dictionary_gene_expression(path_out_exp)
    if flag_filter:
        dict_gene_name_filt = filter_dictionary_gene_name(dict_gene_exp, remove_tag)
        dict_gene_exp_filt = filter_dictionary_expression_cut_off(dict_gene_name_filt, expr_cut_off)
    else:
        dict_gene_exp_filt = dict_gene_exp

    if flag_find_median_exp:
        medianexp = get_exp_thres_cutoff(dict_gene_exp_filt)
        print(medianexp)
    else:
        df_feature = make_feature_table(dict_gene_exp_filt, list_sample_id)
        df_meta = read_file_metadata(path_out_meta) 
        df_meta = df_meta.rename(columns={"Project-ID": "SUBJID"})
        df_feature_meta_merged = merge_feature_table_with_metadata(df_feature, df_meta, "inner", "SUBJID")
        save_merged_feature_table(df_feature_meta_merged, path_feature)

# %%
