# %%
import glob
import os
import pandas as pd
from utility import Utility

# %%
WORKDIR = ".."
dict_conversion = {"readLength": "Read_Length", 
                    "Project_Year": "Sampling_Year", 
                    "Sequencing_Date": "Sequencing_Year",
                    "Sequencing_Type_Platform": "Sequencing_Platform", 
                    "Sample_Trait": "Sample_Phenotype"
                    }

util = Utility()
colgene = "Gene_ID"
colsample = "Project-ID"

# %%
def get_list_pca_input_study_groups(WORKDIR, list_study_groups):
    list_pca_input = glob.glob(f"{WORKDIR}/**/*.txt", recursive=True)
    list_pca_input_study_groups = list(filter(lambda x: os.path.basename(x).replace(".txt", "") in list_study_groups, list_pca_input))
    
    return list_pca_input_study_groups

def get_list_pca_exp_input(list_pca_input_study_groups):
    list_pca_exp_input = sorted(list(filter(lambda x: "Expression" in x, list_pca_input_study_groups)))

    return list_pca_exp_input

def get_list_multiple_dfs_exp(list_pca_exp_input, colgene, tag_gene="ENSG"):
    list_multiple_dfs_exp = list()
    for pca_exp_input in list_pca_exp_input:
        df_exp = pd.read_csv(pca_exp_input, sep="\t", index_col=0)
        df_exp.index.name = colgene
        if not str(df_exp.index[0]).startswith(tag_gene):
            df_exp = df_exp.set_index(colgene)
        list_genes = list(df_exp.index)
        if len(list_genes[0].split(".")) > 1:
            list_genes = list(map(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]), list_genes))
            df_exp.index = list_genes
        
        list_multiple_dfs_exp.append(df_exp)

    return list_multiple_dfs_exp

def get_list_genes_intersection(list_multiple_dfs_exp):
    list_genes_multiple_dfs = list(map(lambda x: x.index, list_multiple_dfs_exp))
    # Perform set intersection across all dataframes
    list_genes_intersection = list(set.intersection(*map(set, list_genes_multiple_dfs)))

    return list_genes_intersection

def get_df_exp_merged(list_multiple_dfs_exp, list_genes_intersection):
    list_multiple_dfs_exp_filtered = list(map(lambda x: x.loc[list_genes_intersection, :], list_multiple_dfs_exp))
    df_exp_merged = pd.concat(list_multiple_dfs_exp_filtered, axis=1)

    return df_exp_merged

def get_list_pca_meta_input(list_pca_input_study_groups):
    list_pca_meta_input = sorted(list(filter(lambda x: "Metadata" in x , list_pca_input_study_groups)))

    return list_pca_meta_input

def get_list_multiple_dfs_meta(list_pca_meta_input, dict_conversion, colsample):
    list_multiple_dfs_meta = list()
    for pca_meta_input in list_pca_meta_input:
        df_meta = pd.read_csv(pca_meta_input, sep="\t")
        df_meta_renamed = df_meta.rename(columns=dict_conversion)
        df_meta_drop = util.drop_samples_meta(df_meta_renamed, list_drop_samples=list(), colsample=colsample)
        df_meta_set_idx = df_meta_drop.set_index(colsample)
        df_meta_set_idx["Sample_Age_Group"] = df_meta_set_idx["Sample_Age"].apply(lambda x: str(x)[0] + "0")
        df_meta_set_idx["Sample_Sex"] = df_meta_set_idx["Sample_Sex"].apply(lambda x: "M" if x == "male" else x)
        df_meta_set_idx["Sample_Sex"] = df_meta_set_idx["Sample_Sex"].apply(lambda x: "F" if x == "female" else x)
        list_meta_columns = list(df_meta_set_idx.columns)
        for meta_col in list_meta_columns:
            if ":" in df_meta_set_idx[meta_col].astype(str).to_list()[0]:
                df_meta_set_idx[meta_col] = df_meta_set_idx[meta_col].apply(lambda x: x.split(":")[-1].lstrip())

        list_multiple_dfs_meta.append(df_meta_set_idx)
    
    return list_multiple_dfs_meta

def get_list_meta_intersection(list_multiple_dfs_meta):
    list_meta_multiple_dfs = list(map(lambda x: x.columns, list_multiple_dfs_meta))
    # Perform set intersection across all metadata columns
    list_meta_intersection = list(set.intersection(*map(set, list_meta_multiple_dfs)))

    return list_meta_intersection

def get_df_meta_merged(list_multiple_dfs_meta, list_meta_intersection):
    list_multiple_dfs_meta_filtered = list(map(lambda x: x.loc[:, list_meta_intersection], list_multiple_dfs_meta))
    df_meta_merged = pd.concat(list_multiple_dfs_meta_filtered, axis=0)

    return df_meta_merged

# %%
# list_study_groups = ["healthy_cohort_train_prev_median_filter", "healthy_cohort_valid_prev", "healthy_cohort2", "viral_v1_prev", "anxiety_prev", "suicide_attempt_prev", "depression_prev"]
list_study_groups = [""]
list_pca_input_study_groups = get_list_pca_input_study_groups(WORKDIR, list_study_groups)
list_pca_exp_input = get_list_pca_exp_input(list_pca_input_study_groups)
list_multiple_dfs_exp = get_list_multiple_dfs_exp(list_pca_exp_input, colgene)
if len(list_multiple_dfs_exp) > 1:
    list_genes_intersection = get_list_genes_intersection(list_multiple_dfs_exp)
else:
    list_genes_intersection = list(map(lambda x: list(x.index), list_multiple_dfs_exp))[0]
df_exp_pca_input = get_df_exp_merged(list_multiple_dfs_exp, list_genes_intersection)
df_exp_pca_input = df_exp_pca_input.sort_index(axis=1)
path_exp_pca_input = os.path.join(WORKDIR, "Expression", f"{'_'.join(list_study_groups)}_pca_input.txt")
df_exp_pca_input.to_csv(path_exp_pca_input, sep="\t", index=True)

list_pca_meta_input = get_list_pca_meta_input(list_pca_input_study_groups)
list_multiple_dfs_meta = get_list_multiple_dfs_meta(list_pca_meta_input, dict_conversion, colsample)
if len(list_multiple_dfs_meta) > 1:
    list_meta_intersection = get_list_meta_intersection(list_multiple_dfs_meta)
else:
    list_meta_intersection = list(map(lambda x: list(x.columns), list_multiple_dfs_meta))[0]
df_meta_pca_input = get_df_meta_merged(list_multiple_dfs_meta, list_meta_intersection)
df_meta_pca_input = df_meta_pca_input.sort_index(axis=0)
path_meta_pca_input = os.path.join(WORKDIR, "Metadata", f"{'_'.join(list_study_groups)}_pca_input.txt")
df_meta_pca_input.to_csv(path_meta_pca_input, sep="\t", index=True)
