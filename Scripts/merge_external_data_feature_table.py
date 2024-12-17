# %%
import os

import numpy as np
import pandas as pd

from utility import Utility


# %%
def get_metadata_geo(dir_metadata, GEO_ID):
    path_meta_geo = os.path.join(dir_metadata, f"{GEO_ID}.txt")
    df_meta_geo = pd.read_csv(path_meta_geo, sep="\t")

    return df_meta_geo

def get_all_samples(df_meta_geo, col_sample="Project-ID"):
    samples_all_geo = df_meta_geo[col_sample].to_list()

    return samples_all_geo

def get_metadata_geo_healthy(df_meta_geo, col_phenotype="Sample_characteristics_ch1", list_tag_healthy=["healthy", "control"]):
    df_meta_geo_healthy = df_meta_geo[df_meta_geo[col_phenotype].str.lower().str.split(":").str[-1].str.lstrip().isin(list_tag_healthy)]

    return df_meta_geo_healthy

def get_samples_healthy(df_meta_geo_healthy, col_sample="Project-ID"):
    samples_healthy_geo = df_meta_geo_healthy[col_sample].to_list()

    return samples_healthy_geo

def get_samples_healthy_selected(df_meta_geo_healthy, samples_healthy_geo, ratio):
    np.random.seed(42)
    samples_selected_healthy_geo = np.random.choice(samples_healthy_geo, size=int(len(df_meta_geo_healthy)*ratio), replace=False)

    return samples_selected_healthy_geo

def get_samples_unselected_geo(samples_all_geo, samples_selected_healthy_geo):
    samples_unselected_geo = list(set(samples_all_geo).difference(set(samples_selected_healthy_geo)))

    return samples_unselected_geo

def save_metadata_test(dir_metadata, GEO_ID, samples_unselected_geo, col_sample="Project-ID"):
    meta_geo_test = os.path.join(dir_metadata, f"{GEO_ID}_test.txt")
    df_meta_geo_set_index = df_meta_geo.set_index(col_sample)
    df_meta_geo_test = df_meta_geo_set_index.loc[samples_unselected_geo]
    df_meta_geo_test.to_csv(meta_geo_test, sep="\t", index=True)

def save_metadata_train(dir_metadata, GEO_ID, samples_selected_geo, col_sample="Project-ID"):
    meta_geo_train = os.path.join(dir_metadata, f"{GEO_ID}_train.txt")
    df_meta_geo_set_index = df_meta_geo.set_index(col_sample)
    df_meta_geo_train = df_meta_geo_set_index.loc[samples_selected_geo]
    df_meta_geo_train.to_csv(meta_geo_train, sep="\t", index=True)

def get_feature_table(dir_feat_geo, geo_id, col_sample="Project-ID"):
    path_feat_geo = os.path.join(dir_feat_geo, f"{geo_id}.txt")
    df_feat_geo = pd.read_csv(path_feat_geo, sep="\t").set_index(col_sample)

    return df_feat_geo

def get_feature_table_geo_test(df_feat_geo, samples_unselected_geo):
    df_feat_geo_test = df_feat_geo.loc[samples_unselected_geo]

    return df_feat_geo_test

def get_feature_table_geo_train(df_feat_geo, samples_selected_geo):
    df_feat_geo_train = df_feat_geo.loc[samples_selected_geo]

    return df_feat_geo_train

def save_feature_table_geo_test(dir_feature_table, GEO_ID, df_feat_geo_test):
    feat_geo_test = os.path.join(dir_feature_table, f"{GEO_ID}_test.txt")
    df_feat_geo_test.to_csv(feat_geo_test, sep="\t", index=True)
    
def save_feature_table_geo_train(dir_feature_table, GEO_ID, df_feat_geo_train):
    feat_geo_train = os.path.join(dir_feature_table, f"{GEO_ID}_test.txt")
    df_feat_geo_train.to_csv(feat_geo_train, sep="\t", index=True)
    
def get_genes_included_feature_table(df_feat_geo):
    genes_geo = list(filter(lambda x: "ENSG" in x, list(df_feat_geo.columns)))

    return genes_geo

def get_genes_intersect(list_genes_geo):
    genes_intersect = list(set(list_genes_geo[0]).intersection(*list_genes_geo[1:]))

    return genes_intersect

def get_feature_table_for_merge(path_feat_geo, list_samples_filter_in, col_y, col_sample, list_target_genes=list()):
    util = Utility()
    X, ys = util.get_feature_table(path_feat_geo)
    X_, y = util.filter_input_data(X, ys, col_y)
    X_set_index = X_.set_index(col_sample)
    if list_target_genes != list():
        X_selected_set_idx = util.select_dataset_target_genes(X_set_index, list_target_genes, col_sample)
    else:
        X_selected_set_idx = X_set_index
    X_selected = X_selected_set_idx.reset_index(drop=False)
    df_feature = pd.concat([X_selected, y], axis=1)
    df_feature_set_index = df_feature.set_index(col_sample)
    if list_samples_filter_in is None:
        df_feature_filtered = df_feature_set_index
    else:
        df_feature_filtered = df_feature_set_index.loc[list_samples_filter_in]

    return df_feature_filtered

def merge_feature_tables(list_dfs_feature_filtered):
    df_feat = pd.concat(list_dfs_feature_filtered, axis=0)

    return df_feat

# %%
col_y = "Sample_Age"
col_sample="Project-ID"
WORKDIR = ".."
dir_metadata = f"{WORKDIR}/Metadata"
dir_feature_table = f"{WORKDIR}/FeatureTable"
list_geo_id = ["GSE134080", "GSE270454"]
train_name = "healthy_cohort_train_prev"

list_genes_geo = list()
for geo_id in list_geo_id:
    df_feature_geo = get_feature_table(dir_feature_table, geo_id)
    genes_geo = get_genes_included_feature_table(df_feature_geo)
    list_genes_geo.append(genes_geo)

genes_intersect = get_genes_intersect(list_genes_geo)

list_df_feat_for_merge = list()
for geo_id in list_geo_id:
    path_feat_geo = os.path.join(dir_feature_table, f"{geo_id}.txt")
    df_meta_geo = get_metadata_geo(dir_metadata, geo_id)
    samples_all_geo = get_all_samples(df_meta_geo)
    df_meta_geo_healthy = get_metadata_geo_healthy(df_meta_geo)
    samples_healthy_geo = get_samples_healthy(df_meta_geo_healthy)
    samples_selected_healthy_geo = get_samples_healthy_selected(df_meta_geo_healthy, samples_healthy_geo, ratio=0.9)
    df_feature_geo = get_feature_table(dir_feature_table, geo_id)    
    samples_unselected_geo = get_samples_unselected_geo(samples_all_geo, samples_selected_healthy_geo)
    save_metadata_test(dir_metadata, geo_id, samples_unselected_geo)
    save_metadata_train(dir_metadata, geo_id, samples_selected_healthy_geo)
    df_feature_geo_test = get_feature_table_geo_test(df_feature_geo, samples_unselected_geo)
    df_feature_geo_train = get_feature_table_geo_train(df_feature_geo, samples_selected_healthy_geo)
    save_feature_table_geo_test(dir_feature_table, geo_id, df_feature_geo_test)
    save_feature_table_geo_train(dir_feature_table, geo_id, df_feature_geo_train)
    df_feat_for_merge = get_feature_table_for_merge(path_feat_geo, samples_selected_healthy_geo, col_y, col_sample, list_target_genes=genes_intersect).reset_index(drop=False)
    list_df_feat_for_merge.append(df_feat_for_merge)

# %%
geo_names = "_".join(list_geo_id)
filename = f"{train_name}_{geo_names}_mixed_ethnicity.txt"
df_feat_merged = merge_feature_tables(list_df_feat_for_merge).set_index("Project-ID")
df_feat_original = get_feature_table(dir_feature_table, train_name)
list_targetd_genes = list(set(list(df_feat_merged.columns)).intersection(set(list(df_feat_original.columns))))
list_targetd_genes.remove(col_y)
list_new_columns = list_targetd_genes + [col_y]
df_feat_merged = df_feat_merged[list_new_columns]
df_feat_original = df_feat_original[list_new_columns]
df_feat_original_merged = pd.concat([df_feat_original, df_feat_merged], axis=0)
feat_new = os.path.join(dir_feature_table, filename)
df_feat_original_merged.to_csv(feat_new, sep="\t", index=True)

# %%
# dir_metadata = f"{WORKDIR}/Metadata"
# dir_feature_table = f"{WORKDIR}/FeatureTable"
# geo_id = "GSE177044"
# col_phenotype="Project-ID"
# tag_healthy="Control"
# geo_id = "GSE134080"
# col_phenotype="Sample_characteristics_ch1"
# tag_healthy="healthy"

# path_feat_geo = os.path.join(dir_feature_table, f"{geo_id}.txt")
# df_meta_geo = get_metadata_geo(dir_metadata, geo_id)
# samples_all_geo = get_all_samples(df_meta_geo)
# df_meta_geo_healthy = get_metadata_geo_healthy(df_meta_geo, col_phenotype=col_phenotype, tag_healthy=tag_healthy)
# samples_healthy_geo = get_samples_healthy(df_meta_geo_healthy)
# samples_selected_healthy_geo = get_samples_healthy_selected(df_meta_geo_healthy, samples_healthy_geo, ratio=0.9)
# df_feature_geo = get_feature_table(dir_feature_table, geo_id)    
# samples_unselected_geo = get_samples_unselected_geo(samples_all_geo, samples_selected_healthy_geo)
# save_metadata_test(dir_metadata, geo_id, samples_unselected_geo)
# df_feature_geo_test = get_feature_table_geo_test(df_feature_geo, samples_unselected_geo)
# save_feature_table_geo_test(dir_feature_table, geo_id, df_feature_geo_test)
# df_feat_for_merge = get_feature_table_for_merge(path_feat_geo, samples_selected_healthy_geo, col_y, col_sample).reset_index(drop=False)
# df_feat_for_merge.to_csv(os.path.join(dir_feature_table, f"{geo_id}_train.txt"), sep="\t", index=False)
# save_metadata_train(dir_metadata, geo_id, samples_selected_healthy_geo)

# %%
