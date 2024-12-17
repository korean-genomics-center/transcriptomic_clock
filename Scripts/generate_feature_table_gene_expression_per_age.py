#%%
import os
import warnings
from joblib import Parallel, delayed
from utility import Utility
warnings.filterwarnings("ignore")

# %%
WORKDIR = ".."
DIR_EXP = f"{WORKDIR}/Expression"
DIR_META = f"{WORKDIR}/Metadata"
DIR_FEAT = f"{WORKDIR}/FeatureTable"
os.makedirs(DIR_FEAT, exist_ok=True)
colsample = "Project-ID"
colgene = "Gene_ID"

# %%
util = Utility()

def make_merged_feature_table(group, colsample, colgene):
    filename = f"{group}.txt"
    path_feature = os.path.join(DIR_FEAT, filename)
    path_exp = os.path.join(DIR_EXP, filename)
    path_meta = os.path.join(DIR_META, filename)
    list_sample_id = util.get_list_sample_id(path_exp, colgene)
    dict_gene_exp = util.get_dictionary_gene_expression(path_exp)
    df_feature = util.make_feature_table(dict_gene_exp, list_sample_id, colsample)
    df_meta = util.read_metadata(path_meta, colsample) 
    df_feature_meta_merged = util.merge_feature_table_with_metadata(df_feature, df_meta, "inner", colsample)
    util.save_merged_feature_table(df_feature_meta_merged, path_feature)
    print(f"{group} DONE!")

list_groups = [
    "healthy_cohort2"
    # "healthy_cohort_train_prev_median_filter_healthy_cohort_valid_prev_healthy_cohort2_viral_v1_prev_anxiety_prev_suicide_attempt_prev_depression_prev_pca_input"
    # "viral_v1_long",
    # "viral_v2_long",
    # "viral_v3_long",
    # "healthy_cohort_train_prev",
    # "healthy_cohort_valid_prev",
    # "viral_v1_prev", 
    # "viral_v2_prev", 
    # "viral_v3_prev",
    # "viral_recov_prev", 
    # "suicide_attempt_prev", 
    # "anxiety_prev", 
    # "depression_prev",
    # "cardio",
    # "GSE107993",
    # "GSE107994",
    # "GSE119117",
    # "GSE134080",
    # "GSE160299",
    # "GSE177044",
    # "GSE202625",
    # "GSE221911",
    # "GSE267625",
    # "GSE270454",
    # "GSE273149",
    # "GTEx"
]

# %%
with Parallel(n_jobs=30) as parallel:
    parallel(delayed(make_merged_feature_table)(group, colsample, colgene) for group in list_groups)


# %%
