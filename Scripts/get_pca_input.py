# %%
import os
import pandas as pd
from pcarna import PCARna

# %%
corr_thres = 0.3
pval_thres = 0.05
list_study_groups = ["healthy_train", "healthy", "healthy_bgi", "stress", "viral_v1", "viral_v2", "viral_v3", "viral_recov", "anxiety", "suicide_attempt", "depression"]
study_group = list_study_groups[0]
dir_train_scaled_exp = f"./LASSO_INFO/standardized/healthy_illumina/corr_{corr_thres}_above_pval_{pval_thres}_below"
file_features = os.path.join(dir_train_scaled_exp,"regression_results.txt")

pca = PCARna("", "", list_drop_samples=[], list_select_samples=[])
list_features = pca.get_important_features(file_features)
list_header = ["SUBJID"] + list_features + ["Sample_Age"]
path_train_scaled_exp = f"{dir_train_scaled_exp}/testing_dataset_selected_{study_group}.tsv"

df_train_scaled_exp = pd.read_csv(path_train_scaled_exp, sep="\t")
df_trained_scaled_exp_all_age_genes_only = df_train_scaled_exp[list_header]
list_subj_id = df_trained_scaled_exp_all_age_genes_only["SUBJID"].to_list()

metadata = "Metadata/Metadata_RNA.tsv"
df_meta = pd.read_csv(metadata, sep="\t")
df_meta_set_idx = df_meta.set_index("Project-ID")
df_meta_filtered = df_meta_set_idx.loc[list_subj_id,:]
df_meta_filtered_reidx = df_meta_filtered.reset_index(drop=False)

file_pca_exp_input = os.path.join(dir_train_scaled_exp, "PCA", f"pca_input_exp_{study_group}.tsv")
file_pca_meta_input = os.path.join(dir_train_scaled_exp, "PCA", f"pca_input_meta_{study_group}.tsv")
dir_pca = os.path.dirname(file_pca_exp_input)
os.makedirs(dir_pca, exist_ok=True)

df_trained_scaled_exp_all_age_genes_only.to_csv(file_pca_exp_input, sep="\t", index=False)
df_meta_filtered_reidx.to_csv(file_pca_meta_input, sep="\t", index=False)

# %%
