 #%%
import os
import warnings
import numpy as np
warnings.filterwarnings("ignore")
from joblib import Parallel, delayed
from missing import Impute
from utility import Utility


# %%
def standardize_testing_data(train_group, test_group, impute_group, corr_thres, pval_thres, colsample, col_y, workdir, save_train=False):
    DIR_CORR = os.path.join(workdir, "Correlation", train_group)
    DIR_FEAT = os.path.join(workdir, "FeatureTable")
    DIR_LASSO = os.path.join(workdir, "LASSO_INFO", train_group)

    file_corr = f"{DIR_CORR}/corr_{corr_thres}_above_pval_{pval_thres}_below/linear_regression_expression_Sample_Age_train.tsv"
    path_train_dataset = f"{DIR_FEAT}/{train_group}.txt"
    path_test_dataset = f"{DIR_FEAT}/{test_group}.txt"
    outdir = f"{DIR_LASSO}/corr_{corr_thres}_above_pval_{pval_thres}_below"
    os.makedirs(outdir, exist_ok=True)

    util = Utility()
    list_target_genes = util.get_list_selected_genes(file_corr)
    util.save_selected_genes(outdir, list_target_genes)

    X_train, ys_train = util.get_feature_table(path_train_dataset)
    X_train_, y_train = util.filter_input_data(X_train, ys_train, col_y)
    X_train_set_idx = X_train_.set_index(colsample)
    X_train_selected_set_idx = util.select_dataset_target_genes(X_train_set_idx, list_target_genes, colsample)
    
    X_test, ys_test = util.get_feature_table(path_test_dataset)
    X_test_, y_test = util.filter_input_data(X_test, ys_test, col_y)
    X_test_set_idx = X_test_.set_index(colsample)
    X_test_selected_set_idx = util.select_dataset_target_genes(X_test_set_idx, list_target_genes, colsample)

    path_impute_dataset = f"{DIR_FEAT}/{impute_group}.txt"
    X_impute, ys_impute = util.get_feature_table(path_impute_dataset)
    X_impute_, y_impute = util.filter_input_data(X_impute, ys_impute, col_y)
    X_impute_set_idx = X_impute_.set_index(colsample)
    X_test_set_idx = X_test_.set_index(colsample)
    X_impute_selected_set_idx = util.select_dataset_target_genes(X_impute_set_idx, list_target_genes, colsample)
    imputer = Impute(X_impute_selected_set_idx, X_test_selected_set_idx, list_target_genes)
    if imputer.is_necessary_impute():
        X_imputed_set_idx = imputer.knn_method()
        X_imputed_filtered_set_idx = imputer.return_dataframe_filtered_imputed_samples_only(X_imputed_set_idx)
        X_test_selected_set_idx = X_imputed_filtered_set_idx
        
    X_train_reset_idx_std, X_test_reset_idx_std = util.standardize_expression(X_train_selected_set_idx, X_test_selected_set_idx, colsample)

    if save_train:
        util.save_input_dataset(X_train_reset_idx_std, y_train, f"{train_group}_std_scaled.txt", outdir, index=False)
    
    util.save_input_dataset(X_test_reset_idx_std, y_test, f"{test_group}_std_scaled.txt", outdir, index=False)

# %%
list_corr_thres = list(map(lambda x: round(float(x), 2), list(np.arange(0.29, 0.4, 0.01))))

# %%
util = Utility()
WORKDIR = ".."
train_group = "healthy_cohort_train_prev_added_back_10_80_90s_median_filter"
list_test_groups = [
    # "healthy_cohort_train_prev_added_back_10_80_90s_median_filter"
    # "healthy_cohort_train_prev_median_filter_healthy_cohort_valid_prev_healthy_cohort2_viral_v1_prev_anxiety_prev_suicide_attempt_prev_depression_prev_pca_input"
    # "healthy_cohort_train_prev_median_filter_healthy_cohort_valid_prev_healthy_cohort2_viral_v1_prev_viral_v2_prev_viral_v3_prev_viral_recov_prev_anxiety_prev_suicide_attempt_prev_depression_prev_pca_input"
    # "viral_v1_long", 
    # "viral_v2_long", 
    # "viral_v3_long", 
    # "GTEx",
    # "healthy_cohort_valid_prev", 
    "healthy_cohort2", 
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
    # "GSE273149"
]

impute_group = train_group
pval_thres = 0.05
colsample="Project-ID"
col_y = "Sample_Age"

for corr_thres in list_corr_thres:
    with Parallel(n_jobs=30) as parallel:
         parallel(delayed(standardize_testing_data)(train_group, test_group, impute_group, corr_thres, pval_thres, colsample, col_y, WORKDIR, save_train=False) for test_group in list_test_groups)

# %%
util = Utility()
train_group = "healthy_cohort_train_prev_added_back_10_80_90s_median_filter"
test_group = "healthy_cohort_valid_prev"
impute_group = "healthy_cohort_train_prev_added_back_10_80_90s_median_filter"
pval_thres = 0.05
colsample="Project-ID"
col_y = "Sample_Age"
WORKDIR = ".."

for corr_thres in list_corr_thres:
    standardize_testing_data(train_group, test_group, impute_group, corr_thres, pval_thres, colsample, col_y, WORKDIR, save_train=True)