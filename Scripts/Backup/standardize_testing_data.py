 #%%
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
from joblib import Parallel, delayed
from utility import Utility

# %%
def standardize_testing_data(WORKDIR, train_group, test_group, corr_thres, pval_thres, std=True, reg=True, col_y = "Sample_Age"):
    if std:
        path_train_dataset = f"{WORKDIR}/LinearRegression/unstandardized/{train_group}/training_dataset.tsv"
    else:
        path_train_dataset = f"{WORKDIR}/LinearRegression/standardized/{train_group}/training_dataset.tsv"
    if std:
        path_test_dataset = f"{WORKDIR}/FeatureTable/{test_group}.txt"
    else:
        path_test_dataset = f"{WORKDIR}/LinearRegression/standardized/{train_group}/testing_dataset_{test_group}.tsv"

    outdir = f"{WORKDIR}/LASSO_INFO/standardized/{train_group}/corr_{corr_thres}_above_pval_{pval_thres}_below"
    os.makedirs(outdir, exist_ok=True)

    util = Utility()
    X_train, ys_train = util.get_feature_table(path_train_dataset)
    X_train_, y_train = util.filter_input_data(X_train, ys_train, col_y)
    X_test, ys_test = util.get_feature_table(path_test_dataset)
    X_test_, y_test = util.filter_input_data(X_test, ys_test, col_y)
    if reg:
        file_linear_reg = f"{WORKDIR}/LinearRegression/standardized/{train_group}/corr_{corr_thres}_above_pval_{pval_thres}_below/linear_regression_expression_Sample_Age_train.tsv"
        list_target_genes = util.get_list_selected_genes(file_linear_reg)
        util.save_selected_genes(outdir, list_target_genes)
    else:
        list_target_genes = list(X_train.columns)[1:]

    if std:
        try:
            X_train_select, X_test_select = util.select_dataset_target_genes(X_train_, X_test_, list_target_genes, reset_index=False)
            X_train_std, X_test_std = util.standardize_expression(X_train_select, X_test_select, list_target_genes)
            util.save_input_dataset(X_train_std, X_test_std, y_train, y_test, f"testing_dataset_selected_{test_group}.tsv", outdir, index=False)
        except Exception as e:
            print(e)
    
    else:
        X_train_select, X_test_select = util.select_dataset_target_genes(X_train_, X_test_, list_target_genes, reset_index=True)
        util.save_input_dataset(X_train_select, X_test_select, y_train, y_test, f"testing_dataset_selected_{test_group}.tsv", outdir, index=False)

    return f"{corr_thres} Done!"

# %%
util = Utility()
list_corr_thres = list(map(lambda x: round(float(x), 2), list(np.arange(0.10, 0.70, 0.01))))
pval_thres = 0.05
WORKDIR = "./"
train_group = "healthy_illumina"
# list_test_group = ["healthy_bgi", "stress", "viral_v1", "viral_v2", "viral_v3", "viral_v4", "viral_recov", "suicide_attempt", "anxiety", "depression"]
list_test_group = ["viral_v1_long", "viral_v2_long", "viral_v3_long"]
for test_group in list_test_group:
    print(test_group)
    with Parallel(n_jobs=20) as parallel:
        list_status = parallel(delayed(standardize_testing_data)(WORKDIR, train_group, test_group, corr_thres, pval_thres, std=True, reg=True) for corr_thres in list_corr_thres)
        for status in list_status:
            print(status)

# %%
