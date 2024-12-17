# %%
import os
import numpy as np
from utility import Utility
from transcriptomic_clock import Clock
from joblib import Parallel, delayed

# %%
def make_transcriptomic_clock(train_group, test_group, col_sample_id, col_y, flag_gene_selection, opt_mode, num_threads, max_alphas, corr_thres, work_dir):
    DIR_CORR = os.path.join(work_dir, "Correlation", train_group)
    file_corr = f"{DIR_CORR}/corr_{corr_thres}_above_pval_{pval_thres}_below/linear_regression_expression_Sample_Age_train.tsv"
    DIR_LASSO = os.path.join(work_dir, "LASSO_INFO", train_group, f"corr_{corr_thres}_above_pval_{pval_thres}_below")
    file_training_data = os.path.join(DIR_LASSO, f"{train_group}_std_scaled.txt")
    file_testing_data = os.path.join(DIR_LASSO, f"{test_group}_std_scaled.txt")
    dir_trained_clock = DIR_LASSO
    os.makedirs(dir_trained_clock, exist_ok=True)
    
    try:
        util = Utility()
        X_train, y_train = util.get_feature_table(file_training_data)
        X_test, y_test = util.get_feature_table(file_testing_data)
        clock = Clock(dir_trained_clock, file_testing_data, col_sampleid=col_sample_id, col_y=col_y)
        print("creating input arrays...")
        if flag_gene_selection == "model":
            list_selected_genes = util.get_list_selected_genes(file_corr)
            util.save_selected_genes(dir_trained_clock, list_selected_genes)

        elif flag_gene_selection == "list":
            util.save_selected_genes(dir_trained_clock, list_selected_genes)

        else:
            list_selected_genes = list(X_train.columns)[1:]
            util.save_selected_genes(dir_trained_clock, list_selected_genes)

        array_array_train_gene_exp = clock.get_array_array_gene_exp(X_train, list_selected_genes)
        array_array_train_variable = clock.get_array_array_variable(y_train)
        array_array_test_gene_exp  = clock.get_array_array_gene_exp(X_test, list_selected_genes)
        array_array_test_variable = clock.get_array_array_variable(y_test)

        print("tuning hyperparameters...")
        if opt_mode == "grad":
            model_tune= clock.tune_hyperparameters_lasso_gradient_decscent(array_array_train_gene_exp, array_array_train_variable, num_threads=num_threads, max_alphas=max_alphas, random_state=1)
            util.dump_pickle(model_tune, dir_trained_clock, f"model_tune_{clock.col_y}.pk.gz")
        elif opt_mode == "info":
            alpha_aic_, alpha_bic_ = clock.tune_hyperparameters_lasso_aic_bic(array_array_train_gene_exp, array_array_train_variable, dir_trained_clock, verbose=True)
        else:
            print(f"opt_mode {opt_mode} not supported")
                
        print("storing hyperparameters...")
        if opt_mode == "grad":
            dict_hyperparam = clock.get_dictionary_hyperparameters_lasso_gradient_decscent(model_tune, dir_trained_clock)
        elif opt_mode == "info":
            # dict_hyperparam = clock.get_dictionary_hyperparameters_lasso_aic_bic(alpha_aic_, "AIC", dir_trained_clock)
            dict_hyperparam = clock.get_dictionary_hyperparameters_lasso_aic_bic(alpha_bic_, "BIC", dir_trained_clock)
        else:
            print(f"opt_mode {opt_mode} not supported")
        util.dump_pickle(dict_hyperparam, dir_trained_clock, f"hyperparameters_{clock.col_y}.pk.gz")

        print("training the model...")
        model_train = clock.fit_model_lasso(array_array_train_gene_exp, array_array_train_variable, dict_hyperparam)
        util.dump_pickle(model_train, dir_trained_clock, f"model_train_{clock.col_y}.pk.gz")

        print("done!")

    except Exception as e:
        print(e)

# %%
list_corr_thres = list(map(lambda x: round(float(x), 2), list(np.arange(0.29, 0.4, 0.01))))
pval_thres = 0.05
num_threads = 20
max_alphas = 1000
opt_mode = "info"
col_sample_id = "Project-ID"
col_y = "Sample_Age"
flag_gene_selection = "model"
WORKDIR = ".."
train_group = "healthy_cohort_train_prev_added_back_10_80_90s_median_filter"
test_group = "healthy_cohort_valid_prev"

with Parallel(n_jobs=20) as parallel:
    parallel(delayed(make_transcriptomic_clock)(train_group, test_group, col_sample_id, col_y, flag_gene_selection, opt_mode, num_threads, max_alphas, corr_thres, WORKDIR) for corr_thres in list_corr_thres)

# %%
