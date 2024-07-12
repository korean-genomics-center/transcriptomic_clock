# %%
import os
import numpy as np
from utility import Utility
from transcriptomic_clock import Clock

# %%
list_corr_thres = list(map(lambda x: round(float(x), 2), list(np.arange(0.10, 0.70, 0.01))))
for corr_thres in list_corr_thres:
    try:
        pval_thres = 0.05
        num_threads = 20
        max_alphas = 1000
        opt_mode = "info"
        flag_gene_selection = "model"
        list_selected_genes = []
        WORKDIR = "/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240622"
        dir_trained_model = f"{WORKDIR}/LASSO_INFO/standardized/healthy_illumina/corr_{corr_thres}_above_pval_{pval_thres}_below"
        file_training_data = f"{WORKDIR}/LinearRegression/standardized/healthy_illumina/training_dataset.tsv"
        file_testing_data = f"{WORKDIR}/LinearRegression/standardized/healthy_illumina/testing_dataset_healthy.tsv"
        file_linear_reg = f"{WORKDIR}/LinearRegression/standardized/healthy_illumina/corr_{corr_thres}_above_pval_{pval_thres}_below/linear_regression_expression_Sample_Age_train.tsv"
        os.makedirs(dir_trained_model, exist_ok=True)

        util = Utility()
        X_train, y_train = util.get_feature_table(file_training_data)
        X_test, y_test = util.get_feature_table(file_testing_data)
        clock = Clock(dir_trained_model, file_testing_data, col_sampleid="SUBJID", col_y="Sample_Age")
        print("creating input arrays...")
        if flag_gene_selection == "model":
            list_selected_genes = util.get_list_selected_genes(file_linear_reg)
            util.save_selected_genes(dir_trained_model, list_selected_genes)

        elif flag_gene_selection == "list":
            util.save_selected_genes(dir_trained_model, list_selected_genes)

        else:
            list_selected_genes = list(X_train.columns)[1:]
            util.save_selected_genes(dir_trained_model, list_selected_genes)

        array_array_train_gene_exp = clock.get_array_array_gene_exp(X_train, list_selected_genes)
        array_array_train_variable = clock.get_array_array_variable(y_train)
        array_array_test_gene_exp  = clock.get_array_array_gene_exp(X_test, list_selected_genes)
        array_array_test_variable = clock.get_array_array_variable(y_test)

        print("tuning hyperparameters...")
        if opt_mode == "grad":
            model_tune= clock.tune_hyperparameters_lasso_gradient_decscent(array_array_train_gene_exp, array_array_train_variable, num_threads=num_threads, max_alphas=max_alphas, random_state=1)
            util.dump_pickle(model_tune, dir_trained_model, f"model_tune_{clock.col_y}.pk.gz")
        elif opt_mode == "info":
            alpha_aic_, alpha_bic_ = clock.tune_hyperparameters_lasso_aic_bic(array_array_train_gene_exp, array_array_train_variable, dir_trained_model, verbose=True)
        else:
            print(f"opt_mode {opt_mode} not supported")
                
        print("storing hyperparameters...")
        if opt_mode == "grad":
            dict_hyperparam = clock.get_dictionary_hyperparameters_lasso_gradient_decscent(model_tune, dir_trained_model)
        elif opt_mode == "info":
            # dict_hyperparam = clock.get_dictionary_hyperparameters_lasso_aic_bic(alpha_aic_, "AIC", dir_trained_model)
            dict_hyperparam = clock.get_dictionary_hyperparameters_lasso_aic_bic(alpha_bic_, "BIC", dir_trained_model)
        else:
            print(f"opt_mode {opt_mode} not supported")
        util.dump_pickle(dict_hyperparam, dir_trained_model, f"hyperparameters_{clock.col_y}.pk.gz")

        print("training the model...")
        model_train = clock.fit_model_lasso(array_array_train_gene_exp, array_array_train_variable, dict_hyperparam)
        util.dump_pickle(model_train, dir_trained_model, f"model_train_{clock.col_y}.pk.gz")

        print("done!")

    except Exception as e:
        print(e)
# %%
# if __name__ == "__main__":
#     group = ["healthy_illumina"]
#     list_corr_thres = list(map(lambda x: round(float(x), 2), list(np.arange(0, 0.7, 0.01))))
#     gene_tag = "ENS"
#     flag_gene_selection = "model"
#     num_threads = 20
#     list_selected_genes = []
#     col_y = "Sample_Age"
#     pval_thres = 0.05
#     max_alphas = 500
#     opt_mode = "info"
#     for corr_thres in list_corr_thres:
#         try: 
#             dir_trained_model = f"/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240225/LASSO_INFO/standardized/{group}/corr_{corr_thres}_pval_{pval_thres}"
#             training_data = f"/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240225/LinearRegression/standardized/{group}/training_dataset.tsv"
#             testing_data = f"/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240225/LinearRegression/standardized/{group}/testing_dataset.tsv"
#             linear_reg = f"/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240225/LinearRegression/standardized/{group}/corr_{corr_thres}_pval_{pval_thres}/glm_{col_y}.tsv"
#             os.makedirs(dir_trained_model, exist_ok=True)
            
#             main(training_data, testing_data, gene_tag, flag_gene_selection, linear_reg, list_selected_genes, dir_trained_model, col_y, num_threads, max_alphas, opt_mode)
#         except Exception as e:
#             print(e)
#             print(f"No data for {corr_thres}. Skipped...")

# %%
