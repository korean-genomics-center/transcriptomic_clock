#%%
import argparse
import os
import warnings

import age_gene_expression_correlation as agexpr
import numpy as np
from utility import Utility

warnings.filterwarnings("ignore")

def main(feature_table, dir_linear_reg, outfilename, col_y, intermediatefilename, save_intermediate, num_threads, corr_thres, pval_thres, corr_direct, pval_direct, pval_correction):
    util = Utility()
    X_train, y_train = util.get_feature_table(feature_table)
    list_gene_name = util.get_list_gene_name(X_train)
    N = len(list_gene_name)
    array_array_train_gene_exp = np.array(X_train.T.iloc[1:, :].values).astype(np.float64)
    array_array_train_variable = np.array(np.vstack([np.ravel(y_train[col_y])]*N)).astype(np.float64)

    print("training...")
    if not os.path.exists(os.path.join(dir_linear_reg, intermediatefilename)):
        models = agexpr.multiprocess_training(array_array_train_variable, array_array_train_gene_exp, num_threads=num_threads)
        if save_intermediate:
            print("saving...")
            agexpr.save_linear_regression_expression_models(models, dir_linear_reg, intermediatefilename)

    models = agexpr.load_linear_regression_expression_models(dir_linear_reg, intermediatefilename)

    print("testing...")
    list_dict_linear_regression_expression = agexpr.multiprocess_testing(models, array_array_train_variable, array_array_train_gene_exp, num_threads=num_threads)
    df_linear_reg = agexpr.make_df_linear_regression_results(list_dict_linear_regression_expression, list_gene_name, corr_thres, pval_thres, corr_direct=corr_direct, pval_direct=pval_direct, pval_correction=pval_correction)

    print("writing...")
    agexpr.write_linear_regression_statistics(dir_linear_reg, df_linear_reg, col_y, outfilename)
        
    print("drawing...")
    list_sort_indices = list(df_linear_reg.index)
    list_array_exp_raw = list(array_array_train_gene_exp)
    sorted_list_array_exp_raw = [list_array_exp_raw[i] for i in list_sort_indices]
    list_array_exp_pred = list(map(lambda x: x["y_pred"], list_dict_linear_regression_expression))
    sorted_list_array_exp_pred = [list_array_exp_pred[i] for i in list_sort_indices]
    agexpr.plot_linear_regression_results(array_array_train_variable[0], sorted_list_array_exp_raw, sorted_list_array_exp_pred, df_linear_reg, col_y, dir_linear_reg, corr_thres, pval_thres, outfilename, corr_direct=corr_direct, pval_direct=pval_direct, pval_correction=pval_correction)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--infeat", type=str, default="./FeatureTable/healthy.txt")
    parser.add_argument("--outdir", type=str, default="./Correlation")
    parser.add_argument("--outfilename", type=str, default="test")
    parser.add_argument("--col_y", type=str, default="Sample_Age")
    parser.add_argument("--intermediatefilename", type=str, default="correlation_intermediate.pk.gz")
    parser.add_argument("--save_intermediate", type=bool, default=True)
    parser.add_argument("--num_threads", type=int, default=5)
    parser.add_argument("--corr_thres", type=float, default=0.01)
    parser.add_argument("--pval_thres", type=float, default=0.05)
    parser.add_argument("--corr_direct", type=str, default="above")
    parser.add_argument("--pval_direct", type=str, default="below")
    parser.add_argument("--pval_correction", type=bool, default=True)


    args = parser.parse_args()
    dict_args = vars(args)

    dir_linear_reg = os.path.join(dict_args['outdir'], f"corr_{dict_args['corr_thres']}_{dict_args['corr_direct']}_pval_{dict_args['pval_thres']}_{dict_args['pval_direct']}")
    os.makedirs(dir_linear_reg, exist_ok=True)

    main(
        dict_args["infeat"],
        dir_linear_reg,
        dict_args["outfilename"],
        dict_args["col_y"],
        dict_args["intermediatefilename"],
        dict_args["save_intermediate"],
        dict_args["num_threads"], 
        dict_args["corr_thres"], 
        dict_args["pval_thres"],
        dict_args["corr_direct"],
        dict_args["pval_direct"],
        dict_args["pval_correction"]
        )
