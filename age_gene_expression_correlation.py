import math
import os
from multiprocessing import Pool
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import fdrcorrection
from utility import Utility

def train_linear_regression_model(array_variable, array_gene_exp) -> dict:
    x = np.array(array_variable).astype(np.float64)
    x_reshape = x.reshape(-1,1)
    y = np.array(array_gene_exp).astype(np.float64)
    model = LinearRegression().fit(x_reshape, y)

    return model

def multiprocess_training(array_array_variable: np.ndarray, array_array_gene_exp: np.ndarray, num_threads=10) -> list:
    with Pool(processes=num_threads) as pool:
        input_arrays = zip(array_array_variable, array_array_gene_exp)
        models = pool.starmap(train_linear_regression_model, input_arrays)
    
    return models

def save_linear_regression_expression_models(models: object, outdir: str, intermediatefilename: str):
    util = Utility()
    util.dump_pickle(models, outdir, intermediatefilename)

def load_linear_regression_expression_models(outdir: str, intermediatefilename: str):
    util = Utility()
    models = util.load_pickle(outdir, intermediatefilename)

    return models

def test_linear_regression_model(model: object, array_variable:np.ndarray, array_gene_exp: np.ndarray) -> dict:
    x = array_variable
    y = array_gene_exp
    x_reshape = x.reshape(-1, 1)
    slope = model.coef_[0]
    rsq = model.score(x_reshape, y)
    y_pred = model.predict(x_reshape)

    corr, pval = pearsonr(x, y)
    median = np.median(y)
    dict_reg_results = {"slope":slope, "rsq":rsq, "corr":corr, "pval":pval, "median": median, "y_pred":y_pred}

    return dict_reg_results

def multiprocess_testing(models:object, array_array_variable: np.ndarray, array_array_gene_exp: np.ndarray, num_threads=10) -> list:
    with Pool(processes=num_threads) as pool:
        input_model_and_arrays = zip(models, array_array_variable, array_array_gene_exp)
        list_dict_linear_regression_expression = pool.starmap(test_linear_regression_model, input_model_and_arrays)
    
    return list_dict_linear_regression_expression

def correct_pval_multitest(df_linear_reg, pval_thres=0.05, pval_direct="below"):
    qval_thres = pval_thres
    df_linear_reg_copy = df_linear_reg.copy()
    pvals = df_linear_reg_copy["pval"].to_list()
    sigs, pvals = fdrcorrection(pvals=pvals, alpha=qval_thres)
    df_linear_reg_copy["testsig"] = sigs
    df_linear_reg_copy["pval"] = pvals
    
    if pval_direct == "below":
        df_linear_reg_corrected = df_linear_reg_copy[df_linear_reg_copy["testsig"] == True]
    elif pval_direct == "above":
        df_linear_reg_corrected = df_linear_reg_copy[df_linear_reg_copy["testsig"] == False]
    else:
        print(f"unsupported mode: {pval_direct}")

    return df_linear_reg_corrected
    
def make_df_linear_regression_results(list_dict_res: list, list_gene_name: list, corr_thres, pval_thres, corr_direct="above", pval_direct="below", pval_correction=False) -> pd.DataFrame:
    df_linear_reg = pd.DataFrame.from_records(list_dict_res, index=list_gene_name)
    df_linear_reg = df_linear_reg.reset_index(drop=False).rename(columns={"index": "gene_id"})
    df_linear_reg = df_linear_reg.drop(columns=["y_pred"])
    df_linear_reg["abs_corr"] = df_linear_reg["corr"].apply(abs)
    if corr_direct == "above":
        df_linear_reg_filt = df_linear_reg[abs(df_linear_reg["corr"]) > corr_thres]
        df_linear_reg_filt = df_linear_reg_filt.sort_values(by=["abs_corr"], ascending=False)
    elif corr_direct == "below":
        df_linear_reg_filt = df_linear_reg[abs(df_linear_reg["corr"]) < corr_thres]
        df_linear_reg_filt = df_linear_reg_filt.sort_values(by=["abs_corr"], ascending=True)
    else:
        print(f"unsupported mode: {corr_direct}")
    if pval_correction:
        df_linear_reg_filt = correct_pval_multitest(df_linear_reg_filt, pval_thres, pval_direct)
    else:
        if pval_direct == "below":
            df_linear_reg_filt = df_linear_reg_filt[df_linear_reg_filt["pval"] < pval_thres]
        elif pval_direct == "above":
            df_linear_reg_filt = df_linear_reg_filt[df_linear_reg_filt["pval"] > pval_thres]
        else:
            print(f"unsupported mode: {pval_direct}")

    return df_linear_reg_filt

def write_linear_regression_statistics(dir_linear_reg:str, df_linear_reg: pd.DataFrame, varname: str, outfilename: str) -> None:
    outfilepath = f"{dir_linear_reg}/linear_regression_expression_{varname}_{outfilename}.tsv"
    outfiledir = os.path.dirname(outfilepath)
    os.makedirs(outfiledir, exist_ok=True)

    df_linear_reg.to_csv(outfilepath, sep="\t", index=False)

    return None

def plot_linear_regression_results(array_variable: np.array, list_array_exp_raw: list, list_array_exp_pred: list, df_linear_reg: pd.DataFrame, varname: str, outdir: str, corr_thres: float, pval_thres: float, outfilename: str, corr_direct="above" , pval_direct="below", pval_correction=False) -> None:    
    df_linear_reg = df_linear_reg.head(25)
    list_gene_id = df_linear_reg["gene_id"].to_list()
    list_array_rsq = df_linear_reg["rsq"].to_list()
    list_array_corr = df_linear_reg["corr"].to_list()
    list_array_pval = df_linear_reg["pval"].to_list()
    num_genes = len(list_gene_id)

    plt.rcParams["font.size"] = 14
    fig, axes = plt.subplots(math.ceil(math.sqrt(num_genes)), math.ceil(math.sqrt(num_genes)), figsize=(20, 20), sharey=False, sharex=True)
    list_contents = zip(list_gene_id, list_array_exp_raw, list_array_exp_pred, list_array_rsq, list_array_corr, list_array_pval, axes.flatten())
    for gene_id, array_exp, array_exp_pred, rsq, corr, pval, ax in list_contents:
        ax.scatter(array_variable, array_exp, color="white", edgecolor="grey")
        ax.plot(array_variable, array_exp_pred, color="red", linewidth=3)
        plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
        plt.setp(ax.set_title(f"{'_'.join(gene_id.split('_')[1:])}\n({gene_id.split('_')[0]})"), fontsize=plt.rcParams["font.size"]+4)
        xticks = list(range(10, 101, 10))
        plt.setp(ax.set_xticks(xticks))
        plt.setp(ax.set_xticklabels(xticks, rotation=90))

        if pval_correction:
            plt.setp(ax.annotate(f"$R^2$:{round(float(rsq),3)}\nr:{round(float(corr),3)}\nq:{round(float(pval),3)}", xy=(0.10, 0.70), xycoords='axes fraction'), weight="bold", fontsize=15)
        else:
            plt.setp(ax.annotate(f"$R^2$:{round(float(rsq),3)}\nr:{round(float(corr),3)}\np:{round(float(pval),3)}", xy=(0.10, 0.70), xycoords='axes fraction'), weight="bold", fontsize=15)


    if corr_direct == "above":
        corr_ineq = ">"
    elif corr_direct == "below":
        corr_ineq = "<"
    else: 
        print("Unsupported mode")
    
    if pval_direct == "above":
        pval_ineq = ">"
    elif pval_direct == "below":
        pval_ineq = "<"
    else:
        print("Unsupported mode")

    fig.subplots_adjust(top=0.95, bottom=0.08, left=0.10, right=0.95, hspace=0.4, wspace=0.3)

    # Set figure title and axis labels manually
    fig.suptitle(f"Age-Expression Correlation ($r_{{gene}}$ {corr_ineq} {str(corr_thres)}; $q$ {pval_ineq} {str(pval_thres)}; N = {str(len(array_variable))})", y=0.98, fontsize=30)
    fig.supylabel("Normalized Expression", fontsize=25)
    fig.supxlabel("Chronological Age", fontsize=25)

    outfigpath = f"{outdir}/linear_regression_expression_{varname}_{outfilename}.png"
    plt.tight_layout()
    plt.savefig(outfigpath, dpi=600)
    plt.show()
    plt.close()