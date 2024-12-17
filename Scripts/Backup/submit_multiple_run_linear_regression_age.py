import os
import subprocess

import numpy as np
from joblib import Parallel, delayed


def run(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, corr_thres, pval_thres, corr_direct, pval_direct, pval_correction):\

    cmd = f"python ./run_linear_regression_age_gene_expression.py \
            --infeat {infeat}\
            --outdir {outdir}\
            --outfilename {outfilename}\
            --col_y {col_y}\
            --intermediatefilename {save_filename}\
            --save_intermediate {save_intermediate}\
            --num_threads {num_threads}\
            --corr_thres {corr_thres}\
            --pval_thres {pval_thres}\
            --corr_direct {corr_direct}\
            --pval_direct {pval_direct}\
            --pval_correction {pval_correction}"

    with open(os.path.join(outdir, "log_file.txt"), mode="w")  as fstdout:
        subprocess.run(cmd, shell=True, stdout=fstdout, stderr=fstdout)

def multiprocess_run(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, list_corr_thres, pval_thres, corr_direct, pval_direct, pval_correction):
    with Parallel(n_jobs=10) as parallel:
        parallel(delayed(run)(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads=num_threads, corr_thres=corr_thres, pval_thres=pval_thres, corr_direct=corr_direct, pval_direct=pval_direct, pval_correction=pval_correction) for corr_thres in list_corr_thres)

def main(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, list_corr_thres, pval_thres, corr_direct, pval_direct, pval_correction):
    multiprocess_run(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, list_corr_thres, pval_thres, corr_direct, pval_direct, pval_correction)

if __name__ == "__main__":
    # infeat = "./LinearRegression/standardized/healthy_illumina/training_dataset.tsv"
    # infeat = "./LinearRegression/standardized/healthy_illumina/testing_dataset_healthy.tsv"
    infeat = "./LinearRegression/standardized/healthy_illumina/testing_dataset_selected_healthy_bgi.tsv"
    outdir = "./LinearRegression/standardized/healthy_illumina"
    os.makedirs(outdir, exist_ok=True)
    # outfilename = "train"
    # outfilename = "test"
    outfilename = "bgi"
    col_y = "Sample_Age"
    save_filename = f"linear_regression_analysis_intermediate_{outfilename}.pk.gz"
    save_intermediate = True
    num_threads = 30
    # list_corr_thres = list(map(lambda x: round(float(x), 2), list(np.arange(0.10, 0.70, 0.01))))
    list_corr_thres = [0]
    # list_corr_thres = [0.05]
    # pval_thres = 0.05
    pval_thres = 1
    corr_direct = "above" 
    # corr_direct = "below"
    pval_direct = "below"
    # pval_direct = "above"
    pval_correction = True
    main(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, list_corr_thres, pval_thres, corr_direct, pval_direct, pval_correction)
