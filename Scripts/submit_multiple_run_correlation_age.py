import os
import subprocess

import numpy as np
from joblib import Parallel, delayed


def run(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, corr_thres, pval_thres, corr_direct, pval_direct, pval_correction):\

    dir_script = os.path.dirname(os.path.realpath(__file__))
    file_script = os.path.join(dir_script, "run_correlation_age_gene_expression.py")
    cmd = f"python {file_script} \
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

    with open(os.path.join(outdir, f"correlation_{outfilename}.log"), mode="w")  as fstdout:
        subprocess.run(cmd, shell=True, stdout=fstdout, stderr=fstdout)

def multiprocess_run(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, list_corr_thres, pval_thres, corr_direct, pval_direct, pval_correction):
    with Parallel(n_jobs=10) as parallel:
        parallel(delayed(run)(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads=num_threads, corr_thres=corr_thres, pval_thres=pval_thres, corr_direct=corr_direct, pval_direct=pval_direct, pval_correction=pval_correction) for corr_thres in list_corr_thres)

def main(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, list_corr_thres, pval_thres, corr_direct, pval_direct, pval_correction):
    multiprocess_run(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, list_corr_thres, pval_thres, corr_direct, pval_direct, pval_correction)

if __name__ == "__main__":
    group = "healthy_cohort_train_prev_added_back_10_80_90s_median_filter"
    WORKDIR = ".."
    dir_feat = f"{WORKDIR}/FeatureTable"
    infeat = f"{dir_feat}/{group}.txt"
    outdir = f"{WORKDIR}/Correlation/{group}"
    os.makedirs(outdir, exist_ok=True)
    outfilename = "train"
    col_y = "Sample_Age"
    save_filename = f"correlation_analysis_intermediate_{outfilename}.pk.gz"
    save_intermediate = True
    num_threads = 30
    list_corr_thres = list(map(lambda x: round(float(x), 2), list(np.arange(0.29, 0.4, 0.01))))
    pval_thres = 0.05
    corr_direct = "above"
    pval_direct = "below"
    pval_correction = True
    main(infeat, outdir, outfilename, col_y, save_filename, save_intermediate, num_threads, list_corr_thres, pval_thres, corr_direct, pval_direct, pval_correction)