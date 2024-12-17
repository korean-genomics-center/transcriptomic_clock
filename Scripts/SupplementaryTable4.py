# %%
import os
import pandas as pd

# %%

WORKDIR = ".."
DIR_DEG = os.path.join(WORKDIR, "DEG")
DIR_TABLE = os.path.join(WORKDIR, "Tables")
file_reg_res = f"{WORKDIR}/LASSO_INFO/healthy_cohort_train_prev_median_filter/corr_0.35_above_pval_0.05_below/regression_results.txt"
list_feat = list()
with open(file_reg_res, mode="r") as fr:
    _ = fr.readline()
    for line in fr:
        record = line.rstrip().split()
        feat = record[0]
        list_feat.append(feat)

list_stage = ["Acute", "Mid", "Late"]
list_phase = ["V1", "V2", "V3"]

outfile = os.path.join(DIR_TABLE, "Table_S4.txt")
with open(outfile, mode="w") as fw:
    for phase, stage in zip(list_phase, list_stage):
        file_deg = f"{DIR_DEG}/DEG_Viral_{phase}_healthy.txt"
        with open(file_deg, mode="r") as fr:
            for line in fr:
                record = line.rstrip().split()
                ensembl_gene = record[0]
                gene = ensembl_gene.split(".")[0] + "_" + "_".join(ensembl_gene.split("_")[1:])
                if gene in list_feat:
                    stage_info = [stage]
                    new_record = stage_info + record
                    new_record[1] = gene
                    deg_info = "\t".join(new_record) + "\n"
                    fw.write(deg_info)
