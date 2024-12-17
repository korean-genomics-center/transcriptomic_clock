# %%
import os
import pandas as pd

# %%
list_cohorts = ["anxiety", "attempt", "depression", "v1", "v2", "v3", "viral_recov"]
for cohort in list_cohorts:
    path_exp = f"/BiO/Access/kyungwhan1998/transcriptomic_clock/LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below/sample_statistics_{cohort}.tsv"
    df_exp = pd.read_csv(path_exp, sep="\t")
    list_subjects = df_exp["SUBJID"].to_list()
    
    path_sample_list = f"/BiO/Research/GeroExpressome/Resources/Miscellaneous/sample_list_{cohort}_prev.txt"
    with open(path_sample_list, mode="w") as fw:
        for subj in list_subjects:
            fw.write(subj + "\n")