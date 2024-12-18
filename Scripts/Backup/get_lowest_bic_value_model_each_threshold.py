# %%
import os
import numpy as np
import pandas as pd

dict_model_alpha_bic = dict()
for i in np.arange(0.30, 0.62, 0.01):
    i = round(float(i), 2)
    try:
        list_alpha = list()
        list_bic = list()
        with open(f"./LASSO_INFO/standardized/healthy_illumina/corr_{i}_above_pval_0.05_below/hyperparam_aic_bic_stat.txt", mode="r") as fr:
            _skipheader = fr.readline()
            for line in fr:
                record = line.rstrip("\n").split("\t")
                alpha = record[1]
                bic = record[-1]
                list_alpha.append(alpha)
                list_bic.append(bic)
        
        dict_alpha_bic = dict(zip(list_alpha, list_bic))
        
        lowest_alpha_bic = {k: v for k, v in dict_alpha_bic.items() if v == min(dict_alpha_bic.values())}
        
        dict_model_alpha_bic[i] = [list(lowest_alpha_bic.keys())[0], list(lowest_alpha_bic.values())[0]]

    except Exception as e:
        print(e)

df = pd.DataFrame.from_dict(dict_model_alpha_bic, orient="index").reset_index(drop=False)
df.columns = ["threshold", "alpha", "BIC"]
df.to_csv("./Supplementary_Tables/Table_S6.txt", sep="\t", index=False)
        
# %%
