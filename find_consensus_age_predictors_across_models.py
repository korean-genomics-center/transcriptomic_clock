# %%
import os
import numpy as np
import pandas as pd

dict_model_consensus_markers = dict()
for i in np.arange(0.30, 0.36, 0.01):
    i = round(float(i), 2)
    try:
        list_marker = list()
        list_coef = list()
        with open(f"./LASSO_INFO/standardized/healthy_illumina/corr_{i}_above_pval_0.05_below/regression_results.txt", mode="r") as fr:
            _skipheader = fr.readline()
            for line in fr:
                record = line.rstrip("\n").split("\t")
                marker = record[0]
                coef = record[1]
                list_marker.append(marker)
                list_coef.append(coef)
        
        dict_model_consensus_markers[i] = list_marker

    except Exception as e:
        print(e)

print(dict_model_consensus_markers)

from itertools import chain
from collections import Counter
list_markers_all_thres = list(chain(*list(dict_model_consensus_markers.values())))
Counter(list_markers_all_thres)
# %%
