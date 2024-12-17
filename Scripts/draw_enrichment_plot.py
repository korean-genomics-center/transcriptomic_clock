## courtesy of yoonsung (Yoonsung1203)

# %%
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

def adjust_phrase_length(value, length_limit):
    if len(value) <= length_limit:
        return value
    import math
    ind_best = 0
    len_front = -math.inf
    list_values_split = value.split()
    for ind in range(len(list_values_split)+1):
        values_front = list_values_split[:ind]
        values_back = list_values_split[ind:]
        len_front_diff = len(' '.join(values_front)) - length_limit
        if len_front_diff > 0:
            break
        elif len_front_diff > len_front:
            ind_best = ind
            len_front = len_front_diff
    value_front = ' '.join(list_values_split[:ind_best])
    value_back = ' '.join(list_values_split[ind_best:])
    if len(value_back) > length_limit:
        value_back = adjust_phrase_length(value_back, length_limit)
    return value_front + '\n' + value_back


def filter_table_with_fdr_fe(table, n_max_show, num_minimum_gene, col_fdr = "Enrichment FDR", col_fe = "Fold Enrichment"):
    import math
    table = table[table["nGenes"] > num_minimum_gene]
    table["Neg_Log10_FDR"] = table[col_fdr].apply(lambda x : -math.log10(x))
    table_fdr_sig = table[table[col_fdr] < 0.05]
    table_fdr_sig["FDR_sig"] = 1
    table_fdr_nonsig = table[table[col_fdr] >= 0.05]
    table_fdr_nonsig["FDR_sig"] = 0
    

    list_tables_to_sort_fe = list()
    if table_fdr_sig.shape[0] > n_max_show:
        table_fdr_sig = table_fdr_sig.sort_values(by = col_fe, ascending = False)
        list_tables_to_sort_fe.append(table_fdr_sig.iloc[:n_max_show, :])    
    else:
        list_tables_to_sort_fe.append(table_fdr_sig)
        table_fdr_nonsig = table_fdr_nonsig.sort_values(by = col_fdr)
        list_tables_to_sort_fe.append(table_fdr_nonsig.iloc[:min(n_max_show-table_fdr_sig.shape[0], table_fdr_nonsig.shape[0]), :])    
    
    table_joined = pd.concat(list_tables_to_sort_fe)
    table_joined = table_joined.sort_values(by = ["FDR_sig", col_fe], ascending = False).reset_index(drop = True)
    return table_joined