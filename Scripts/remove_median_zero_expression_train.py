# %%
import os
from utility import Utility
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
colsample = "Project-ID"
col_y = "Sample_Age"
exp_cut_off = 20
util = Utility()
WORKDIR = ".."
DIR_FEAT = os.path.join(WORKDIR, "FeatureTable")
file_feat = os.path.join(DIR_FEAT, "healthy_cohort_train_prev_added_back_10_80_90s.txt")

df_exp, df_metadata = util.get_feature_table(file_feat)
df_exp_, df_colage = util.filter_input_data(df_exp, df_metadata, col_y)
df_exp_gene_id_indexed = df_exp_.set_index(colsample)
list_genes = list(df_exp_gene_id_indexed.columns)
series_median = df_exp_gene_id_indexed.median(axis=0)

# # %%
# plt.rcParams["font.size"]=14
# df_exp_median = df_exp_gene_id_indexed.median(axis=0)[df_exp_gene_id_indexed.median(axis=0) < 1000]
# plt.hist(df_exp_median, bins=np.linspace(min(df_exp_median), max(df_exp_median), num=1000), color="grey")
# plt.xlim(0, 100)
# plt.ylabel("Number of Genes", fontsize=plt.rcParams["font.size"]+2)
# plt.xlabel("Median Expression", fontsize=plt.rcParams["font.size"]+2)
# plt.axvline(x=20, color="firebrick", linestyle="dotted")
# plt.show()
# plt.close()

# %%
list_genes_median_filter = list(series_median[series_median > exp_cut_off].index)
df_exp_gene_id_indexed_median_filtered = df_exp_gene_id_indexed[list_genes_median_filter]
df_exp_gene_id_indexed_median_filtered = df_exp_gene_id_indexed_median_filtered[list_genes_median_filter]
df_exp_gene_id_indexed_median_filtered = df_exp_gene_id_indexed_median_filtered.reset_index(drop=False)
# series_non_zero_median_exp = df_exp_gene_id_indexed_median_filtered.median(axis=0)
# median_non_zero_median_exp = np.median(df_exp_gene_id_indexed_median_filtered.median(axis=0))
# list_genes_new_median_filter = list(series_non_zero_median_exp[series_non_zero_median_exp > median_non_zero_median_exp].index)
# df_exp_gene_id_indexed_dbl_median_filtered = df_exp_gene_id_indexed_median_filtered[list_genes_new_median_filter]
# df_exp_gene_id_indexed_dbl_median_filtered = df_exp_gene_id_indexed_dbl_median_filtered.reset_index(drop=False)

# %%
df_feat_median_filtered = pd.concat([df_exp_gene_id_indexed_median_filtered, df_colage], axis=1)
file_out = os.path.join(DIR_FEAT, f"{file_feat.replace('.txt', '')}_median_filter.txt")

df_feat_median_filtered.to_csv(file_out, sep="\t", index=False)
