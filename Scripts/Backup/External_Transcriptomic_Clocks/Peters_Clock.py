# %%
class RNAData():
    def __init__(self, path_exp, path_meta, colgene = "Project-ID", colsample="Project-ID", colage="Sample_Age", colsex="Sample_Sex"):
        from utility import Utility
        util = Utility()
        self.path_exp = path_exp
        self.path_meta = path_meta
        df_rna = util.read_expmtx(self.path_exp)
        # df_rna_column_drop = df_rna.drop(columns=["Description"])
        df_rna[colgene] = df_rna[colgene].apply(lambda x: x.split(".")[0])
        df_rna_dup_row_drop = df_rna.drop_duplicates(subset=[colgene])
        self.rna = df_rna_dup_row_drop.set_index(colgene)
        df_meta_not_selected = util.read_metadata(self.path_meta, colsample=colsample)
        df_meta_selected = df_meta_not_selected[[colsample, colage, colsex]]
        meta_colheader = ["id", "age", "sex"]
        df_meta_selected.columns = meta_colheader
        df_meta_selected["sex"] = df_meta_selected["sex"].apply(lambda x: 1 if x == "M" else 2)
        self.metadata = df_meta_selected.set_index(meta_colheader[0])

WORKDIR = ".."
DIR_EXP = f"{WORKDIR}/Expression"
DIR_META = f"{WORKDIR}/Metadata"
path_exp = f"{DIR_EXP}/healthy_cohort_test_prev.txt"
# path_exp = f"{DIR_EXP}/healthy_cohort2.txt"
path_meta = f"{DIR_META}/healthy_cohort_test_prev.txt"
# path_meta = f"{DIR_META}/healthy_cohort2.txt"
data = RNAData(path_exp, path_meta)

# %%
from biolearn.data_library import DataLibrary

data = DataLibrary().get("GSE134080").load()

# %%
from biolearn.model_gallery import ModelGallery
peterclock = ModelGallery().get("TranscriptomicPredictionModel", imputation_method="none")
predicted = peterclock.predict(data)
predicted_values = predicted['Predicted']
true_values = data.metadata['age']

# %%
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import numpy as np

correlation, _ = pearsonr(true_values, predicted_values)
print(f"Pearson Correlation Coefficient: {correlation}")
plt.figure(figsize=(5, 10))
# sns.regplot(x=true_values, y=predicted_values, scatter_kws={'alpha':0.5})
corr, pval = pearsonr(true_values, predicted_values)
mae = np.mean(abs(predicted_values - true_values))
print(mae)
plt.scatter(x=true_values, y=predicted_values, alpha=0.5, s=30, color="white", linewidths=1, edgecolor="k")
plt.plot([9, 91], [9, 91], color="firebrick", linestyle="dotted", linewidth=1, label="Ground truth", zorder=5)
plt.legend(loc='center', bbox_to_anchor=(1.3, 0.5), fontsize=8).set_zorder(30)
plt.ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+5)
plt.xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+5)
plt.grid(axis="both")
plt.ylim(9, )
plt.xlim(9, 91)
plt.show()
plt.close()
# %%
