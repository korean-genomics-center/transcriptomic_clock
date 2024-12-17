#%%
import os
import pandas as pd
from biolearn.model_gallery import ModelGallery

# %%
class RNAData():
    def __init__(self, path_exp, path_meta, colgene="Gene_ID", colsample="Project-ID", colage="Sample_Age", colsex="Sample_Sex"):
        from utility import Utility
        util = Utility()
        self.path_exp = path_exp
        self.path_meta = path_meta
        df_rna = util.read_expmtx(self.path_exp)
        # df_rna_column_drop = df_rna.drop(columns=["Description"])
        if colgene != list(df_rna.columns)[0]:
            df_rna.index.name = colgene
            df_rna = df_rna.reset_index(drop=False)
        df_rna[colgene] = df_rna[colgene].apply(lambda x: x.split(".")[0])
        df_rna_dup_row_drop = df_rna.drop_duplicates(subset=[colgene])
        self.rna = df_rna_dup_row_drop.set_index(colgene)
        df_meta_not_selected = util.read_metadata(self.path_meta, colsample=colsample)
        df_meta_selected = df_meta_not_selected[[colsample, colage, colsex]]
        meta_colheader = ["id", "age", "sex"]
        df_meta_selected.columns = meta_colheader
        df_meta_selected["sex"] = df_meta_selected["sex"].apply(lambda x: 1 if x == "M" else 2)
        self.metadata = df_meta_selected.set_index(meta_colheader[0])

# %%
WORKDIR = ".."
DIR_EXP = f"{WORKDIR}/Expression"
DIR_META = f"{WORKDIR}/Metadata"
DIR_DATA = f"{WORKDIR}/Data"
list_cohorts = ["healthy_cohort_valid_prev", "healthy_cohort2", "GSE134080", "GTEx"]

for cohort in list_cohorts:
    if cohort != "GSE134080":
        path_exp = f"{DIR_EXP}/{cohort}.txt"
        path_meta = f"{DIR_META}/{cohort}.txt"
        data = RNAData(path_exp, path_meta)
    else:
        from biolearn.data_library import DataLibrary
        data = DataLibrary().get(cohort).load()

    # Model predictions
    peterclock = ModelGallery().get("TranscriptomicPredictionModel", imputation_method="none")
    rnaage = peterclock.predict(data)
    rnaage.index.name = "id"
    chronage = data.metadata["age"]
    predicted = pd.concat([rnaage, chronage], axis=1)
    predicted.columns = ["RNAAge", "ChronAge"]
    predicted.to_csv(os.path.join(DIR_DATA, f"{cohort}_prediction_result_PeterClock.txt"), sep="\t", index=True)
