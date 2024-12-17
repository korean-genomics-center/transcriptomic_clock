# %%
import os
import json
import pandas as pd
from racpy import RNAAgeCalc
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# %%
def get_expdata(path_exp, list_gene_excl, colgene="Gene_ID", tagensembl="ENS"):
    expdata = pd.read_csv(path_exp, sep=",")
    first_col = list(expdata.columns)[0]
    expdata = expdata.rename(columns={first_col: colgene})
    expdata = expdata[~expdata[colgene].str.startswith(tagensembl)]
    expdata = expdata[~expdata[colgene].isin(list_gene_excl)]
    expdata_drop_dup = expdata.drop_duplicates(subset=[colgene], keep="first")
    expdata_set_idx = expdata_drop_dup.set_index(colgene)
    
    return expdata_set_idx

def get_metadata(path_meta, colsampleid="Project-ID", colage="Sample_Age"):
    metadata = pd.read_csv(path_meta, sep="\t")
    chronage = metadata[[colsampleid, colage]]

    return chronage

# %%
WORKDIR = ".."
dir_exp = f"{WORKDIR}/Data"
dir_meta = f"{WORKDIR}/Metadata"
dir_excl_gene = f"{WORKDIR}/Data"
file_excl_gene = os.path.join(dir_excl_gene, "gene_list_excl_RNAAgeClock.txt")
list_cohorts = ["healthy_cohort_valid_prev", "healthy_cohort2", "GSE134080"]

path_dict_name_conversion_json = f"{WORKDIR}/Data/dictionary_name_conversion.json"
with open (path_dict_name_conversion_json, mode="rb") as fb1:
    dict_name_conversion = json.load(fb1)

list_title = list(map(lambda x: dict_name_conversion.get(x, x), list_cohorts))

# %%
# Plot configuration
cm = 1 / 2.54
font_path = f"{WORKDIR}/Arial.ttf"  # Check this path or provide a fallback
try:
    prop = fm.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = prop.get_name()
except FileNotFoundError:
    print("Font file not found. Using default font.")

plt.rcParams["font.size"] = 7
width = 10 * cm
height = 8 * cm
figsize = (width * len(list_cohorts), height)  # Adjust figsize for the number of cohorts
plot_linewidth = 3

# Initialize plots
fig, axs = plt.subplots(nrows=1, ncols=len(list_cohorts), figsize=figsize, squeeze=False)
axs = axs.flatten()

for cohort in list_cohorts:
    with open(file_excl_gene, mode="r") as fr:
        list_gene_excl = list(map(lambda x: x.rstrip(), fr.readlines()))

    path_exp = os.path.join(dir_exp, f"{cohort}_rawcount.csv")
    exprdata = get_expdata(path_exp, list_gene_excl, colgene="Gene_ID", tagensembl="ENS")

    path_meta = os.path.join(dir_meta, f"{cohort}.txt")
    chronage = get_metadata(path_meta, colsampleid="Project-ID", colage="Sample_Age")
    
    if not os.path.exists(os.path.join(dir_exp, f"{cohort}_prediction_result_RenClock_GTExAge.txt")):
        renclock = RNAAgeCalc(tissue = "blood", exprtype = "count", idtype = "symbol", stype = "Caucasian", signature = "GTExAge")
        predicted = renclock.predict_age(exprdata=exprdata, chronage=chronage)
        predicted.to_csv(os.path.join(dir_exp, f"{cohort}_prediction_result_RenClock_GTExAge.txt"), sep="\t", index=True)

    else:
        predicted = pd.read_csv(os.path.join(dir_exp, f"{cohort}_prediction_result_RenClock_GTExAge.txt"), sep="\t")
