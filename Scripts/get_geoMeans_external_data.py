# %%
import pandas as pd
import os

# %%
WORKDIR = ".."
file_geoMeans = f"{WORKDIR}/Data/KOGIC/KOGIC_geoMeans.txt"
dir_out = f"{WORKDIR}/GeoMeans"
os.makedirs(dir_out, exist_ok=True)

# %% [GSE134080]
df_gm = pd.read_csv(file_geoMeans, sep="\t")
if len(list(df_gm.index)[0].split("_")) > 1:
    df_gm.index = list(map(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]), list(df_gm.index)))
expData =  f"{WORKDIR}/Data/GEO/GSE134080/GSE134080_count_preprocessed.txt"
df_exp = pd.read_csv(expData, sep="\t")
list_exp_genes = list(df_exp["Gene_ID"])
dict_gm = dict(zip(df_gm.index, df_gm["geoMeans"]))
list_geoMeans = list(map(lambda x: float(dict_gm.get(x, 0.0)), list_exp_genes))
dict_gene_geoMeans = dict(zip(list_exp_genes, list_geoMeans))
df_geoMeans_external = pd.DataFrame.from_dict(dict_gene_geoMeans.items())
df_geoMeans_external.columns = ["Gene_ID", "geoMeans"]
df_geoMeans_external_gene_id_indexed = df_geoMeans_external.set_index("Gene_ID")
df_geoMeans_external_gene_id_indexed.to_csv(f"{dir_out}/GSE134080_geoMeans.txt", sep="\t", index=False)

# %% [GSE273149]
# df_gm = pd.read_csv(file_geoMeans, sep="\t")
# if len(list(df_gm.index)[0].split("_")) > 1:
#     df_gm.index = list(map(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]), list(df_gm.index)))
# expData = f"{WORKDIR}/Data/GEO/GSE273149/GSE273149_count_preprocessed.txt"
# df_exp = pd.read_csv(expData, sep="\t")
# list_exp_genes = list(df_exp["Gene_ID"])
# dict_gm = dict(zip(df_gm.index, df_gm["geoMeans"]))
# list_geoMeans = list(map(lambda x: float(dict_gm.get(x, 0.0)), list_exp_genes))
# dict_gene_geoMeans = dict(zip(list_exp_genes, list_geoMeans))
# df_geoMeans_external = pd.DataFrame.from_dict(dict_gene_geoMeans.items())
# df_geoMeans_external.columns = ["Gene_ID", "geoMeans"]
# df_geoMeans_external_gene_id_indexed = df_geoMeans_external.set_index("Gene_ID")
# df_geoMeans_external_gene_id_indexed.to_csv(f"{dir_out}/GSE273149_geoMeans.txt", sep="\t", index=False)

# # %% [GSE119117] 
# df_gm = pd.read_csv(file_geoMeans, sep="\t")
# if len(list(df_gm.index)[0].split("_")) > 1:
#     df_gm.index = list(map(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]), list(df_gm.index)))
# expData = f"{WORKDIR}/Data/GEO/GSE119117/GSE119117_counts_matrix_preprocessed.txt"
# df_exp = pd.read_csv(expData, sep="\t")
# list_exp_genes = list(df_exp["Gene_ID"])
# dict_gm = dict(zip(df_gm.index, df_gm["geoMeans"]))
# list_geoMeans = list(map(lambda x: float(dict_gm.get(x, 0.0)), list_exp_genes))
# dict_gene_geoMeans = dict(zip(list_exp_genes, list_geoMeans))
# df_geoMeans_external = pd.DataFrame.from_dict(dict_gene_geoMeans.items())
# df_geoMeans_external.columns = ["Gene_ID", "geoMeans"]
# df_geoMeans_external_gene_id_indexed = df_geoMeans_external.set_index("Gene_ID")
# df_geoMeans_external_gene_id_indexed.to_csv(f"{dir_out}/GSE119117_geoMeans.txt", sep="\t", index=False)

# %%
def save_geoMeans_for_normalizing_external_data(file_geoMeans, expData, outfilename):
    df_gm = pd.read_csv(file_geoMeans, sep="\t")
    if len(list(df_gm.index)[0].split("_")) > 1:
        df_gm.index = list(map(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]), list(df_gm.index)))
    df_exp = pd.read_csv(expData, sep="\t")
    list_exp_genes = list(df_exp["Gene_ID"])
    dict_gm = dict(zip(df_gm.index, df_gm["geoMeans"]))
    list_geoMeans = list(map(lambda x: float(dict_gm.get(x, 0.0)), list_exp_genes))
    dict_gene_geoMeans = dict(zip(list_exp_genes, list_geoMeans))
    df_geoMeans_external = pd.DataFrame.from_dict(dict_gene_geoMeans.items())
    df_geoMeans_external.columns = ["Gene_ID", "geoMeans"]
    df_geoMeans_external_gene_id_indexed = df_geoMeans_external.set_index("Gene_ID")
    df_geoMeans_external_gene_id_indexed.to_csv(f"{dir_out}/{outfilename}_geoMeans.txt", sep="\t", index=False)
