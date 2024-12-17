# %%
import os
import subprocess

from joblib import Parallel, delayed

# %%
WORKDIR = ".."
dir_in = f"{WORKDIR}/Data"
dir_out = f"{WORKDIR}/Expression"
dir_files = f"{WORKDIR}/GeoMeans"

# %%
dict_params_1 = {
"expData":f"{dir_in}/GEO/GSE134080/GSE134080_count_preprocessed.txt",
"metaData":f"{dir_in}/GEO/GSE134080/GSE134080_metadata_preprocessed.txt",
"fileout":f"{dir_out}/GSE134080.txt",
"design":"~1",
"id_column":"Project.ID",
"gene_column":"Gene_ID",
"geoMean_column":"geoMeans",
"file_geoMeans":os.path.join(dir_files, "GSE134080_geoMeans.txt"),
"file_sizeFactors":os.path.join(dir_files, "GSE134080_sizeFactors.txt")}

dict_params_2 = {"expData":f"{dir_in}/GEO/GSE273149/GSE273149_count_preprocessed.txt",
"metaData":f"{dir_in}/GEO/GSE273149/GSE273149_metadata_preprocessed.txt",
"fileout":f"{dir_out}/GSE273149.txt",
"design":"~1",
"id_column":"Project.ID",
"gene_column":"Gene_ID",
"geoMean_column":"geoMeans",
"file_geoMeans":os.path.join(dir_files, "GSE273149_geoMeans.txt"),
"file_sizeFactors":os.path.join(dir_files, "GSE273149_sizeFactors.txt")}

dict_params_3 = {"expData":f"{dir_in}/GEO/GSE119117/GSE119117_counts_matrix_preprocessed.txt",
"metaData":f"{dir_in}/GEO/GSE119117/GSE119117_metadata_preprocessed.txt",
"fileout":f"{dir_out}/GSE119117.txt",
"design":"~1",
"id_column":"Project.ID",
"gene_column":"Gene_ID",
"geoMean_column":"geoMeans",
"file_geoMeans":os.path.join(dir_files, "GSE119117_geoMeans.txt"),
"file_sizeFactors":os.path.join(dir_files, "GSE119117_sizeFactors.txt")}

list_dict_params = [dict_params_1, dict_params_2, dict_params_3]
                    
def run(dict_params):
    cmd = f'Rscript normalize_external_expression_data.R --expData {dict_params["expData"]} --metaData {dict_params["metaData"]} --fileout {dict_params["fileout"]} --design {dict_params["design"]} --id_column {dict_params["id_column"]} --gene_column {dict_params["gene_column"]} --geoMean_column {dict_params["geoMean_column"]} --file_geoMeans {dict_params["file_geoMeans"]} --file_sizeFactors {dict_params["file_sizeFactors"]}'
    subprocess.run(cmd, shell=True)

with Parallel(n_jobs=30) as parallel:
    parallel(delayed(run)(dict_params) for dict_params in list_dict_params)


# %%
