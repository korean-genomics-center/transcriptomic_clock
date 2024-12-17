# %%
from inmoose.pycombat import pycombat_seq
import pandas as pd
import numpy as np
import os

target_filename = "GSE119117"
WORKDIR = ".."
path_exp = f"{WORKDIR}/Expression/{target_filename}.txt"
path_meta = f"{WORKDIR}/Metadata/{target_filename}_pca_input.txt"

df_exp = pd.read_csv(path_exp, sep="\t", index_col=[0])
df_meta = pd.read_csv(path_meta, sep="\t", index_col=[0])
list_batch_id = df_meta["Sequencing_Batch"].to_list()
df_exp_norm = pycombat_seq(df_exp, list_batch_id)

outfilename = target_filename + "_batch_removed"
outfile = os.path.join(WORKDIR, "Expression", f"{outfilename}.txt")
df_exp_norm.to_csv(outfile, sep="\t", index=True)
# %%
