#%%
import os
import warnings
from utility import Utility
import plot
warnings.filterwarnings("ignore")

# %%
util = Utility()
WORKDIR = ".."
group = "healthy_cohort1"
DIR_EXP = os.path.join(WORKDIR, "Expression")
DIR_META = os.path.join(WORKDIR, "Metadata")
colsample = "Project-ID"
colgene = "Gene_ID"

# %%
path_metadata = os.path.join(DIR_META, f"{group}.txt")
df_meta = util.read_metadata(path_metadata, colsample=colsample)
df_meta_added_age_group = util.add_df_sample_age_group(df_meta)
list_age_group = util.get_list_column_val(df_meta_added_age_group, column="Sample_Age_Group")
list_age = util.get_list_column_val(df_meta_added_age_group, column="Sample_Age")
plot.draw_barplot(group, list_age_group, DIR_META, group)

# %%
num_samples_selected = 55
list_age_group_selected = [20]
file_suffix = util.get_file_suffix_age_group(list_age_group_selected)
list_samples = util.get_list_column_val(df_meta_added_age_group, column=colsample)
dict_samples_age_group = util.get_dictionary_sample_age_group(list_samples, list_age_group)
selected_samples = util.get_randomly_selected_samples_certain_age_group(dict_samples_age_group, list_age_group_selected, num_samples_selected=num_samples_selected)
with open(os.path.join(DIR_META, f"list_samples_removed_{file_suffix}.txt"), mode="w") as fw:
    for sample in selected_samples:
        fw.write(sample + "\n")
df_meta_filtered_in = util.filter_out_df_excl_samples(df_meta_added_age_group, selected_samples, colsample)
list_age_group_modified = util.get_list_column_val(df_meta_filtered_in, column="Sample_Age_Group")
newfilename = f"{group}_removed_{file_suffix}s"
plot.draw_barplot(group, list_age_group_modified, DIR_META, newfilename)

# %%
list_sample_selected = util.get_list_column_val(df_meta_filtered_in, column=colsample)

path_metadata = os.path.join(DIR_META, f"{group}.txt")
df_meta_filt = util.read_metadata(path_metadata, colsample=colsample)
df_meta_selected = util.select_samples_meta(df_meta, list_sample_selected, colsample)
path_metadata_filt = os.path.join(DIR_META, f"{newfilename}.txt")
util.save_metadata(df_meta_selected, path_metadata_filt)

path_exp_mtx = os.path.join(DIR_EXP, f"{group}.txt")
df_exp = util.read_expmtx(path_exp_mtx)
df_exp_selected = util.select_samples_df(df_exp, list_sample_selected, colgene)
path_exp_filt = os.path.join(DIR_EXP, f"{newfilename}.txt")
util.save_expmtx(df_exp_selected, path_exp_filt)

# %%
list_age_group_removed = [10, 80, 90]
file_suffix = util.get_file_suffix_age_group(list_age_group_removed)
df_meta_filtered_out = util.remove_df_certain_age_group(df_meta_filtered_in, list_age_group_removed)
list_sample_filtered_out = util.get_sampleid_removed_age_group(df_meta_filtered_out, colsample)
with open(os.path.join(DIR_META, f"list_samples_removed_{file_suffix}.txt"), mode="w") as fw:
    for sample in list_sample_filtered_out:
        fw.write(sample + "\n")
df_meta_dbl_filtered_in = util.filter_out_df_excl_samples(df_meta_filtered_in, list_sample_filtered_out, colsample)
list_age_group = util.get_list_column_val(df_meta_dbl_filtered_in, column="Sample_Age_Group")
list_age = util.get_list_column_val(df_meta_dbl_filtered_in, column="Sample_Age")
newfilename += f"_selected_out_{file_suffix}s"
plot.draw_barplot(group, list_age_group, DIR_META, newfilename)

# %%
list_sample_selected = util.get_list_column_val(df_meta_dbl_filtered_in, column=colsample)

path_metadata = os.path.join(DIR_META, f"{group}.txt")
df_meta_filt = util.read_metadata(path_metadata, colsample=colsample)
df_meta_selected = util.select_samples_meta(df_meta, list_sample_selected, colsample)
path_metadata_filt = os.path.join(DIR_META, f"{newfilename}.txt")
util.save_metadata(df_meta_selected, path_metadata_filt)

path_exp_mtx = os.path.join(DIR_EXP, f"{group}.txt")
df_exp = util.read_expmtx(path_exp_mtx)
df_exp_selected = util.select_samples_df(df_exp, list_sample_selected, colgene)
path_exp_filt = os.path.join(DIR_EXP, f"{newfilename}.txt")
util.save_expmtx(df_exp_selected, path_exp_filt)
