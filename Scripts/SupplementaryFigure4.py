# %%
import os
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pcarna import PCARna

# %%
def add_infection_stage_info(x):
    if str(x).endswith("V1") and x.split("-")[1][0] == "C":
        return "COVID-19 (Acute)"
    elif x.split("-")[1][0] == "R":
        return "COVID-19 (Conval.)"
    elif x.split("-")[1][0] == "L":
        return "COVID-19 (Conval.)"
    elif str(x).endswith("V2"):
        return "COVID-19 (Mid)"
    elif str(x).endswith("V3"):
        return "COVID-19 (Late)"
    else:
        pass
    
# %%
list_drop_samples = list()
list_select_samples = list()
pca = PCARna(list_drop_samples=list_drop_samples, list_select_samples=list_select_samples)
train_group = "healthy_cohort_train_prev_median_filter"
# test_group = "healthy_cohort_train_prev_median_filter_healthy_cohort_valid_prev_healthy_cohort2_viral_v1_prev_anxiety_prev_suicide_attempt_prev_depression_prev_pca_input"
test_group = "healthy_cohort_train_prev_median_filter_healthy_cohort_valid_prev_healthy_cohort2_viral_v1_prev_viral_v2_prev_viral_v3_prev_viral_recov_prev_anxiety_prev_suicide_attempt_prev_depression_prev_pca_input"
WORKDIR = ".."
DIR_EXP = f"{WORKDIR}/Expression"
DIR_META = f"{WORKDIR}/Metadata"
DIR_LASSO = f"{WORKDIR}/LASSO_INFO/{train_group}/corr_0.35_above_pval_0.05_below"    
list_meta_columns = ["Sample_Age_Group",
                    "Sample_Phenotype",
                    "Sample_Sex",
                    "Sampling_Year",
                    "Sequencing_Year",
                    "Sequencing_Platform",
                    "Read_Length"]

train_data = train_group
feat_train = os.path.join(DIR_LASSO, f"{train_group}_std_scaled.txt")
df_feat_train = pd.read_csv(feat_train, sep="\t").set_index("Project-ID")
df_feat_train = df_feat_train.drop(columns=["Sample_Age"])
exp_all = os.path.join(DIR_LASSO, f"{test_group}_std_scaled.txt")
meta_all = os.path.join(DIR_META, f"{test_group}.txt")
df_exp_all = pd.read_csv(exp_all, sep="\t", index_col=[0])
df_exp_all = df_exp_all.drop(columns=["Sample_Age"])
df_meta_all = pd.read_csv(meta_all, sep="\t")
list_new_trait = list()
for sample_id, trait in zip(df_meta_all["Project-ID"], df_meta_all["Sample_Phenotype"]):
    if trait == "ViralInfection":
        new_trait = add_infection_stage_info(sample_id)
    else:
        new_trait = trait
    list_new_trait.append(new_trait)
    
df_meta_all["Sample_Phenotype"] = list_new_trait
pca_obj, pca_exp = pca.run_pca(df_feat_train, df_exp_all, std=False)
df_pca_meta = pca.merge_pcadata_and_metadata(pca_exp, df_meta_all)

# %%
pcx_num = 1
pcy_num = 2
cm = 1/2.54
path = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 10
plot_linewidth = 2
width = 25*cm
height = 45*cm
figsize = (width, height)

fig, axes = pca.draw_pc_biplot(df_pca_meta,
                               list_meta_columns,
                               pca_obj,
                               pcx_num=1,
                               pcy_num=2, 
                               dict_legend_ncol={},
                               plot_linewidth=plot_linewidth, 
                               figsize=figsize,
                               xlim=(-150, 100),
                               ylim=(-30, 100))

fig.tight_layout(pad=0, h_pad=1.5, w_pad=1.5)
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_4.tiff", dpi=300, bbox_inches='tight')
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_4.png", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

# %%
