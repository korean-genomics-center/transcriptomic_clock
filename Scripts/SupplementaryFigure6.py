# %%
import json
import os

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

from transcriptomic_clock import Clock

# %%
WORKDIR = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock"
cm = 1/2.54
path = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 7
width = 11*cm
height = 22*cm

path_dict_name_conversion_json = f"{WORKDIR}/Data/dictionary_name_conversion.json"
with open (path_dict_name_conversion_json, mode="rb") as fb1:
    dict_name_conversion = json.load(fb1)

list_test_group = list(dict_name_conversion.keys())

dir_trained_model = f"{WORKDIR}/LASSO_INFO/healthy_cohort_train_prev_median_filter/corr_0.35_above_pval_0.05_below"

# Create a figure to hold all subplots
fig, axs = plt.subplots(4, 2, figsize=(width, height), constrained_layout=True)

# Flatten the axs array to iterate over it easily
axs = axs.flatten()

# Iterate over each test group
for i, testname in enumerate(list_test_group[4:]):
    try:
        file_testing_data = f"{dir_trained_model}/{testname}_std_scaled.txt"
        os.makedirs(dir_trained_model, exist_ok=True)
        clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
        ax = clock.draw_age_prediction_scatterplot(ax=axs[i], fontsize=plt.rcParams["font.size"])
        ax.set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"]+1)
        ax.set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"]+1)
        title = dict_name_conversion[testname]
        figletter=chr(ord('A') + i)
        ax.text(-0.25, 1.1, 
                s=figletter, 
                fontdict={"fontsize":plt.rcParams["font.size"]+2, 
                          "fontweight":"bold"}, 
                ha='left', 
                va='top', 
                transform=ax.transAxes)
        ax.legend().set_visible(False)
        axs[i].set_title(f"{title}", fontsize=plt.rcParams["font.size"]+1, weight="bold", pad=0.5)
        ax.set_xlim(10, 90)
        ax.set_ylim(10, int(ax.get_yticklabels()[-1].get_text())+5)
    except Exception as e:
        print(f"Error processing {testname}: {str(e)}")

handles, labels = axs[0].get_legend_handles_labels()
axs[-1].legend(handles, labels, loc='center', fontsize=plt.rcParams["font.size"]+1)
axs[-1].axis('off')
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_6.png", dpi=300)
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_6.tiff", dpi=300)
plt.show()
plt.close()

# %%
