# %%
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter

# %%
WORKDIR = ".."
cm = 1/2.54
path = f"{WORKDIR}/Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 7
width = 12 * cm
height = 20 * cm
figsize = (width, height)
plot_linewidth = 1

projects = ["healthy_cohort_valid_prev", "healthy_cohort2", "viral_v1_prev", "viral_v2_prev", "viral_v3_prev", "anxiety_prev", "suicide_attempt_prev", "depression_prev"]
dict_conversion = {"healthy_cohort_valid_prev": "Healthy Cohort 1",
                    "healthy_cohort2": "Healthy Cohort 2", 
                    "anxiety_prev": "Anxiety Disorder", 
                    "depression_prev": "Major Depression", 
                    "suicide_attempt_prev": "Suicide Attempt",
                    "viral_v1_prev": "Acute Phase",
                    "viral_v2_prev": "Mid Phase",
                    "viral_v3_prev": "Late Phase",
                    "viral_recov": "Convalescent"}

# Create the figure and subplots
fig, axes = plt.subplots(4, 2, figsize=figsize)
axes = axes.flatten()

for i, proj in enumerate(projects):
    ax = axes[i]
    path_metadata = f"{WORKDIR}/FeatureTable/{proj}.txt"
    df_meta = pd.read_csv(path_metadata, sep="\t")
    list_age = list(map(int, df_meta["Sample_Age"]))
    mean_age = np.mean(list_age)
    total_samples = len(list_age)
    
    ax.hist(list_age, weights=np.ones(len(list_age)) / len(list_age), color='grey', alpha=0.5)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    
    ax.axvline(mean_age, color='tab:red', linestyle='--', linewidth=plot_linewidth)
    ax.text(mean_age+1.0, ax.get_ylim()[1] * 0.5, f'Mean Age\n={mean_age:.2f}', color='tab:red', fontsize=plt.rcParams["font.size"]+1, ha='left', weight='bold')

    title = dict_conversion[proj]
    ax.set_title(f"{title} (N={total_samples})", fontsize=plt.rcParams["font.size"] + 1, pad=0.1, weight="bold")
    ax.set_xlabel('Chronological Age (years)', fontsize=plt.rcParams["font.size"]+1)
    ax.set_ylabel('Sample Proportion (%)', fontsize=plt.rcParams["font.size"]+1)
    ax.set_xlim(15, 80)

    figletter = chr(ord('A') + i)
    ax.text(-0.35, 1.1, figletter, transform=ax.transAxes,
        fontsize=plt.rcParams["font.size"]+2, fontweight='bold', va='top', ha='right')

plt.tight_layout(rect=[0, 0, .95, .95])
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_5.tiff", dpi=300, bbox_inches="tight")
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_5.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()
# %%
