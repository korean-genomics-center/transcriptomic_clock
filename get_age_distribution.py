# %%
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter

# Font and figure properties
cm = 1/2.54
path = "./Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
width = 10 * cm
height = 20 * cm
figsize = (width, height)
plot_linewidth = 1

projects = ["healthy_illumina", "healthy_bgi", "viral_v1", "viral_v2", "viral_v3", "anxiety", "suicide_attempt", "depression"]
dict_conversion = {"healthy_illumina": "Healthy Cohort 1",
                    "healthy_bgi": "Healthy Cohort 2", 
                    "anxiety": "Anxiety Disorder", 
                    "depression": "Major Depression", 
                    "suicide_attempt": "Suicide Attempt",
                    "viral_v1": "Acute Phase",
                    "viral_v2": "Mid Phase",
                    "viral_v3": "Late Phase",
                    "viral_recov": "Convalescent"}

# Create the figure and subplots
fig, axes = plt.subplots(4, 2, figsize=figsize)
axes = axes.flatten()
letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

for i, proj in enumerate(projects):
    ax = axes[i]
    path_metadata = f"./FeatureTable/{proj}.txt"
    df_meta = pd.read_csv(path_metadata, sep="\t")
    list_age = list(map(int, df_meta["Sample_Age"]))
    mean_age = np.mean(list_age)
    total_samples = len(list_age)
    
    # Plot histogram
    ax.hist(list_age, weights=np.ones(len(list_age)) / len(list_age), color='grey', alpha=0.5)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    
    # Plot mean age as vertical line
    ax.axvline(mean_age, color='tab:red', linestyle='--', linewidth=plot_linewidth)
    ax.text(mean_age+1.0, ax.get_ylim()[1] * 0.5, f'Mean Age\n={mean_age:.2f}', color='tab:red', fontsize=plt.rcParams["font.size"]+1, ha='left', weight='bold')

    # Set title and labels
    title = dict_conversion[proj]
    ax.set_title(f"{title} (N={total_samples})", fontsize=plt.rcParams["font.size"] + 1, pad=0.1, weight="bold")
    ax.set_xlabel('Chronological Age (years)', fontsize=plt.rcParams["font.size"]+1)
    ax.set_ylabel('Sample Proportion (%)', fontsize=plt.rcParams["font.size"]+1)
    ax.set_xlim(15, 80)
    # Add figure letter
    fig.text(0.02 + (i % 2) * 0.48, 0.95 - (i // 2) * 0.24, letters[i], ha='center', va='center', fontsize=plt.rcParams["font.size"] + 2, weight='bold')

plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to make space for the letters
plt.savefig("Figures/Extended_Data_Fig_5.pdf", dpi=300, bbox_inches="tight")
plt.savefig("Figures/Extended_Data_Fig_5.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()
# %%
