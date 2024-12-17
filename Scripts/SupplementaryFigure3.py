# %%
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.stats import pearsonr
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression

# %%
WORKDIR = ".."
DIR_EXP = f"{WORKDIR}/Expression"
DIR_META = f"{WORKDIR}/Metadata"
DIR_DATA = f"{WORKDIR}/Data"

list_figletters = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
list_cohorts = ["healthy_cohort_valid_prev", "healthy_cohort2", "GSE134080"]
list_clocks = ["KoreanBloodClock", "RenClockGTExAge", "PeterClock"]
dict_clocks = {"KoreanBloodClock": "Korean Blood Clock (This Study)", "PeterClock": "Peters Clock (2015)", "RenClockGTExAge": "Ren Clock (2020)"}
list_file_cohort_clock_pairs = list()
for clock in list_clocks:
    for cohort in list_cohorts:
        cohort_clock_pair = f"{cohort}_prediction_result_{clock}.txt"
        file_cohort_clock_pair = os.path.join(DIR_DATA, cohort_clock_pair)
        list_file_cohort_clock_pairs.append(file_cohort_clock_pair)

# %%
# Plot configuration
cm = 1 / 2.54
font_path = f"{WORKDIR}/Arial.ttf"
try:
    prop = fm.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = prop.get_name()
except FileNotFoundError:
    print("Font file not found. Using default font.")

plt.rcParams["font.size"] = 7
width = 3 * cm
height = 28 * cm
figsize = (width * len(list_file_cohort_clock_pairs), height)
plot_linewidth = 3

nrows = 3
ncols = int(len(list_file_cohort_clock_pairs)/nrows)

fig, axs = plt.subplots(nrows=nrows ,ncols=ncols, figsize=figsize, squeeze=False)
axs = axs.flatten()

path_dict_name_conversion_json = f"{DIR_DATA}/dictionary_name_conversion.json"
with open (path_dict_name_conversion_json, mode="rb") as fb1:
    dict_name_conversion = json.load(fb1)

for idx, (file_cohort_clock_pair, ax) in enumerate(zip(list_file_cohort_clock_pairs, axs)):
        predicted = pd.read_csv(file_cohort_clock_pair, sep="\t")
        basename_cohort_clock_pair = os.path.basename(file_cohort_clock_pair)
        cohort_name = "_".join(basename_cohort_clock_pair.split("_")[:-3])
        title_name = dict_name_conversion.get(cohort_name, cohort_name)

        true_values = predicted["ChronAge"]
        predicted_values = predicted["RNAAge"]

        corr, pval = pearsonr(true_values, predicted_values)
        mae = np.mean(abs(predicted_values - true_values))
        r2 = r2_score(true_values, predicted_values)

        # Linear regression line
        true_values_reshape = np.array(true_values).reshape(-1, 1)
        model = LinearRegression().fit(true_values_reshape, predicted_values)
        baseline_transcriptomic_age = model.predict(true_values_reshape)

        # Scatter plot and regression line
        ax.scatter(true_values, predicted_values, color="white", edgecolor="k", linewidth=0.2, label="Data points", s=20, alpha=1, zorder=2)
        ax.plot(true_values, baseline_transcriptomic_age, color="darkred", linestyle="solid", linewidth=1, label="Regression line", zorder=10)
        ax.plot([min(true_values), max(true_values)], [min(true_values), max(true_values)], color="royalblue", linestyle="dotted", linewidth=1, label="Ground truth", zorder=5)

        # Labels and annotations
        ax.set_ylabel("Predicted Age (years)", fontsize=plt.rcParams["font.size"] + 2)
        ax.set_xlabel("Chronological Age (years)", fontsize=plt.rcParams["font.size"] + 2)
        # ax.grid(axis="both", alpha=0.5)
        ax.annotate(f"Pearson's r: {corr:.2f}\nMAE: {mae:.1f}\n$R^2$: {r2:.1f}",
                    xy=(0.05, 0.73),
                    xycoords='axes fraction',
                    fontsize=plt.rcParams["font.size"] + 2,
                    zorder=3)
        ax.tick_params(axis='both', which='major', labelsize=plt.rcParams["font.size"] + 1)
        ax.set_yticks(range(20, 91, 20))
        ax.set_title(f"{title_name} (N={len(predicted)})", fontsize=plt.rcParams["font.size"] + 2, weight="bold")
        ax.set_xlim(9, 91)
        ax.set_ylim(9, 91)
        ax.legend(loc='lower right', fontsize=plt.rcParams["font.size"]).set_zorder(30)
        ax.text(-0.2, 1.1, list_figletters[idx], transform=ax.transAxes,
                fontsize=plt.rcParams["font.size"] + 4, fontweight='bold', va='top', ha='left')
        # By Yoonsung1203
        # ax.add_patch(plt.Rectangle((0.03, 0.72), 0.4, 0.24, facecolor="white", clip_on=False, linewidth = 1, edgecolor="k", transform = ax.transAxes, zorder=2, alpha=0.7))
        if idx in [1, 4, 7]:
            clock_name = basename_cohort_clock_pair.split("_")[-1].replace(".txt", "")
            text_name = dict_clocks.get(clock_name, clock_name)
            axs[idx].text(0.5, 1.15, s=text_name, horizontalalignment='center', verticalalignment='bottom', transform = axs[idx].transAxes, zorder=10, fontdict={"size":plt.rcParams["font.size"]+6, "weight": "bold"})
            width = 3.6
            height = 0.005
            # By Yoonsung1203
            axs[idx].add_patch(plt.Rectangle((0.5-width/2,1+0.13+height),width, height,facecolor='k',
                              clip_on=False,linewidth = 0, transform = axs[idx].transAxes))

plt.subplots_adjust(wspace=0.3, hspace=0.5)
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_3.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{WORKDIR}/Figures/Supplementary_Figure_3.tiff", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# %%
