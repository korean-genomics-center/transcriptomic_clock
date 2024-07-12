#%%
import os
import warnings
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from utility import Utility
warnings.filterwarnings("ignore")

# %%
def get_stratify_bins(y, nbins):
    step = (max(y) - min(y)) / nbins
    print(nbins)
    print(step)
    bins = np.arange(min(y), max(y)+step/2, step)
    bins[0] -= step/2
    bins[-1] += step/2
    print(max(y))
    print(min(y))
    print(bins)
    stratify = np.digitize(y, bins=bins, right=True)
    
    return stratify

def stratify_split_train_test_data(X: pd.DataFrame, y: pd.DataFrame, stratify=None, train_ratio=0.70, random_state=1) -> pd.DataFrame:
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=stratify, test_size= 1 - train_ratio, random_state=random_state)
    
    util = Utility()
    X_train_reidx = util.reset_indices(X_train)
    X_test_reidx = util.reset_indices(X_test)
    y_train_reidx = util.reset_indices(y_train)
    y_test_reidx = util.reset_indices(y_test)
    
    return X_train_reidx, X_test_reidx, y_train_reidx, y_test_reidx

def draw_histogram_train_test(y_train, y_test, outdir, **kwargs):
    from matplotlib.ticker import PercentFormatter
    fig, axes = plt.subplots(1,2, sharey=True)
    axes[0].hist(y_train, label="Train", color="skyblue", weights=np.ones(len(y_train)) / len(y_train), **kwargs)
    axes[1].hist(y_test, label="Test", color="orange", weights=np.ones(len(y_test)) / len(y_test), **kwargs)
    axes[0].legend(loc="best")
    axes[1].legend(loc="best")
    axes[0].set_xlim(19, 81)
    axes[1].set_xlim(19, 81)
    fig.supylabel("Proportion", fontsize=12)
    fig.supxlabel("Sample Age", fontsize=12)
    fig.suptitle("Histograms of Stratified Split Between Train and Test by Age")
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.savefig(os.path.join(outdir, "distribution_train_test_split.tiff"), dpi=600)
    plt.show()
    plt.close()

# %%
std = False
# std = True
train_ratio = 0.8
col_y = "Sample_Age"
nbins = 6
WORKDIR = "/BiO/Research/RNASeqReference/Workspace/KU10K/Results/AgingTranscriptomePaper1/healthy_illumina_DESeq2_norm_20240622"
testfilename = "testing_dataset_healthy.tsv"
train_ratio = train_ratio
path_feature = f"{WORKDIR}/FeatureTable/healthy_illumina.txt"

if std:
    dir_train_test = f"{WORKDIR}/LinearRegression/standardized/healthy_illumina"
else:
    dir_train_test = f"{WORKDIR}/LinearRegression/unstandardized/healthy_illumina"

os.makedirs(dir_train_test, exist_ok=True)
util = Utility()
X, ys = util.get_feature_table(path_feature)
X_, y = util.filter_input_data(X, ys, col_y)
stratify = get_stratify_bins(y, nbins)
X_train, X_test, y_train, y_test = stratify_split_train_test_data(X_, y, stratify=stratify, train_ratio=train_ratio, random_state=42)
draw_histogram_train_test(y_train, y_test, dir_train_test, bins=31, alpha=0.6)
if std:
    X_train_, X_test_ = util.select_dataset_target_genes(X_train, X_test, list_target_genes=None, reset_index=False)
    X_train_std, X_test_std = util.standardize_expression(X_train_, X_test_, list_target_genes=None)
    util.save_input_dataset(X_train_std, X_test_std, y_train, y_test, testfilename, dir_train_test, index=False)
else:
    util.save_input_dataset(X_train, X_test, y_train, y_test, testfilename, dir_train_test, index=False)

print(len(X_train))
print(len(X_test))
# %%
