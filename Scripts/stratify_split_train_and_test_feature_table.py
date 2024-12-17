#%%
import math
import warnings
from utility import Utility
import plot
warnings.filterwarnings("ignore")

# %%
train_ratio = 0.8
col_y = "Sample_Age"
colsample = "Project-ID"
nbins = 6
WORKDIR = ".."
DIR_META = f"{WORKDIR}/Metadata"
DIR_FEAT = f"{WORKDIR}/FeatureTable"
group = "healthy_cohort1_removed_20s_selected_out_10_80_90s"
path_feature = f"{DIR_FEAT}/{group}.txt"
trainfilename = f"{group}_random_split_{math.ceil(train_ratio*100)}.txt"
testfilename = f"{group}_random_split_{math.ceil((1-train_ratio)*100)}.txt"

# %%
util = Utility()
X, ys = util.get_feature_table(path_feature)
X_, y = util.filter_input_data(X, ys, col_y)
stratify = util.get_stratify_bins(y, nbins)
X_train, X_test, y_train, y_test = util.stratify_split_train_test_data(X_, y, stratify=stratify, train_ratio=train_ratio, random_state=42)
plot.draw_histogram_train_test(y_train, y_test, DIR_META, bins=31, alpha=0.6)
util.save_input_dataset(X_train, y_train,trainfilename, DIR_FEAT, index=False)
util.save_input_dataset(X_test, y_test, testfilename, DIR_FEAT, index=False)
# %%
print(len(X_train))
print(len(X_test))
# %%
