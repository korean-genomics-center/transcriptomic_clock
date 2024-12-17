# %%
import gzip
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
# from sklearn.ensemble import RandomForestRegressor
# https://scikit-learn.org/0.17/auto_examples/linear_model/plot_lasso_model_selection.html
from sklearn.linear_model import Lasso, LassoCV, LassoLarsIC, LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
# from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedKFold


# %%
class Clock():
    def __init__(self, dir_trained_model, file_testing_data=None, col_sampleid="Project-ID", col_y="Sample_Age"):
        self.dir_trained_model = dir_trained_model
        os.makedirs(self.dir_trained_model, exist_ok=True)
        self.file_testing_data = file_testing_data
        self.col_sampleid = col_sampleid
        self.col_y = col_y

    def find_no_match_gene(self, list_target_samples_gene_name, list_query_gene_name):
        list_diff = list(set(list_query_gene_name).difference(set(list_target_samples_gene_name)))
        diff_genes = ",".join(list_diff)
        num_diff_genes = len(list_diff)
        if num_diff_genes == 0:
            print("all genes are found in the original gene list")
        else:
            print(f"{num_diff_genes} genes are not found in the original gene list")

        assert num_diff_genes == 0, print(f"non-matching genes: {diff_genes}")

    def get_array_array_gene_exp(self, X, list_selected_gene_name):
        array_array_gene_exp = np.array(X[list_selected_gene_name].values).astype(np.float64)

        return array_array_gene_exp

    def get_array_array_variable(self,y):
        array_array_variable = np.array(np.ravel(y.values)).astype(np.float64)

        return array_array_variable

    def tune_hyperparameters_lasso_gradient_decscent(self, X, y, n_splits=10, n_repeats=3, random_state=42, tol=1e-2, num_threads=50, max_alphas=500):
        alphas = np.arange(0, int(max_alphas)+1, 10)
            
        cv = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
        model = LassoCV(alphas=alphas, cv=cv, tol=tol, n_jobs=num_threads, verbose=True).fit(X, y)

        return model

    def get_dictionary_hyperparameters_lasso_gradient_decscent(self, model: object, dir_lasso: str) -> dict:
        opt_alpha = model.alpha_
        dict_hyperparam = {"alpha":opt_alpha}

        df_hyperparam = pd.DataFrame.from_dict([dict_hyperparam])
        outfilehyperparam = os.path.join(dir_lasso, "hyperparameters.txt")
        os.makedirs(os.path.dirname(outfilehyperparam), exist_ok=True)
        df_hyperparam.to_csv(outfilehyperparam, sep="\t", index=True)

        return dict_hyperparam

    def plot_ic_criterion(self, model, name, color):
        alpha_ = model.alpha_
        alphas_ = model.alphas_
        criterion_ = model.criterion_
        plt.plot(-np.log10(alphas_), criterion_, '--', color=color,
                linewidth=3, label='%s criterion' % name)
        plt.axvline(-np.log10(alpha_), color=color, linewidth=3,
                    label='alpha: %s estimate' % name)
        plt.xlabel('minuslog10(alpha)', fontdict={"fontsize": 16})
        plt.ylabel('criterion', fontdict={"fontsize": 16})

    def tune_hyperparameters_lasso_aic_bic(self, X, y, dir_lasso, verbose=True):
        model_bic = LassoLarsIC(criterion='bic')
        model_bic.fit(X, y)
        alpha_bic_ = model_bic.alpha_

        model_aic = LassoLarsIC(criterion='aic')
        model_aic.fit(X, y)
        alpha_aic_ = model_aic.alpha_
        
        if verbose:
            outstatpath = os.path.join(dir_lasso, "hyperparam_aic_bic_stat.txt")
            outfigpathtiff = os.path.join(dir_lasso, "hyperparam_aic_bic_figure.tiff")
            outfigpathpng = os.path.join(dir_lasso, "hyperparam_aic_bic_figure.png")
            dict_stat = {"alphas": model_aic.alphas_,
                        "AIC criterion": model_aic.criterion_,
                        "BIC criterion": model_bic.criterion_}
            df_stat = pd.DataFrame.from_dict(dict_stat, orient="columns")
            df_stat.to_csv(outstatpath, sep="\t")
            plt.figure(figsize=(10, 10))
            self.plot_ic_criterion(model_aic, "AIC", "b")
            self.plot_ic_criterion(model_bic, "BIC", "r")
            plt.legend()
            plt.savefig(outfigpathtiff, dpi=600)
            plt.savefig(outfigpathpng, dpi=600)
            plt.show()
            plt.close()
            
        return alpha_aic_, alpha_bic_
        
    def get_dictionary_hyperparameters_lasso_aic_bic(self, alpha, criterion, dir_lasso) -> dict:
        if criterion == "AIC":
            dict_hyperparam = {"alpha":alpha}
        elif criterion == "BIC":
            dict_hyperparam = {"alpha":alpha}
        else:
            print("unknown criterion")

        df_hyperparam = pd.DataFrame.from_dict([dict_hyperparam])
        outfilehyperparam = os.path.join(dir_lasso, "hyperparameters.txt")
        os.makedirs(os.path.dirname(outfilehyperparam), exist_ok=True)
        df_hyperparam.to_csv(outfilehyperparam, sep="\t", index=True)

        return dict_hyperparam

    def get_testing_dataset(self):
        df_testing_dataset = pd.read_csv(self.file_testing_data, sep="\t")

        return df_testing_dataset

    def get_X(self):
        df_testing_dataset = self.get_testing_dataset()
        list_col_header = list(df_testing_dataset.columns)
        list_col_header.remove(self.col_sampleid)
        list_col_header.remove(self.col_y)
        X = df_testing_dataset[list_col_header].to_numpy()

        return X

    def get_y(self):
        df_testing_dataset = self.get_testing_dataset()
        y = df_testing_dataset[self.col_y].to_numpy()
        
        return y

    def get_trained_model(self):
        path_trained_model = os.path.join(self.dir_trained_model, f"model_train_{self.col_y}.pk.gz")
        with gzip.open(path_trained_model, mode="rb") as fr:
            model_train = pickle.load(fr)

        return model_train

    def get_ingenenames(self):
        file_ingenename = os.path.join(self.dir_trained_model, "features_ingenename.txt")
        list_ingenenames = list()
        with open(file_ingenename, mode="r") as fr:
            for line in fr:
                record = line.rstrip("\n").split("\t")
                ingenename = record[0]
                list_ingenenames.append(ingenename)
        
        return list_ingenenames

    def get_predicted_y(self, impute="zero"):
        model_train = self.get_trained_model()
        X = self.get_X()
        if impute == "zero":
            X = np.nan_to_num(X, nan=0)
        yhat = model_train.predict(X)

        return yhat

    def fit_model_lasso(self, X_train, y_train, dict_hyperparam):
        model = Lasso(alpha=dict_hyperparam["alpha"]).fit(X_train, y_train)
            
        return model

    def get_prediction_error(self):
        yhat = self.get_predicted_y()
        y = self.get_y()
        error = yhat - y

        return error
    
    def get_abs_error(self):
        error = self.get_prediction_error()
        abs_error = abs(error)

        return abs_error

    def calculate_mae(self):
        yhat = self.get_predicted_y()
        y = self.get_y()
        mae = np.mean(abs(yhat - y))

        return mae

    def calculate_rmse(self):
        yhat = self.get_predicted_y()
        y = self.get_y()
        rmse = np.sqrt(mean_squared_error(y, yhat))
        
        return rmse

    def calculate_r2(self):
        yhat = self.get_predicted_y()
        y = self.get_y()
        r2 = r2_score(y, yhat)

        return r2

    def calculate_pearsonr(self):
        yhat = self.get_predicted_y()
        y = self.get_y()
        corr, pval = pearsonr(y, yhat)

        return corr, pval

    def evalute_prediction_accuracy(self):
        mae = self.calculate_mae()
        rmse = self.calculate_rmse()
        r2 = self.calculate_r2()
        corr, pval = self.calculate_pearsonr()
        
        dict_eval = {"mae": mae, "rmse": rmse, "r2": r2, "corr": corr, "pval": pval}

        return dict_eval
    
    def make_dataframe_error_statistics_per_sample(self):
        df_testing_dataset = self.get_testing_dataset()
        sampleid = df_testing_dataset[self.col_sampleid]
        y = self.get_y()
        yhat = self.get_predicted_y()
        error = self.get_prediction_error()
        abs_error = self.get_abs_error()
        
        dict_stats = {self.col_sampleid: sampleid, "y": y, "yhat": yhat, "error": error, "abs_error": abs_error}
        df_stats = pd.DataFrame(dict_stats)
        df_stats_sorted = df_stats.sort_values(by="abs_error", ascending=False)
        
        return df_stats_sorted

    def draw_age_prediction_scatterplot(self, ax, fontsize=7):
        real_age = self.get_y()
        real_age_reshape = np.array(real_age).reshape(-1, 1)
        observed_transcriptomic_age = self.get_predicted_y()
        model = LinearRegression().fit(real_age_reshape, observed_transcriptomic_age)
        baseline_transcriptomic_age = model.predict(real_age_reshape)
        ax.scatter(real_age, observed_transcriptomic_age, color="white", edgecolor="k", linewidth=0.2, label="Data points", s=5, alpha=1, zorder=2)
        ax.plot(real_age, baseline_transcriptomic_age, color="darkred", linestyle="solid", linewidth=1, label="Regression line", zorder=10)
        ax.plot([min(real_age), max(real_age)], [min(real_age), max(real_age)], color="royalblue", linestyle="dotted", linewidth=1, label="Ground truth", zorder=5)
        dict_eval = self.evalute_prediction_accuracy()
        round_mae = round(float(dict_eval["mae"]), 2)
        round_corr = round(float(dict_eval["corr"]), 2)
        round_r2 = round(float(dict_eval["r2"]), 2)
        ax.annotate(f"Pearson's r:{round_corr}\nMAE:{round_mae}\n$R^2$:{round_r2}", 
                    xy=(0.05, 0.75), 
                    xycoords='axes fraction', 
                    fontsize=fontsize)
        ax.tick_params(axis='both', which='major', labelsize=fontsize)
        ax.set_xlim(left=15, right=75)
        ax.set_ylim(bottom=15)
        ax.legend(loc='lower right', fontsize=fontsize).set_zorder(30)

        return ax

    def get_indices_nonzero_efficient(self, model_train):
        ind_nonzero_coef = list(filter(lambda x : model_train.coef_[x] != 0, range(len(model_train.coef_))))

        return ind_nonzero_coef

    def get_dict_selected_features(self):
        model_train = self.get_trained_model()
        list_gene_name = self.get_ingenenames()
        ind_nonzero_coef = self.get_indices_nonzero_efficient(model_train)
        list_selected_features = list(map(lambda x : list_gene_name[x], ind_nonzero_coef))
        list_nonzero_coef = list(map(lambda x: list(model_train.coef_)[x], ind_nonzero_coef))
        dict_selected_features = {"feat":list_selected_features, "coef":list_nonzero_coef}

        return dict_selected_features

    def make_dataframe_clock_features(self):
        dict_selected_features = self.get_dict_selected_features()
        df_features = pd.DataFrame.from_dict(dict_selected_features)
        df_features.columns = ["Selected_Features", "Regression_Coefficient"]

        df_features["GeneSymbol"] = df_features["Selected_Features"].apply(lambda x: "_".join(x.split("_")[1:]))
        df_features["Direction"] = df_features["Regression_Coefficient"].apply(lambda x: "Positive" if float(x) > 0 else "Negative")
        df_features["Abs_Regression_Coefficient"] = df_features["Regression_Coefficient"].apply(lambda x: abs(float(x)))
        df_features_sorted = df_features.sort_values(by=["Abs_Regression_Coefficient"], ascending=False)
        df_features_sorted = df_features_sorted.drop(columns=["Abs_Regression_Coefficient"])
        outfeaturename = os.path.join(self.dir_trained_model, "regression_results.txt")
        df_features_sorted.to_csv(outfeaturename, sep="\t", index=False)

        return df_features_sorted

    def draw_feature_importance_barplot(self, ax, fontsize=18):
        df_features_sorted = self.make_dataframe_clock_features()
        sns.barplot(data=df_features_sorted, 
                    x="Regression_Coefficient", 
                    y="GeneSymbol", 
                    orient="h", 
                    hue="Direction", 
                    palette="RdBu", 
                    ax = ax,
                    legend=True)

        # Formatting ticks and labels
        ax.set_yticks(range(df_features_sorted.shape[0]), df_features_sorted["GeneSymbol"], fontstyle="italic", fontsize=fontsize)

        return ax

    def get_annotation_prediction_evaluation_plot(self):
        list_gene_name = self.get_ingenenames()
        dict_selected_features = self.get_dict_selected_features()
        list_selected_features = dict_selected_features["feat"]
        mae = self.calculate_mae()
        corr = self.calculate_pearsonr()

        return {"NumFeatures": len(list_selected_features), "NumGenes": len(list_gene_name), "MAE":mae, "corr":corr}

    def summarize_performance_statistics(self, df_annot_age_eval_plot_all):
        col = list(df_annot_age_eval_plot_all.columns)
        mae = df_annot_age_eval_plot_all.loc["MAE"]
        corr = df_annot_age_eval_plot_all.loc["corr"]
        numgene = df_annot_age_eval_plot_all.loc["NumGenes"]
        list_num_input_genes = list(map(int, numgene.to_list()))
        numfeat = df_annot_age_eval_plot_all.loc["NumFeatures"]
        list_num_features = list(map(int, numfeat.to_list()))

        return col, mae, corr, list_num_input_genes, list_num_features

    def draw_plot_performance_profile(self,
                                      col_train, 
                                      col_test,  
                                      mae_train, 
                                      mae_test, 
                                      corr_train, 
                                      corr_test, 
                                      list_num_features_train, 
                                      list_num_input_genes_train, 
                                      plot_linewidth=3.5,
                                      figsize=(30, 10)):

        fig, ax1 = plt.subplots(figsize=figsize)

        # Extract the correlation values from the tuples
        corr_train_values = [x[0] for x in corr_train]
        corr_test_values = [x[0] for x in corr_test]

        # Ensure col_test is an iterable with correct tick positions
        ax1.set_xticks(col_test)

        # Modify col_train for proper labeling
        col_train_modif = list(map(lambda x: str(x).split('.')[-1], col_train))
        col_train_modif = list(map(lambda x: x + "0" if len(x) < 2 else x, col_train_modif))

        # Ensure the number of labels matches the number of ticks
        xtick_labels = list(map(lambda x: f".{x[0]}\n\n{x[1]}\n\n{x[2]}", zip(col_train_modif, list_num_features_train, list_num_input_genes_train)))

        # Adjust the number of ticks and labels to be equal
        if len(col_test) > len(xtick_labels):
            col_test = col_test[:len(xtick_labels)]
        elif len(xtick_labels) > len(col_test):
            xtick_labels = xtick_labels[:len(col_test)]

        # Ensure xtick_labels are strings
        xtick_labels = [str(label) for label in xtick_labels]
        color = 'darkred'
        ax1.set_xticks(col_test)
        ax1.set_xticklabels(xtick_labels, fontsize=plt.rcParams["font.size"])
        ax1.plot(col_train, mae_train, color=color, linestyle="dashed", linewidth=plot_linewidth)
        ax1.plot(col_test, mae_test, color=color, linestyle="solid", linewidth=plot_linewidth)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.set_yticks(np.linspace(0, 20, 5))
        ax1.set_yticklabels(np.linspace(0, 20, 5), fontsize=plt.rcParams["font.size"])
        ax1.set_ylim(4.5, 20)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'darkblue'
        ax2.plot(col_train, corr_train_values, color=color, linestyle="dashed", linewidth=plot_linewidth)
        ax2.plot(col_test, corr_test_values, color=color, linestyle="solid", linewidth=plot_linewidth)
        ax2.set_yticks(np.linspace(0, 1, 5))
        ax2.set_yticklabels(np.linspace(0, 1, 5), fontsize=plt.rcParams["font.size"])
        ax2.set_ylim(0, 1)
        ax2.tick_params(axis='y', labelcolor=color)

        return fig, ax1, ax2
# %%
