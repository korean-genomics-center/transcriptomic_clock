# %%
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import os
import gzip
import pickle

# %%
class Utility():
    def __init__(self):
        self.colsample="SUBJID"
        self.gene_tag = "ENSG"
        self.remove_genes = ["PAR_Y", "CTC-338M12.4"]
        self.verbose = False
        self.returnX = False
        self.returnY = True
        self.sep = "\t"

    def dump_pickle(self, data: str, dir_out: str, filename: str) ->  str:
        outfilename = os.path.join(dir_out, filename)
        os.makedirs(os.path.dirname(outfilename), exist_ok=True)
        with gzip.open(outfilename, 'wb') as f:
            pickle.dump(data, f)

        return None

    def load_pickle(self, dir_out: str, filename: str) -> str:
        loadfilename = os.path.join(dir_out, filename)
        with gzip.open(loadfilename,'rb') as f:
            data = pickle.load(f)
        
        return data
    
    def reset_indices(self, df: pd.DataFrame) -> pd.DataFrame:
        df_reidx = df.reset_index(drop=True)

        return df_reidx

    def get_end_idx_gene_list(self, df_feature: pd.DataFrame, gene_tag: str) -> int:
        list_header = list(df_feature.columns)
        list_gene_idx = list()
        for idx, column in enumerate(list_header):
            if str(column).startswith(gene_tag):
                list_gene_idx.append(idx)
        end_idx_gene = int(max(list_gene_idx))

        return end_idx_gene

    def lookup_available_variables(self, df_dataset: pd.DataFrame) -> list:
        list_dataset = df_dataset.columns.to_list()
        
        return list_dataset

    def overview_variables(self, X: pd.DataFrame, ys: pd.DataFrame) -> list:
        list_x = self.lookup_available_variables(X)
        print(f"Xs: {list_x}")

        list_y = self.lookup_available_variables(ys)
        print(f"Ys: {list_y}")

        if self.returnX:
            return list_x

        if self.returnY:
            return list_y
        
        if (self.returnX and self.returnY):
            return list_x, list_y

    def get_list_selected_genes(self, file_linear_reg):
        list_selected_genes = list()
        with open(file_linear_reg, mode="r") as fr:
            _skiprow = fr.readline()
            for line in fr:
                record = line.rstrip("\n").split("\t")
                gene = record[0]
                list_selected_genes.append(gene)

        return list_selected_genes

    def save_selected_genes(self, dir_out, list_selected_genes):
        fileingenename = os.path.join(dir_out, "features_ingenename.txt")
        with open(fileingenename, mode="w") as fw:
            for gene in list_selected_genes:
                fw.write(gene + "\n")
    
    def get_feature_table(self, feature_table) -> pd.DataFrame:
        df_feat = pd.read_csv(feature_table, sep=self.sep, low_memory=False)
        end_idx_gene = self.get_end_idx_gene_list(df_feat, self.gene_tag)
        X = df_feat.iloc[:, :(end_idx_gene+1)]
        for remove_tag in self.remove_genes:
            list_filt_colnames = list(filter(lambda x: remove_tag not in x, list(X.columns[1:])))
        X_filtered = X[[X.columns[0]] + list_filt_colnames]
        list_new_colnames = [X_filtered.columns[0]] + list(map(lambda x: x.split(".")[0] + "_" + "_".join(x.split("_")[1:]) if len(x.split("."))>1 else x, list_filt_colnames))
        X_filtered.columns = list_new_colnames
        ys = df_feat.iloc[:, (end_idx_gene+1):]

        if self.verbose:
            self.overview_variables(X, ys)

        return X_filtered, ys

    def filter_input_data(self, X, ys, col_y):
        list_idx_null = list(np.where(~pd.isnull(ys[col_y]))[0])
        y = ys[col_y].iloc[list_idx_null]
        X_ = X.iloc[list_idx_null, :]
        print(f"'{col_y}' selected")

        return X_, y

    def get_list_gene_name(self, X):
        list_gene_name = list(X.iloc[:, 1:].columns)

        return list_gene_name
    
    def get_target_genes(self, path_target):
        with open(path_target, mode="r") as fr:
            list_target_genes = fr.readlines()
        
        list_target_genes = list(map(lambda x: x.rstrip("\n"), list_target_genes))

        return list_target_genes

    def select_dataset_target_genes(self, X_train_, X_test_, list_target_genes, reset_index=False):
        if list_target_genes == []:
            list_target_genes = None
        
        X_train_set_idx = X_train_.set_index(self.colsample)
        if list_target_genes is not None:
            X_train_set_idx = X_train_set_idx[list_target_genes]

        X_test_set_idx = X_test_.set_index(self.colsample)
        if list_target_genes is not None:
            X_test_set_idx = X_test_set_idx[list_target_genes]
      
        if reset_index:
            X_train_reset_idx = X_train_set_idx.reset_index(drop=False).rename(columns={"index": self.colsample})    
            X_test_reset_idx = X_test_set_idx.reset_index(drop=False).rename(columns={"index": self.colsample})

            return X_train_reset_idx, X_test_reset_idx

        return X_train_set_idx, X_test_set_idx

    def standardize_expression(self, X_train_set_idx, X_test_set_idx,list_target_genes):        
        scaler = StandardScaler().fit(X_train_set_idx)

        list_train_samples = list(X_train_set_idx.index)
        list_train_genes = list(X_train_set_idx.columns)

        list_test_samples = list(X_test_set_idx.index)
        list_test_genes = list(X_test_set_idx.columns)

        if list_target_genes == []:
            list_target_genes = None

        if list_target_genes is not None:
            X_train_set_idx_std = pd.DataFrame(scaler.transform(X_train_set_idx), index=list_train_samples, columns=list_target_genes)
        else:
            X_train_set_idx_std = pd.DataFrame(scaler.transform(X_train_set_idx), index=list_train_samples, columns=list_train_genes)
        

        if list_target_genes is not None:
            X_test_set_idx_std = pd.DataFrame(scaler.transform(X_test_set_idx), index=list_test_samples, columns=list_target_genes)
        else:
            X_test_set_idx_std = pd.DataFrame(scaler.transform(X_test_set_idx), index=list_test_samples, columns=list_test_genes)
    
        X_train_reset_idx_std = X_train_set_idx_std.reset_index(drop=False).rename(columns={"index": self.colsample})    
        X_test_reset_idx_std = X_test_set_idx_std.reset_index(drop=False).rename(columns={"index": self.colsample})
    
        return X_train_reset_idx_std, X_test_reset_idx_std

    def save_input_dataset(self, X_train, X_test, y_train, y_test, testfilename, outdir, **kwargs):
        training_dataset = pd.concat([X_train, y_train], axis=1)
        training_dataset.to_csv(os.path.join(outdir,"training_dataset.tsv"), sep=self.sep, **kwargs)

        testing_dataset = pd.concat([X_test, y_test], axis=1)
        testing_dataset.to_csv(os.path.join(outdir, testfilename), sep=self.sep, **kwargs)
# %%
