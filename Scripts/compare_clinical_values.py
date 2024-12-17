# %%
import numpy as np
import seaborn as sns
from matplotlib.cm import get_cmap
import statsmodels.api as sm
from scipy.stats import pearsonr

def group_patients_by_crp_level(dict_crp):
    dict_crp_group = dict()
    for sample, crp in dict_crp.items():
        if crp <= 1:
            group = "Low"
        elif crp > 1:
            group = "High"
        else:
            group = "Unknown"

        dict_crp_group[sample] = group

    return dict_crp_group        

def get_column_crp_group(df_crf, dict_crp_group, colsample):
    list_crp_group = list(map(lambda x: dict_crp_group.get(x, None), df_crf[colsample]))

    return list_crp_group

def draw_regplot(data, x, y, hue=None, palette=None, **kwargs):
    regplots = list()
    if hue is not None:
        levels = data[hue].unique()
        default_colors = get_cmap('tab10')
        palette = {k: default_colors(i) for i, k in enumerate(levels)}
        for key in levels:
            regplots.append(
                sns.regplot(
                    x=x,
                    y=y,
                    data=data[data[hue]==key],
                    color=palette[key],
                    **kwargs
                )
            )
    else:
        sns.regplot(
                    x=x,
                    y=y,
                    data=data,
                    color=palette,
                    **kwargs
                    )
        
    return regplots

def get_clinical_and_error(df_merge, clinical="Neutrophil", error="error"):
    df_merge = df_merge[["Project-ID", clinical, error]]
    df_merge = df_merge.dropna(axis=0)
    clinical_values = df_merge[clinical].to_list()
    error_values = df_merge[error].to_list()

    return clinical_values, error_values
    
def get_clinical_and_error_divided_by_CRP_group(df_merge, group="CRP_Group", clinical="Neutrophil", error="error"):
    df_merge = df_merge[["Project-ID", group, clinical, error]]
    df_merge = df_merge.dropna(axis=0)
    clincal_by_group = df_merge.groupby(group)[clinical].apply(list)
    error_by_group = df_merge.groupby(group)[error].apply(list)
    if len(clincal_by_group) == 2:
        list_clinical_high = clincal_by_group.loc["High"]
        list_clinical_low = clincal_by_group.loc["Low"]
        list_error_high = error_by_group.loc["High"]
        list_error_low = error_by_group.loc["Low"]
        
        return list_clinical_high, list_clinical_low, list_error_high, list_error_low
    
    else:
        raise Exception

def calculate_correlation(x, y):
    r, p = pearsonr(x, y)

    return r, p

def calculate_regression_slope(x, y):
    X = sm.add_constant(x)
    Y = np.reshape(y, (-1, 1))
    results = sm.OLS(Y, X).fit()
    slope = results.params[1]
    pvalue = results.pvalues[1]
    
    return slope, pvalue

def compare_regression_slopes_by_CRP_group(df_merge, var_x = "Neutrophil", var_y = "error", group="CRP_Group"):
    df_merge = df_merge[["Project-ID", group, var_x, var_y]]
    df_merge = df_merge.dropna(axis=0)
    df_merge['group_encoded'] = (df_merge[group] == "Low").astype(int)
    df_merge['interaction'] = df_merge['group_encoded'] * df_merge[var_x]
    X = sm.add_constant(df_merge[[var_x, 'group_encoded', 'interaction']])
    model = sm.OLS(df_merge[var_y], X).fit()
    
    return model
