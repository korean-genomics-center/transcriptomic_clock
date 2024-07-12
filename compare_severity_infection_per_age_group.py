# %%
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import Counter
import matplotlib.font_manager as fm
from matplotlib.cm import get_cmap
from matplotlib import gridspec

def group_patients_by_crp_level(dict_crp):
    dict_crp_group = dict()
    for sample, crp in dict_crp.items():
        if crp <= 1:
            group = "Low"
        elif crp > 1:
            group = "High"
        
        dict_crp_group[sample] = group

    return dict_crp_group        

def get_column_crp_group(df_crf, dict_crp_group):
    list_crp_group = list(map(lambda x: dict_crp_group.get(x, None), df_crf["Individual"]))

    return list_crp_group

def hue_regplot(data, x, y, hue, palette=None, **kwargs):
    regplots = list()
    levels = data[hue].unique()
    if palette is None:
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
    
    return regplots

# %%
