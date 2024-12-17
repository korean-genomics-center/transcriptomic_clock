# %%
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from transcriptomic_clock import Clock

# %%
def make_df_perfomance_evaluation(workdir, testing_data, list_corr_thres):
    dict_annot_age_eval_plot = dict()
    for corr_thres in list_corr_thres:
        dir_trained_model = f"{workdir}/corr_{corr_thres}_above_pval_0.05_below"
        file_testing_data = f"{dir_trained_model}/{testing_data}"
        clock = Clock(dir_trained_model, file_testing_data, col_sampleid="Project-ID", col_y="Sample_Age")
        print(f"----{corr_thres}----")
        try:
            clock.make_dataframe_clock_features()
            dict_annot_corr_thres = clock.get_annotation_prediction_evaluation_plot()
            dict_annot_age_eval_plot[corr_thres] = dict_annot_corr_thres
        except Exception as e:
            print(e)
            print(f"No data for {corr_thres}... Skipping...")
    df_annot_age_eval_plot = pd.DataFrame(dict_annot_age_eval_plot)

    return df_annot_age_eval_plot

# %% [draw_plot_performance_profile]
cm = 1/2.54
path = "./Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
width = 10*cm
height = 5*cm
figsize = (width, height)
plot_linewidth = 3

WORKDIR = "./LASSO_INFO/standardized/healthy_illumina"
dir_trained_model = f"{WORKDIR}/corr_0.3_above_pval_0.05_below"
list_corr_thres = list(map(lambda x: round(float(x), 2), list(np.arange(0.3, 0.51, 0.01))))
clock = Clock(dir_trained_model,  )
df_annot_age_eval_plot_train = make_df_perfomance_evaluation(WORKDIR, "training_dataset_std_scaled.txt", list_corr_thres)
df_annot_age_eval_plot_test1 = make_df_perfomance_evaluation(WORKDIR, f"testing_dataset_selected_healthy_std_scaled.txt", list_corr_thres)
col_train, mae_train, corr_train, list_num_input_genes_train, list_num_features_train = clock.summarize_performance_statistics(df_annot_age_eval_plot_train)
col_test, mae_test, corr_test, list_num_input_genes_test, list_num_features_test = clock.summarize_performance_statistics(df_annot_age_eval_plot_test1)
fig, ax1, ax2 = clock.draw_plot_performance_profile(col_train, 
                                    col_test, 
                                    mae_train, 
                                    mae_test, 
                                    corr_train, 
                                    corr_test, 
                                    list_num_features_train, 
                                    list_num_input_genes_train,
                                    plot_linewidth=plot_linewidth, 
                                    figsize=figsize)
# Annotate
ax1.annotate('Thresholds\n\nNo. LASSO Features\n\nNo. Input Genes', 
             xy=(0.01, 1.3), 
             xycoords=ax1.get_yaxis_transform(),
             xytext=(0, 0), 
             textcoords="offset points", 
             ha="right", 
             va="center", 
             fontsize=plt.rcParams["font.size"], 
             weight="bold")

# Grid and labels
ax1.grid(axis="both", linestyle="dashed", linewidth=plot_linewidth-0.5)
ax1.set_ylabel('MAE', color='darkred', fontsize=plt.rcParams["font.size"]+2, labelpad=10)
ax2.set_ylabel("Pearson's r", color='darkblue', fontsize=plt.rcParams["font.size"]+2, labelpad=10)  # we already handled the x-label with ax1

# Add vertical lines
plt.axvline(x=0.30, color="grey", linewidth=plot_linewidth)
plt.axvline(x=0.37, color="grey", linewidth=plot_linewidth)

# Legend
legend_elements = [
                Line2D([0], [0], color='k', linestyle="dotted", linewidth=plot_linewidth, label='Training'),
                Line2D([0], [0], color='k', linestyle="solid", linewidth=plot_linewidth, label='Validation')
                ]
plt.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(1, 1), prop={'size': plt.rcParams["font.size"]}, frameon=False)

# Final adjustments
plt.margins(x=0.1)
plt.tight_layout()

# Show the plot
plt.show()
plt.close()
# %%
