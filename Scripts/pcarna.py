# %%
import math
import os
import string

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
from matplotlib.colors import rgb2hex
from matplotlib.patches import Ellipse
from matplotlib.pyplot import get_cmap
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

class PCARna():
    def __init__(self, list_drop_samples=[], list_select_samples=[]):
        self.list_drop_samples = list_drop_samples
        self.list_select_samples = list_select_samples
        self.colsampleid = "ID"
        self.sep = "\t"
        self.n_components = 20

    def shorten_meta_info(self, df_meta):
        df_meta_shortened = df_meta.copy()
        for col in list(df_meta_shortened.columns)[3:]:
            meta_info = df_meta_shortened[col]
            new_meta_info = list()
            for elem in meta_info:
                if len(str(elem)) > 20:
                    new_elem = str(elem)[:21] + ".."
                else:
                    new_elem = str(elem)
                new_meta_info.append(new_elem)
            df_meta_shortened[col] = new_meta_info
        
        return df_meta_shortened

    def run_pca(self, df_train, df_test=None, std=True):
        train_values = df_train.values
        if std:
            scaler = StandardScaler().fit(df_train)
            train_values = scaler.transform(df_train)
        pca = PCA(n_components=self.n_components).fit(train_values)
        if df_test is None:
            pca_exp = pca.transform(train_values)
        else:
            test_values = df_test.values
            pca_exp = pca.transform(test_values)

        return pca, pca_exp

    def draw_scree_plot(self, pca, outdir, test_name):
        # Draw scree plot
        y = np.cumsum(pca.explained_variance_ratio_)
        x = list(1, range(len(y)), 1)
        plt.plot(x, y)
        plt.xlabel('number of components')
        plt.ylabel('cumulative explained variance')
        plt.show()
        outfigname = os.path.join(outdir, f"pca_scree_plot_{test_name}.png")
        plt.savefig(outfigname, dpi=600)
        plt.close()

    def merge_pcadata_and_metadata(self, pca_exp, df_meta, colsampleid="Project-ID"):
        # merge pcadata and metadata
        df_meta["Sample_Age_Group"] = df_meta["Sample_Age"].apply(lambda x: int(str(x)[0] + "0"))
        if colsampleid not in list(df_meta.columns):
            list_samples = list(df_meta.index)
        else:
            list_samples = df_meta[colsampleid].to_list()
        list_PC = list(map(lambda x: f"PC{x}", list(range(1, self.n_components+1, 1))))
        df_pca_exp = pd.DataFrame(pca_exp, columns=list_PC, index=list_samples).reset_index(drop=False).rename(columns={"index":colsampleid})

        df_pca_meta = pd.merge(df_pca_exp, df_meta, how="inner", on=colsampleid)

        return df_pca_meta

    def get_loading_vector(self, df_exp_dropna, pca, outdir, project_info, colsampleid="ID"):
        # Get loading vector
        list_feature_name = df_exp_dropna[colsampleid].to_list()
        list_pc = [np.transpose(pca.components_[i]) for i in range(pca.n_components_)]

        df_load = pd.DataFrame()
        df_load["Gene_Symbol"] = list_feature_name
        for ind, pc in enumerate(list_pc, 1):
            df_load[f"PC{ind}"] = pc

        df_load.to_csv(os.path.join(outdir, f"pca_loading_vector_{project_info}.tsv"))

        return df_load

    def extract_top_pc_features(self, df_exp_transposed, pca, num_top=50):
        # Extract genes involved in making components
        # https://stackoverflow.com/questions/50796024/feature-variable-importance-after-a-pca-analysis
        # name of features
        feature_names = list(df_exp_transposed.columns)
        # number of components
        n_pcs= pca.components_.shape[0]
        # get the index of the most important feature on EACH component
        most_important = [np.abs(pca.components_[i]).argsort()[::-1][:num_top] for i in range(n_pcs)]

        for _, list_idx in enumerate(list(map(list, most_important)), 1):
            list_features_top = list()
            for idx in list_idx:
                feat = feature_names[idx]
                list_features_top.append(feat)
        
        return list_features_top

    def confidence_ellipse(self, x, y, ax, n_std=3.0, facecolor='none', **kwargs):
        """
        Create a plot of the covariance confidence ellipse of *x* and *y*.

        Parameters
        ----------
        x, y : array-like, shape (n, )
            Input data.

        ax : matplotlib.axes.Axes
            The axes object to draw the ellipse into.

        n_std : float
            The number of standard deviations to determine the ellipse's radiuses.

        **kwargs
            Forwarded to `~matplotlib.patches.Ellipse`

        Returns
        -------
        matplotlib.patches.Ellipse
        """
        if x.size != y.size:
            raise ValueError("x and y must be the same size")

        cov = np.cov(x, y)
        pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
        # Using a special case to obtain the eigenvalues of this
        # two-dimensional dataset.
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)
        ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                        facecolor=facecolor, **kwargs)

        # Calculating the standard deviation of x from
        # the squareroot of the variance and multiplying
        # with the given number of standard deviations.
        scale_x = np.sqrt(cov[0, 0]) * n_std
        mean_x = np.mean(x)

        # calculating the standard deviation of y ...
        scale_y = np.sqrt(cov[1, 1]) * n_std
        mean_y = np.mean(y)

        transf = transforms.Affine2D() \
            .rotate_deg(45) \
            .scale(scale_x, scale_y) \
            .translate(mean_x, mean_y)

        ellipse.set_transform(transf + ax.transData)
        return ax.add_patch(ellipse)
    
    def rgb_to_cmyk(self, r, g, b, CMYK_SCALE = 100, RGB_SCALE = 1):
        if (r, g, b) == (0, 0, 0):
            # black
            return 0, 0, 0, CMYK_SCALE

        # rgb [0,255] -> cmy [0,1]
        c = 1 - r / RGB_SCALE
        m = 1 - g / RGB_SCALE
        y = 1 - b / RGB_SCALE

        # extract out k [0, 1]
        min_cmy = min(c, m, y)
        c = (c - min_cmy) / (1 - min_cmy)
        m = (m - min_cmy) / (1 - min_cmy)
        y = (y - min_cmy) / (1 - min_cmy)
        k = min_cmy

        # rescale to the range [0,CMYK_SCALE]
        return c * CMYK_SCALE, m * CMYK_SCALE, y * CMYK_SCALE, k * CMYK_SCALE

    def cmyk_to_rgb(self, c, m, y, k, cmyk_scale = 100, rgb_scale=1, alpha = 1):
        r = rgb_scale * (1.0 - c / float(cmyk_scale)) * (1.0 - k / float(cmyk_scale))
        g = rgb_scale * (1.0 - m / float(cmyk_scale)) * (1.0 - k / float(cmyk_scale))
        b = rgb_scale * (1.0 - y / float(cmyk_scale)) * (1.0 - k / float(cmyk_scale))
        return r, g, b, alpha

    def draw_pc_biplot(self,
                        df_pca_input, 
                        list_meta_columns, 
                        pca_obj, 
                        pcx_num, 
                        pcy_num, 
                        ncol=2, 
                        dict_legend_ncol=dict(),
                        plot_linewidth=1,
                        figsize=(5, 5),
                        xlim=(-20, 80),
                        ylim=(-30, 20),
                        colormap="Spectral"):
    
        pcx = f"PC{pcx_num}"
        pcy = f"PC{pcy_num}"
        nrow = math.ceil(len(list_meta_columns) / ncol)
        fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=figsize)
        axes_flat = axes.flatten()
        
        for i, (meta, ax) in enumerate(zip(list_meta_columns, axes_flat)):
            fig_letter = list(string.ascii_uppercase)[i]

            # Check if the grouping by metadata column has data
            try:
                val_groupby_meta_pcx = df_pca_input.groupby(meta)[pcx].apply(np.array).tolist()
                val_groupby_meta_pcy = df_pca_input.groupby(meta)[pcy].apply(np.array).tolist()
                
                # Ensure the grouped data is not empty
                if not val_groupby_meta_pcx or not val_groupby_meta_pcy:
                    print(f"Warning: No data found for metadata column {meta}. Skipping.")
                    continue  # Skip this metadata column if no data
                
            except Exception as e:
                print(f"Error processing meta column {meta}: {e}")
                continue  # Skip this subplot if data is missing or invalid
            
            # Color map and transformation
            cmap = get_cmap(colormap, len(val_groupby_meta_pcx))
            list_colors = [cmap(ind) for ind in range(len(val_groupby_meta_pcx))]
            dim_level = 10
            list_colors_cmyk = [self.rgb_to_cmyk(rgb[0], rgb[1], rgb[2]) for rgb in list_colors]
            list_colors_cmyk_dim = [(cmyk[0], cmyk[1], cmyk[2], cmyk[3] + dim_level) for cmyk in list_colors_cmyk]
            list_colors_rgb_dim = [self.cmyk_to_rgb(cmyk[0], cmyk[1], cmyk[2], cmyk[3]) for cmyk in list_colors_cmyk_dim]

            # Create legends
            list_legends = df_pca_input.groupby(meta)[pcy].apply(np.array).index.tolist()
            list_legends = [int(float(x)) if str(x).endswith(".0") else x for x in list_legends]
            legend_ncol = dict_legend_ncol.get(meta, 1)

            # Plot each group
            for (x, y, color, legend) in zip(val_groupby_meta_pcx, val_groupby_meta_pcy, list_colors_rgb_dim, list_legends):
                if len(x) < 2 or len(y) < 2:  # Ensure enough points to plot an ellipse
                    print(f"Warning: Not enough data points to plot ellipse for {meta}. Skipping confidence ellipse.")
                    ax.scatter(x, y, c=color, s=plot_linewidth*10, label=legend, alpha=0.3, zorder=-11)
                else:
                    ax.scatter(x, y, c=color, s=plot_linewidth*10, label=legend, alpha=0.3, zorder=-11)
                    self.confidence_ellipse(x, y, ax, n_std=3.0, edgecolor=color, facecolor='none', linewidth=plot_linewidth*2.0, alpha=0.5)
                ax.scatter(np.mean(x), np.mean(y), c=color, edgecolor="k", linewidth=plot_linewidth*0.5, marker="*", s=plot_linewidth*70, alpha=1, zorder=999)
            
            # Set plot limits and labels
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_title(meta, fontsize=plt.rcParams["font.size"]+1, fontweight="bold", pad=2)
            ax.set_ylabel(f"principal component {pcy_num}\n({round(pca_obj.explained_variance_ratio_[pcy_num - 1] * 100, 1)}%)", 
                        fontsize=plt.rcParams["font.size"]+1, labelpad=2)
            ax.set_xlabel(f"principal component {pcx_num}\n({round(pca_obj.explained_variance_ratio_[pcx_num - 1] * 100, 1)}%)", 
                        fontsize=plt.rcParams["font.size"]+1, labelpad=2)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

            # Add legend inside the plot without overlapping with data points
            if len(list_legends) > 0:  # Ensure that legends are added only when valid
                legend = ax.legend(loc="upper right", fontsize=plt.rcParams["font.size"]-2,
                                ncol=legend_ncol, markerscale=plot_linewidth, handletextpad=1.0,
                                bbox_to_anchor=(1.0, 1.0), frameon=True)
                legend.get_frame().set_alpha(1.0)  # Set legend box alpha to make it opaque
                plt.setp(legend.get_texts(), weight='bold')

                # Make legend dots completely opaque
                for handle in legend.legendHandles:
                    handle.set_alpha(1)
            
            # Annotate figure lettering outside the plot, on the far left to avoid overlap
            ax.annotate(fig_letter,
                        xy=(-0.2, 1.04),
                        xycoords='axes fraction',
                        xytext=(0, 0),
                        textcoords='offset points',
                        size=plt.rcParams["font.size"]+5, ha='left', va='center',
                        fontweight="bold", color="black")
        
        # Turn off axes for any remaining subplots if necessary
        for extra_ind in range(len(list_meta_columns), len(axes_flat)):
            axes_flat[extra_ind].axis("off")
        
        return fig, axes

    def get_list_selected_genes(self, infile):
        list_selected_genes = list()
        with open(infile, mode="r") as fr:
            _skiprow = fr.readline()
            for line in fr:
                record = line.rstrip("\n").split("\t")
                gene = record[-1]
                list_selected_genes.append(gene)

        return list_selected_genes

    def save_selected_genes(self, list_selected_genes, outdir):
        fileingenename = os.path.join(outdir, "features_ingenename.txt")
        with open(fileingenename, mode="w") as fw:
            for gene in list_selected_genes:
                fw.write(gene + "\n")

    def get_important_features(self, file_features, colfeat="Selected_Features"):
        list_features = list()
        with open(file_features, mode="r") as fr:
            list_header = fr.readline().rstrip("\n").split("\t")
            idx_feat = list_header.index(colfeat)
            for line in fr:
                record = line.rstrip().split()
                gene = record[idx_feat]
                list_features.append(gene)

        return list_features

    def get_num_original_feature(self, dir_scaled_exp):
        file_original_features = os.path.join(os.path.dirname(dir_scaled_exp), "regression_results.txt")
        list_original_features = self.get_important_features(file_original_features)
        num_original_features = len(list_original_features)

        return num_original_features

    def get_path_file_features(self, dir_scaled_exp):
        list_filename_scaled_exp = os.listdir(dir_scaled_exp)
        list_path_scaled_exp = list(map(lambda x: os.path.join(dir_scaled_exp, x), list_filename_scaled_exp))

        return list_path_scaled_exp

    def get_select_random_features(self, scaled_exp, num_original_features):
        with open(scaled_exp, mode="r") as fr:
            list_header = fr.readline().rstrip("\n").split("\t")
            np.random.seed(1)
            list_rand_features = list(np.random.choice(list_header[1:-1], size=num_original_features, replace=False))
            list_new_header = [list_header[0]] + list_rand_features + [list_header[-1]]

        return list_new_header
# %%
