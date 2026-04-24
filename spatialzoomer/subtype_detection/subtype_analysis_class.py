from .CFSFDP import CFSFDP
from .utils import calculate_p_value, ward_linkage
from .visualization import generate_colors
from .spatial_neighborhood import calculate_neighbor_composition

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import seaborn as sns
from kneed import KneeLocator
from matplotlib.colors import to_hex
import os

class SubclusterAnalysis:
    def __init__(self, adata, save_path, cluster_key='leiden_Raw_res1'):
        self.adata = adata
        self.umap_keys=self._get_umap_keys()
        self.cluster_key = cluster_key
        self.adata.obs[self.cluster_key] = self.adata.obs[self.cluster_key].astype('category')
        self.adata.obs_names = self.adata.obs_names.astype(str)

        self.n_subcluster = None
        self.selected_cluster = None
        self.adata_selected = None
        self.adata_selected_plot = None
        self.optimal_scale = None
        self.var_dict = None
        self.k_dict = None

        self.save_path = save_path
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        self._set_visual_config()
    
    def _get_umap_keys(self):
        obsm_key = list(self.adata.obsm.keys())
        umap_keys = [key for key in obsm_key if key.startswith('X_umap_scale')]
        umap_keys = sorted(umap_keys, key=lambda x: float(x.split('scale')[-1]))
        umap_keys = [umap_key for umap_key in umap_keys if float(umap_key.split('scale')[-1]) <= 10.0]
        return umap_keys
    
    def _set_visual_config(self):
        self.dpi=300
        self.new_colors = [to_hex(c) for c in generate_colors(200)]
        self.colors_3 = ['#CD69C9', '#9ACD32', "#FF7F00"] + self.new_colors
        config = {
            "font.family":'sans-serif',
            "font.size": 7,
            "mathtext.fontset":'stix',
            "font.serif": ['MS Arial'],
        }
        rcParams.update(config)
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
    
    def calculate_metrics(self, nonepi=True, plot_clusters=True, figsize = (10, 4)):
        """
        Plot covariance-scale curves and covariance values.
        input adata: AnnData object, with umaps at raw ("X_umap_Raw") and different scales (for example "X_umap_scale0.01")
        input cluster_key: key of the clustering result in adata.obs
        input save_path: path to save the figures
        """

        ######################### plot umap #########################
        # set as categorical
        cell_type_counts = self.adata.obs[self.cluster_key].value_counts()
        sorted_cell_types = cell_type_counts.index.tolist()
        color_mapping = {cell_type: self.new_colors[i] for i, cell_type in enumerate(sorted_cell_types)}
        sorted_colors = [color_mapping[cell_type] for cell_type in self.adata.obs[self.cluster_key].cat.categories]
        self.adata.uns[self.cluster_key + '_colors'] = sorted_colors
        self.adata.obsm['X_umap'] = self.adata.obsm['X_umap_Raw']

        if plot_clusters:
            fig, axs = plt.subplots(1, 2, figsize=figsize, dpi=self.dpi, 
                                    gridspec_kw={'width_ratios': [1.2, 1.8]})
            sc.pl.umap(self.adata, 
                    color=self.cluster_key, 
                    size=5, 
                    ax=axs[0], 
                    show=False, 
                    title='Cell clusters')
            for collection in axs[0].collections:
                collection.set_rasterized(True) 
            # hide legend
            axs[0].get_legend().remove()

            sq.pl.spatial_scatter(
                self.adata, 
                library_id="spatial", 
                shape=None, 
                color= self.cluster_key,
                wspace=0.1,
                frameon=False,
                size=1,
                figsize=(18, 18),
                dpi=self.dpi,
                outline=False,
                img=False,
                marker='.',
                ax=axs[1],
                title='Cell clusters',
            )
            for collection in axs[1].collections:
                collection.set_rasterized(True)
            plt.tight_layout()
            plt.savefig(self.save_path + "/Spatial_clusters.pdf", bbox_inches='tight', dpi=self.dpi)
            plt.show()
            plt.close()

        
        ######################### plot covariance-scale curves #########################        
        # normalize umap
        self.adata.obsm['scaled_X_umap_Raw'] = (self.adata.obsm['X_umap_Raw'] - self.adata.obsm['X_umap_Raw'].min()) / (self.adata.obsm['X_umap_Raw'].max() - self.adata.obsm['X_umap_Raw'].min())
        for umap_scale in self.umap_keys:
            # X.obsm['umap_scale']归一化
            self.adata.obsm['scaled'+umap_scale] = (self.adata.obsm[umap_scale] - self.adata.obsm[umap_scale].min()) / (self.adata.obsm[umap_scale].max() - self.adata.obsm[umap_scale].min())
        
        # loop over clusters and scales to calculate variance
        nmf_clusters = self.adata.obs[self.cluster_key].unique()
        var_dict = {}
        k_dict = {}
        for nmf_cluster in nmf_clusters:
            var_list = []
            for umap_scale in self.umap_keys:
                X = self.adata[self.adata.obs[self.cluster_key] == nmf_cluster].obsm['scaled'+umap_scale]
                var = np.var(X, axis=0).mean()
                var_list.append(var)
            k = np.polyfit([float(x.split('scale')[-1]) for x in self.umap_keys], var_list, 1)[0]
            var_dict[nmf_cluster] = var_list
            k_dict[nmf_cluster] = k

        if nonepi:
            var_dict = {k: v for k, v in var_dict.items() if "Tumor Cells" not in k 
                            and "Epithelial" not in k 
                            and "Epithelium" not in k
                            and "AT2" not in k 
                            and "Unassigned" not in k 
                            and "Ciliated" not in k
                            and "Malignant" not in k}
            k_dict = {k: v for k, v in k_dict.items() if "Tumor Cells" not in k 
                            and "Epithelial" not in k 
                            and "Epithelium" not in k
                            and "AT2" not in k 
                            and "Unassigned" not in k 
                            and "Ciliated" not in k
                            and "Malignant" not in k}

        k_dict = dict(sorted(k_dict.items(), key=lambda x: x[1], reverse=False))
        clusters = list(k_dict.keys())
        values = list(k_dict.values())


        fig, axs = plt.subplots(1, 2, figsize=(10, 4), dpi=self.dpi)  # 调整画布大小以适应布局

        # plot variance-scale curves for each cluster on the left
        for nmf_cluster, var_list in var_dict.items():
            axs[0].plot(var_list, label=nmf_cluster, color=color_mapping[nmf_cluster])
        # 移除图例和数值
        axs[0].set_xticks(range(len(self.umap_keys)))
        axs[0].set_xticklabels([x.split('scale')[-1] for x in self.umap_keys])  # x轴标签

        axs[0].set_xlabel("Scales")
        axs[0].set_ylabel("Variance of sample distribution")
        axs[0].set_title("Variance-scale curves")

        # bar plot of k values for each cluster on the right
        right_ax = axs[1]
        right_ax.barh(clusters, values, height=0.8, color=[color_mapping[cluster] for cluster in clusters])

        right_ax.set_yticks(range(len(clusters)))
        right_ax.set_yticklabels(clusters)
        right_ax.set_title("Rank by slope values")


        plt.tight_layout()
        plt.savefig(self.save_path + "/Spatially_dependent_clusters_rank.pdf", dpi=self.dpi)
        plt.show()
        plt.close()
        
        print("Suggested cluster: ", clusters[-1])    # the cluster with the largest k value
        self.var_dict = var_dict
        self.k_dict = k_dict

    # def kde_plot(self):
    #     """
    #     Plot KDE plots for each cluster.
    #     input adata: AnnData object, with umaps at raw ("X_umap_Raw") and different scales (for example "X_umap_scale0.01")
    #     input cluster_key: key of the clustering result in adata.obs
    #     input save_path: path to save the figures
    #     """
    #     if self.selected_cluster is None:
    #         print("Select a cluster first. SubclusterAnalysis.selected_cluster = 'cluster_name'")
    #         return
    #     else:
    #         print("Selected cluster:", self.selected_cluster)
        
    #     adata_selected = self.adata[self.adata.obs[self.cluster_key] == self.selected_cluster]
    #     self.adata_selected = adata_selected

    #     fig, axes = plt.subplots(int((len(self.umap_keys)+2)/3), 3, figsize=(9, len(self.umap_keys)), dpi=self.dpi)

    #     for key in self.umap_keys:
    #         # subsample for large datasets during kde plot
    #         if adata_selected.n_obs > 10000:
    #             adata_selected_subsampled = sc.pp.subsample(adata_selected, n_obs=10000, random_state=0, copy=True)
    #         else:
    #             adata_selected_subsampled = adata_selected

    #         X = adata_selected_subsampled.obsm[key]
    #         ax = axes.flatten()[self.umap_keys.index(key)]
    #         ax.scatter(X[:, 0], X[:, 1], color='gray', s=1, marker='.', alpha=0.5)
    #         for collection in ax.collections:
    #             collection.set_rasterized(True)  # 将散点设置为栅格化
    #         # KDE plot
    #         umap_df = pd.DataFrame(X, columns=['UMAP1', 'UMAP2'])
    #         sns.kdeplot(data=umap_df, x='UMAP1', y='UMAP2', ax=ax, bw_adjust=1.5,
    #                     cmap='Blues', fill=True, thresh=0.01, levels=8, alpha=0.7)
    #         ax.set_title(key)
    #     plt.tight_layout()
    #     plt.savefig(self.save_path + "/umaps_kde.pdf", dpi=self.dpi)
    #     plt.show()
    #     plt.close()

    def kde_plot(self):
        """
        Plot KDE plots for each cluster.
        input adata: AnnData object, with umaps at raw ("X_umap_Raw") and different scales (for example "X_umap_scale0.01")
        input cluster_key: key of the clustering result in adata.obs
        input save_path: path to save the figures
        """
        if self.selected_cluster is None:
            print("Select a cluster first. SubclusterAnalysis.selected_cluster = 'cluster_name'")
            return
        else:
            print("Selected cluster:", self.selected_cluster)
        
        adata_selected = self.adata[self.adata.obs[self.cluster_key] == self.selected_cluster]
        self.adata_selected = adata_selected

        n_keys = len(self.umap_keys)
        n_cols = 3
        n_rows = int(np.ceil(n_keys / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(3 * n_cols, 3 * n_rows), dpi=self.dpi)
        axes = axes.flatten()

        for i, key in enumerate(self.umap_keys):
            # subsample for large datasets during kde plot
            if adata_selected.n_obs > 10000:
                adata_selected_subsampled = sc.pp.subsample(adata_selected, n_obs=10000, random_state=0, copy=True)
            else:
                adata_selected_subsampled = adata_selected

            X = adata_selected_subsampled.obsm[key]
            ax = axes[i]
            ax.scatter(X[:, 0], X[:, 1], color='gray', s=1, marker='.', alpha=0.5)
            for collection in ax.collections:
                collection.set_rasterized(True)  # 栅格化以减小 PDF 大小

            # KDE plot
            umap_df = pd.DataFrame(X, columns=['UMAP1', 'UMAP2'])
            sns.kdeplot(data=umap_df, x='UMAP1', y='UMAP2', ax=ax, bw_adjust=1.5,
                        cmap='Blues', fill=True, thresh=0.01, levels=8, alpha=0.7)
            ax.set_title(key)

        # 删除多余的子图
        for j in range(n_keys, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plt.savefig(self.save_path + "/UMAP_with_KDE.pdf", dpi=self.dpi)
        plt.show()
        plt.close()

    def find_optimal_scale(self):
        if self.var_dict is None:
            print("Running calculate_metrics() first...")
            self.calculate_metrics()

        var_list = self.var_dict[self.selected_cluster]
        kl = KneeLocator(
            range(len(var_list)), 
            var_list, 
            curve='concave', 
            direction='increasing', 
            online=True,
            S=0.1,
        )
        optimal_scale = self.umap_keys[kl.elbow]
        optimal_scale = optimal_scale.split('scale')[-1]
        print("Optimal scale for cluster", self.selected_cluster, "is", optimal_scale)
        self.optimal_scale = optimal_scale

    def detect_subclusters(self, n_clusters=2, plot_all_scales=True, calculate_ward=True, plot_spatial=True):
        if self.selected_cluster is None:
            print("Select a cluster first. SubclusterAnalysis.selected_cluster = 'cluster_name'")
            return
        if self.adata_selected is None:
            adata_selected = self.adata[self.adata.obs[self.cluster_key] == self.selected_cluster]
        else:
            adata_selected = self.adata_selected

        self.n_subcluster = n_clusters
        adata_selected.obsm['X_umap'] = adata_selected.obsm['X_umap_scale' + str(self.optimal_scale)]
        clustering_model = CFSFDP(n_clusters=n_clusters)
        res = clustering_model.fit_predict(adata_selected.obsm['X_umap'])
        adata_selected.obs['cluster_label'] = list(res.astype(str))

        adata_selected_plot = adata_selected[adata_selected.obs['cluster_label'] != '-1']
        adata_selected_plot.uns['cluster_label_colors'] = self.colors_3[0:self.n_subcluster] #需要修改，动态调整颜色


        if plot_spatial:
            fig, axs = plt.subplots(1, 2, figsize=(10, 4), dpi=self.dpi, 
                                    gridspec_kw={'width_ratios': [1.2, 1.8]})
            sc.pl.umap(adata_selected_plot, 
                       color=['cluster_label'], 
                       title='Subclusters at scale' + str(self.optimal_scale), 
                       show=False, 
                       ax=axs[0], 
                       size=5)
            for collection in axs[0].collections:
                collection.set_rasterized(True) 
            # hide legend
            axs[0].get_legend().remove()

            sq.pl.spatial_scatter(
                adata_selected_plot, 
                library_id="spatial", 
                shape=None, 
                color= ['cluster_label'],
                wspace=0.1,
                frameon=False,
                size=1,
                figsize=(18, 18),
                dpi=self.dpi,
                outline=False,
                img=False,
                marker='.',
                ax=axs[1],
                title='Cluster_spatial_' + self.selected_cluster + '_scale' + str(self.optimal_scale),
            )
            for collection in axs[1].collections:
                collection.set_rasterized(True)
            plt.tight_layout()
            plt.savefig(self.save_path + "/Spatial_" + str(self.selected_cluster) + "_scale" + str(self.optimal_scale) + ".pdf", dpi=self.dpi)
            plt.show()
            plt.close()

        umap_keys = self.umap_keys
        if calculate_ward:
            ward_list = dict()
            for key in umap_keys:
                X = adata_selected_plot.obsm['scaled'+key]
                y = adata_selected_plot.obs['cluster_label']
                ward_list[key] = ward_linkage(X, y)

        if plot_all_scales:
            n_plots = len(umap_keys)
            n_rows = int((n_plots + 2) / 3)
            fig, axes = plt.subplots(n_rows, 3, figsize=(12, n_rows * 3.6), dpi=self.dpi)
            axes = axes.flatten()
            for idx, key in enumerate(umap_keys):
                adata_selected_plot.obsm['X_umap'] = adata_selected_plot.obsm[key]
                ax = axes[idx]
                sc.pl.umap(adata_selected_plot, 
                        size=5,
                        color=['cluster_label'], 
                        title=key, 
                        ax=ax, 
                        show=False)
                for collection in ax.collections:
                    collection.set_rasterized(True)
                if calculate_ward:
                    ax.set_title(key.replace("X_umap_s", "S") + '\n ward: %.2f' % ward_list[key])
                else:
                    ax.set_title(key)

            for j in range(len(umap_keys), len(axes)):
                fig.delaxes(axes[j])

            plt.tight_layout()
            plt.savefig(self.save_path + "/UMAP_cluster_scale" + str(self.optimal_scale) + ".pdf", dpi=self.dpi)
            plt.show()
            plt.close()
        
        self.adata_selected = adata_selected
        self.adata_selected_plot = adata_selected_plot    # without outliers

    def _deg2df_noise_filtered(self, deg_groups, subcluster_name, deg_current_cluster_df):
        result_df_subcluster = pd.DataFrame({key: deg_groups[key][subcluster_name] for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
        lfc_seleccted_cluster = deg_current_cluster_df.loc[list(result_df_subcluster['names'])]['logfoldchanges']
        result_df_subcluster['logfoldchanges_celltype'] = list(lfc_seleccted_cluster)
        print("saving DEG results for subcluster", subcluster_name)
        result_df_subcluster.to_csv(self.save_path + f"/DEGs_"+str(subcluster_name)+".csv")
        result_df_subcluster = result_df_subcluster[result_df_subcluster['logfoldchanges_celltype'] > 0]
        return result_df_subcluster

    def deg_analysis_subclusters(self, dotplot=True):
        if self.adata_selected_plot is None:
            print("Run detect_subclusters() first...")
            return
        sc.tl.rank_genes_groups(self.adata_selected_plot, 'cluster_label', method='wilcoxon')

        sc.tl.rank_genes_groups(self.adata, self.cluster_key, method='wilcoxon')
        result = self.adata.uns['rank_genes_groups']
        deg_current_cluster_df = pd.DataFrame({key: result[key][self.selected_cluster] for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
        deg_current_cluster_df.index = deg_current_cluster_df['names']
        deg_df_list = []
        for subcluster_name in np.unique(self.adata_selected_plot.obs['cluster_label']):
            result_df_subcluster = self._deg2df_noise_filtered(self.adata_selected_plot.uns['rank_genes_groups'], subcluster_name, deg_current_cluster_df)
            deg_df_list.append(result_df_subcluster)

        if dotplot:
            genes = []
            for deg_df in deg_df_list:
                genes = genes + list(deg_df.iloc[:6,0])
            # Generate the stacked violin plot
            green_cmap = ListedColormap(plt.cm.Greens(np.linspace(0.2, 1, 256)))
            sc.pl.dotplot(
                self.adata_selected_plot,
                var_names=genes,
                groupby='cluster_label',
                #inner='box',  # Add box plots within the violins
                figsize=(4.5, 1.8),  # Larger size for better visibility
                layer=None,
                swap_axes=False,  # Stacks violins horizontally
                cmap=green_cmap,  # Apply the green color map
                show=False,  # Delay showing to customize further
                standard_scale=None,  # Scale the violins by width
            )

            # Adjust legend and other aesthetics
            ax = plt.gca()
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
            plt.tight_layout()  # Adjust layout to avoid overlap
            plt.savefig(self.save_path + "/Dotplot_DEGs.pdf", dpi=self.dpi)
            plt.show()
            plt.close()
    
    def _subcluster_label_merge(self, add_key):
        celltype_labels = self.adata.obs[self.cluster_key]
        self.adata.obs[add_key] = [celltype_labels[i] if i not in self.adata_selected.obs_names else (celltype_labels[i] + '_cluster' + self.adata_selected.obs['cluster_label'][i]) for i in celltype_labels.index]
        print("Saving cell type labels with subcluster labels...")
        pd.DataFrame({'cell_id': self.adata.obs['cell_id'], 'group': self.adata.obs[add_key]}).to_csv(self.save_path + "/Cell_clusters.csv", index=False)
        adata_plot = self.adata[self.adata.obs[add_key] != self.selected_cluster+'_cluster-1']
        return adata_plot

    def neighborhood_analysis(self, k=10):
        if self.adata_selected_plot is None:
            print("Run detect_subclusters() first...")
            self.detect_subclusters()
            return
        
        # emerge subcluster labels into cluster labels
        adata_plot = self._subcluster_label_merge(add_key='celltype_sub')
        component_matrix, cluster_labels = calculate_neighbor_composition(
            adata=adata_plot, 
            k=k,
            celltype_key='celltype_sub')

        env_component_adata = sc.AnnData(X=component_matrix, obs=adata_plot.obs)
        env_component_adata.obs_names = adata_plot.obs_names
        env_component_adata.var_names = cluster_labels

        env_component_adata_selected = env_component_adata[self.adata_selected_plot.obs_names]
        env_component_adata_selected.obs = self.adata_selected_plot.obs.copy()

        var_show = env_component_adata_selected.var_names
        low_cnt_cluster = env_component_adata_selected.var_names[env_component_adata_selected.X.mean(axis=0) < 0.01]
        for c in low_cnt_cluster:
            var_show = np.delete(var_show, np.where(var_show == c))
        var_show = np.delete(var_show, np.where(var_show == 'Unassigned'))

        from matplotlib import rcParams
        rcParams['legend.fontsize'] = 7
        rcParams['axes.labelsize'] = 7
        rcParams['xtick.labelsize'] = 7
        rcParams['ytick.labelsize'] = 7

        sc.pl.heatmap(env_component_adata_selected, 
                    groupby='cluster_label', 
                    cmap='Oranges',
                    figsize=(4.5, 3),
                    swap_axes=True,
                    var_names = var_show,
                    # save=f"_{dataset}_{selected_cluster}_component_matrix_heatmap.pdf"
                    )

    def morphology_analysis(self, morphology_df):
        if self.adata_selected_plot is None:
            print("Run detect_subclusters() first...")
            return
        adata_selected_morphology = self.adata_selected_plot.copy()
        morphology_df.index = morphology_df.index.astype(str)

        adata_selected_morphology = adata_selected_morphology[[str(obs) in morphology_df.index for obs in adata_selected_morphology.obs_names]]
        cluster_label_list = adata_selected_morphology.obs['cluster_label'].copy()
        adata_selected_morphology.obs = morphology_df.loc[adata_selected_morphology.obs_names.astype(str)]
        adata_selected_morphology.obs['cluster_label'] = list(cluster_label_list)
        adata_selected_morphology.obs['cluster_label'] = adata_selected_morphology.obs['cluster_label'].astype("category")

        keys = [col for col in morphology_df.columns if col != 'Cell_Type']
        significant_pvals = dict()
        significant_boxplot = dict()
        for key in keys:
            group0 = adata_selected_morphology[adata_selected_morphology.obs['cluster_label'] == '0', :].obs[key]
            group1 = adata_selected_morphology[adata_selected_morphology.obs['cluster_label'] == '1', :].obs[key]


            p_value = calculate_p_value(group0, group1)
            print(f"p-value for {key}: {p_value:.2e}")

            if p_value < 5e-2:
                plot_data = pd.DataFrame({
                    'cluster_label': ['0']*len(group0) + ['1']*len(group1),
                    key: pd.concat([group0, group1])
                })
                significant_pvals[key] = p_value
                significant_boxplot[key] = plot_data

        n_significant = len(significant_pvals)
        if n_significant == 0:
            print("No significant morphology features found.")
            return
        
        n_row = int((n_significant+2)/3)
        fig, axes = plt.subplots(n_row, 3, figsize=(10, 3*n_row), dpi=self.dpi)
        for i in range(n_significant):
            key = list(significant_pvals.keys())[i]
            ax = axes.flatten()[i]
            plot_data = significant_boxplot[key]
            p_value = significant_pvals[key]
            sns.boxplot(
                x='cluster_label',
                y=key,
                data=plot_data,
                palette=self.colors_3[:self.n_subcluster], 
                width=0.45,                     
                linewidth=1.8,                  
                fliersize=3.5,                  
                notch=False,                    
                showfliers=False,                  
                ax=ax
            )
            ax.set_xlabel('Cluster Label', fontsize=7)
            ax.set_ylabel(key.replace('_', ' ').title(), fontsize=7)
            ax.set_title(f"{key}\n(p={p_value:.2e})", fontsize=7)
        plt.setp(ax.lines, linewidth=1)  # 控制须线粗细
        plt.tight_layout()
        plt.savefig(self.save_path + f"/Morphology_boxplot.pdf", bbox_inches='tight', dpi=self.dpi)
        plt.show()
        plt.close()

