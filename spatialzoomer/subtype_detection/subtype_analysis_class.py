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
import os

class SubclusterAnalysis:
    def __init__(self, ad, save_path, cluster_key='leiden_Raw_res1'):
        self.ad = ad
        self.umap_keys=self._get_umap_keys()
        self.cluster_key = cluster_key
        self.ad.obs[self.cluster_key] = self.ad.obs[self.cluster_key].astype('category')
        self.ad.obs_names = self.ad.obs_names.astype(str)

        self.n_subcluster = None
        self.selected_cluster = None
        self.ad_selected = None
        self.ad_selected_plot = None
        self.optimal_scale = None
        self.var_dict = None
        self.k_dict = None

        self.save_path = save_path
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        self._set_visual_config()
    
    def _get_umap_keys(self):
        obsm_key = list(self.ad.obsm.keys())
        umap_keys = [key for key in obsm_key if key.startswith('X_umap_scale')]
        umap_keys = sorted(umap_keys, key=lambda x: float(x.split('scale')[-1]))
        umap_keys = [umap_key for umap_key in umap_keys if float(umap_key.split('scale')[-1]) <= 10.0]
        return umap_keys
    
    def _set_visual_config(self):
        self.dpi=100
        self.new_colors = generate_colors(200)
        self.colors_3 = ['#CD69C9', '#9ACD32', "#FF7F00"]
        config = {
            "font.family":'sans-serif',
            "font.size": 12,
            "mathtext.fontset":'stix',
            "font.serif": ['MS Arial'],
        }
        rcParams.update(config)
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
    
    def calculate_metrics(self, nonepi=True, plot_clusters=True):
        """
        Plot covariance-scale curves and covariance values.
        input ad: AnnData object, with umaps at raw ("X_umap_Raw") and different scales (for example "X_umap_scale0.01")
        input cluster_key: key of the clustering result in ad.obs
        input save_path: path to save the figures
        """

        ######################### plot umap #########################
        # set as categorical
        cell_type_counts = self.ad.obs[self.cluster_key].value_counts()
        sorted_cell_types = cell_type_counts.index.tolist()
        color_mapping = {cell_type: self.new_colors[i] for i, cell_type in enumerate(sorted_cell_types)}
        sorted_colors = [color_mapping[cell_type] for cell_type in self.ad.obs[self.cluster_key].cat.categories]
        self.ad.uns[self.cluster_key + '_colors'] = sorted_colors
        self.ad.obsm['X_umap'] = self.ad.obsm['X_umap_Raw']

        if plot_clusters:
            fig, axs = plt.subplots(1, 2, figsize=(10, 4), dpi=self.dpi)
            sc.pl.umap(self.ad, 
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
                self.ad, 
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
            plt.savefig(self.save_path + "/spatial_clusters.pdf", bbox_inches='tight', dpi=self.dpi)
            plt.show()
            plt.close()

        
        ######################### plot covariance-scale curves #########################        
        # normalize umap
        self.ad.obsm['scaled_X_umap_Raw'] = (self.ad.obsm['X_umap_Raw'] - self.ad.obsm['X_umap_Raw'].min()) / (self.ad.obsm['X_umap_Raw'].max() - self.ad.obsm['X_umap_Raw'].min())
        for umap_scale in self.umap_keys:
            # X.obsm['umap_scale']归一化
            self.ad.obsm['scaled'+umap_scale] = (self.ad.obsm[umap_scale] - self.ad.obsm[umap_scale].min()) / (self.ad.obsm[umap_scale].max() - self.ad.obsm[umap_scale].min())
        
        # loop over clusters and scales to calculate variance
        nmf_clusters = self.ad.obs[self.cluster_key].unique()
        var_dict = {}
        k_dict = {}
        for nmf_cluster in nmf_clusters:
            var_list = []
            for umap_scale in self.umap_keys:
                X = self.ad[self.ad.obs[self.cluster_key] == nmf_cluster].obsm['scaled'+umap_scale]
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


        fig, axs = plt.subplots(1, 2, figsize=(9, 3), dpi=self.dpi)  # 调整画布大小以适应布局

        # plot variance-scale curves for each cluster on the left
        for nmf_cluster, var_list in var_dict.items():
            axs[0].plot(var_list, label=nmf_cluster, color=color_mapping[nmf_cluster])
        # 移除图例和数值
        axs[0].set_xticks(range(len(self.umap_keys)))
        axs[0].set_xticklabels([x.split('scale')[-1] for x in self.umap_keys])  # x轴标签

        # bar plot of k values for each cluster on the right
        right_ax = axs[1]
        right_ax.barh(clusters, values, color=[color_mapping[cluster] for cluster in clusters])

        right_ax.set_yticks(range(len(clusters)))
        right_ax.set_yticklabels(clusters)
        right_ax.set_title("k values")


        plt.tight_layout()
        plt.savefig(self.save_path + "/var_plot_withk.pdf", dpi=self.dpi)
        plt.show()
        plt.close()
        
        print("Suggested cluster: ", clusters[-1])    # the cluster with the largest k value
        self.var_dict = var_dict
        self.k_dict = k_dict

    def kde_plot(self):
        """
        Plot KDE plots for each cluster.
        input ad: AnnData object, with umaps at raw ("X_umap_Raw") and different scales (for example "X_umap_scale0.01")
        input cluster_key: key of the clustering result in ad.obs
        input save_path: path to save the figures
        """
        if self.selected_cluster is None:
            print("Select a cluster first. SubclusterAnalysis.selected_cluster = 'cluster_name'")
            return
        else:
            print("Selected cluster:", self.selected_cluster)
        
        if self.ad_selected is None:
            ad_selected = self.ad[self.ad.obs[self.cluster_key] == self.selected_cluster]
            self.ad_selected = ad_selected

        fig, axes = plt.subplots(int((len(self.umap_keys)+2)/3), 3, figsize=(9, len(self.umap_keys)), dpi=self.dpi)

        for key in self.umap_keys:
            # subsample for large datasets during kde plot
            if ad_selected.n_obs > 10000:
                ad_selected_subsampled = sc.pp.subsample(ad_selected, n_obs=10000, random_state=0, copy=True)
            else:
                ad_selected_subsampled = ad_selected

            X = ad_selected_subsampled.obsm[key]
            ax = axes.flatten()[self.umap_keys.index(key)]
            ax.scatter(X[:, 0], X[:, 1], color='gray', s=1, marker='.', alpha=0.5)
            for collection in ax.collections:
                collection.set_rasterized(True)  # 将散点设置为栅格化
            # KDE plot
            umap_df = pd.DataFrame(X, columns=['UMAP1', 'UMAP2'])
            sns.kdeplot(data=umap_df, x='UMAP1', y='UMAP2', ax=ax, bw_adjust=1.5,
                        cmap='Blues', fill=True, thresh=0.01, levels=8, alpha=0.7)
            ax.set_title(key)
        plt.tight_layout()
        plt.savefig(self.save_path + "/umaps_kde.pdf", dpi=self.dpi)
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
        if self.ad_selected is None:
            ad_selected = self.ad[self.ad.obs[self.cluster_key] == self.selected_cluster]
        else:
            ad_selected = self.ad_selected

        self.n_subcluster = n_clusters
        ad_selected.obsm['X_umap'] = ad_selected.obsm['X_umap_scale' + str(self.optimal_scale)]
        clustering_model = CFSFDP(n_clusters=n_clusters)
        res = clustering_model.fit_predict(ad_selected.obsm['X_umap'])
        ad_selected.obs['cluster_label'] = list(res.astype(str))

        ad_selected_plot = ad_selected[ad_selected.obs['cluster_label'] != '-1']
        ad_selected_plot.uns['cluster_label_colors'] = self.colors_3[0:2]
        fig, ax = plt.subplots(figsize=(5, 5), dpi=self.dpi)
        sc.pl.umap(ad_selected_plot, color=['cluster_label'], title='subclusters at scale' + str(self.optimal_scale), 
                show=False, ax=ax, size=5)
        plt.savefig(self.save_path + "/umap_subclusters_scale" + str(self.optimal_scale) + ".pdf", dpi=self.dpi)
        plt.show()
        plt.close()

        umap_keys = self.umap_keys
        if calculate_ward:
            ward_list = dict()
            for key in umap_keys:
                X = ad_selected_plot.obsm['scaled'+key]
                y = ad_selected_plot.obs['cluster_label']
                ward_list[key] = ward_linkage(X, y)

        if plot_all_scales:
            fig, axes = plt.subplots(int((len(umap_keys)+2)/3), 3, figsize=(9, len(umap_keys)), dpi=self.dpi)

            for key in umap_keys:
                ad_selected_plot.obsm['X_umap'] = ad_selected_plot.obsm[key]
                ax = axes.flatten()[umap_keys.index(key)]
                sc.pl.umap(ad_selected_plot, 
                        size=5,
                        color=['cluster_label'], title=key, 
                        ax=ax, 
                        show=False)
                for collection in ax.collections:
                    collection.set_rasterized(True)
                if calculate_ward:
                    ax.set_title(key + '\n ward: %.2f' % ward_list[key])
                else:
                    ax.set_title(key)
            
            plt.tight_layout()
            plt.savefig(self.save_path + "/umaps_cluster_scale" + str(self.optimal_scale) + ".pdf", dpi=self.dpi)
            plt.show()
            plt.close()
        
        if plot_spatial:
            fig, ax = plt.subplots(figsize=(6, 4), dpi=self.dpi)
            sq.pl.spatial_scatter(
                ad_selected_plot, 
                library_id="spatial", 
                shape=None, 
                color= ['cluster_label'],
                wspace=0.1,
                frameon=False,
                size=1,
                figsize=(6, 4),
                dpi=self.dpi,
                outline=False,
                img=False,
                marker='.',
                ax=ax,
                title='cluster_spatial_' + self.selected_cluster + '_scale' + str(self.optimal_scale),
            )
            for collection in ax.collections:
                collection.set_rasterized(True)
            plt.tight_layout()
            plt.savefig(self.save_path + "/spatial_" + str(self.selected_cluster) + "_scale" + str(self.optimal_scale) + ".pdf", dpi=self.dpi)
            plt.show()
            plt.close()

        self.ad_selected = ad_selected
        self.ad_selected_plot = ad_selected_plot    # without outliers

    def _deg2df_noise_filtered(self, deg_groups, subcluster_name, deg_current_cluster_df):
        result_df_subcluster = pd.DataFrame({key: deg_groups[key][subcluster_name] for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
        lfc_seleccted_cluster = deg_current_cluster_df.loc[list(result_df_subcluster['names'])]['logfoldchanges']
        result_df_subcluster['logfoldchanges_celltype'] = list(lfc_seleccted_cluster)
        print("saving DEG results for subcluster", subcluster_name)
        result_df_subcluster.to_csv(self.save_path + f"/diff_analysis_"+str(subcluster_name)+".csv")
        result_df_subcluster = result_df_subcluster[result_df_subcluster['logfoldchanges_celltype'] > 0]
        return result_df_subcluster

    def deg_analysis_subclusters(self, dotplot=True):
        if self.ad_selected_plot is None:
            print("Run detect_subclusters() first...")
            return
        sc.tl.rank_genes_groups(self.ad_selected_plot, 'cluster_label', method='wilcoxon')

        sc.tl.rank_genes_groups(self.ad, self.cluster_key, method='wilcoxon')
        result = self.ad.uns['rank_genes_groups']
        deg_current_cluster_df = pd.DataFrame({key: result[key][self.selected_cluster] for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
        deg_current_cluster_df.index = deg_current_cluster_df['names']
        deg_df_list = []
        for subcluster_name in np.unique(self.ad_selected_plot.obs['cluster_label']):
            result_df_subcluster = self._deg2df_noise_filtered(self.ad_selected_plot.uns['rank_genes_groups'], subcluster_name, deg_current_cluster_df)
            deg_df_list.append(result_df_subcluster)

        if dotplot:
            genes = []
            for deg_df in deg_df_list:
                genes = genes + list(deg_df.iloc[:6,0])
            # Generate the stacked violin plot
            green_cmap = ListedColormap(plt.cm.Greens(np.linspace(0.2, 1, 256)))
            sc.pl.dotplot(
                self.ad_selected_plot,
                var_names=genes,
                groupby='cluster_label',
                #inner='box',  # Add box plots within the violins
                figsize=(5, 2),  # Larger size for better visibility
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
            plt.savefig(self.save_path + "/dotplot_degenes.pdf", dpi=self.dpi)
            plt.show()
            plt.close()
    
    def _subcluster_label_merge(self, add_key):
        celltype_labels = self.ad.obs[self.cluster_key]
        self.ad.obs[add_key] = [celltype_labels[i] if i not in self.ad_selected.obs_names else (celltype_labels[i] + '_cluster' + self.ad_selected.obs['cluster_label'][i]) for i in celltype_labels.index]
        print("Saving cell type labels with subcluster labels...")
        pd.DataFrame({'cell_id': self.ad.obs['cell_id'], 'group': self.ad.obs[add_key]}).to_csv(self.save_path + "/cell_groups.csv", index=False)
        ad_plot = self.ad[self.ad.obs[add_key] != self.selected_cluster+'_cluster-1']
        return ad_plot
    
    def neighborhood_analysis(self, k=10):
        if self.ad_selected_plot is None:
            print("Run detect_subclusters() first...")
            self.detect_subclusters()
            return
        
        # emerge subcluster labels into cluster labels
        ad_plot = self._subcluster_label_merge(add_key='celltype_sub')
        component_matrix, cluster_labels = calculate_neighbor_composition(
            adata=ad_plot, 
            k=k,
            celltype_key='celltype_sub')

        env_component_ad = sc.AnnData(X=component_matrix, obs=ad_plot.obs)
        env_component_ad.obs_names = ad_plot.obs_names
        env_component_ad.var_names = cluster_labels

        env_component_ad_selected = env_component_ad[self.ad_selected_plot.obs_names]
        env_component_ad_selected.obs = self.ad_selected_plot.obs.copy()

        var_show = env_component_ad_selected.var_names
        low_cnt_cluster = env_component_ad_selected.var_names[env_component_ad_selected.X.mean(axis=0) < 0.01]
        for c in low_cnt_cluster:
            var_show = np.delete(var_show, np.where(var_show == c))
        var_show = np.delete(var_show, np.where(var_show == 'Unassigned'))
        sc.pl.heatmap(env_component_ad_selected, 
                    groupby='cluster_label', 
                    cmap='Oranges',
                    figsize=(4, 2.5),
                    swap_axes=True,
                    var_names = var_show,
                    # save=f"_{dataset}_{selected_cluster}_component_matrix_heatmap.pdf"
                    )

    def morphology_analysis(self, morphology_df):
        if self.ad_selected_plot is None:
            print("Run detect_subclusters() first...")
            return
        ad_selected_morphology = self.ad_selected_plot.copy()
        morphology_df.index = morphology_df.index.astype(str)

        ad_selected_morphology = ad_selected_morphology[[str(obs) in morphology_df.index for obs in ad_selected_morphology.obs_names]]
        cluster_label_list = ad_selected_morphology.obs['cluster_label'].copy()
        ad_selected_morphology.obs = morphology_df.loc[ad_selected_morphology.obs_names.astype(str)]
        ad_selected_morphology.obs['cluster_label'] = list(cluster_label_list)
        ad_selected_morphology.obs['cluster_label'] = ad_selected_morphology.obs['cluster_label'].astype("category")

        keys = [col for col in morphology_df.columns if col != 'Cell_Type']
        significant_pvals = dict()
        significant_boxplot = dict()
        for key in keys:
            group0 = ad_selected_morphology[ad_selected_morphology.obs['cluster_label'] == '0', :].obs[key]
            group1 = ad_selected_morphology[ad_selected_morphology.obs['cluster_label'] == '1', :].obs[key]
            

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
                palette=self.colors_3[:self.n_subcluster],  # 继承原始颜色方案
                width=0.45,                     # 控制箱体宽度
                linewidth=1.8,                  # 箱线轮廓粗细
                fliersize=3.5,                  # 异常点大小
                notch=False,                    # 是否显示中位数的凹口
                showfliers=False,                  # 是否显示异常值
                ax=ax
            )
            ax.set_xlabel('Cluster Label', fontsize=12)
            ax.set_ylabel(key.replace('_', ' ').title(), fontsize=10)
            ax.set_title(f"{key}\n(p={p_value:.2e})", fontsize=11)
        plt.setp(ax.lines, linewidth=1)  # 控制须线粗细
        plt.tight_layout()
        plt.savefig(self.save_path + f"/morphology_boxplot_{key}.pdf", bbox_inches='tight', dpi=self.dpi)
        plt.show()
        plt.close()

