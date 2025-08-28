from .Preprocessing import Preprocess
from .Spectral_graph_filter_bank import performDR, getKNN, performSGFB
from .Identify_typical_scales import Identify_Typical_Scales
from .Two_step_Clustering import Clustering_raw_signal, Clustering_transformed_signal
from .Simpson_Index import plot_simpson_indices
from .Visualization import plot_clusters

from matplotlib import rcParams
import matplotlib as mpl 
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import squidpy as sq
import os

class MultiscaleAnalysis:
    def __init__(self, adata, runLabel, save_path):
        self.adata = adata
        self.runLabel = runLabel
        self.save_path = save_path
        # if not exist
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)
        
        # multiscale transform
        self.use_rep = None
        self.transformed_signals = None
        self.scales = None
        self.n_neighbors = None

        # typical scales selection
        self.typical_scales_df = None

        # clustering
        self.scales_df_use = None
        self.clustered = False
        self.resolutions = None

        # simpson indices
        self.simpson_spatial_dict = None
        self.simpson_expr_dict = None
        self.ratios_dict = None

        self._set_visual_config()

    def _set_visual_config(self):
        self.dpi=200
        config = {
            "font.family":'sans-serif',
            "font.size": 6,
            "mathtext.fontset":'stix',
            "font.serif": ['MS Arial'],
        }
        rcParams.update(config)

        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42        

    def multiscale_transform(self, use_rep='X_nmf', scales=None, n_neighbors=10):
        self.n_neighbors = n_neighbors
        if scales is None:
            self.scales = [0.01] + np.arange(0.1, 2.1, 0.1).tolist() + np.arange(2.5, 15.5, 0.5).tolist() + np.arange(16, 21, 1).tolist() +  np.arange(25, 55, 5).tolist()
            self.scales = [round(x, 2) for x in self.scales]
        else:
            self.scales = scales

        self.use_rep = use_rep
        if use_rep not in self.adata.obsm:
            raise ValueError(f"Data representation '{use_rep}' not found in `adata.obsm`.")

        # adata = performDR(adata, type = 'NMF', n_components=50)
        self.adata = getKNN(
            self.adata, 
            use_rep = 'spatial', 
            name = 'spatial_knn', 
            n_neighbors = n_neighbors, 
            pattern='auto')

        self.transformed_signals, self.scales = performSGFB(
            self.adata, 
            use_knn = 'spatial_knn', 
            use_rep = use_rep, 
            scales = self.scales, 
            runLabel = self.runLabel, 
            save_path = self.save_path)

    def identify_typical_scales(self, max_clusters=10, min_clusters=3):
        if self.transformed_signals is None:
            raise ValueError("Transformed signals are not available. Please run `multiscale_transform` first.")
        
        _, _, self.typical_scales_df = Identify_Typical_Scales(
            self.transformed_signals, 
            self.scales, 
            max_clusters=max_clusters, 
            min_clusters=min_clusters,
            save_path = self.save_path
            )
        
    def clustering(self, n_clusters_kmeans=10000, resolutions=None, min_scale=None, max_scale=None):
        """
        Perform clustering on both raw and transformed signals.
        """
        if resolutions is None:
            resolutions = [0.4, 0.6, 0.8, 1, 1.2]
        self.resolutions = resolutions
        
        if self.typical_scales_df is None:
            raise ValueError("Typical scales are not available. Please run `identify_typical_scales` first.")

        self.adata = Clustering_raw_signal(
            self.adata, 
            use_rep = self.use_rep,  
            n_clusters_kmeans = n_clusters_kmeans, 
            resolutions = resolutions, 
            title = 'Raw', 
            save_path = self.save_path, 
            runLabel = self.runLabel)
        
        scales_df_use = self.typical_scales_df.copy()
        if min_scale is not None:
            scales_df_use = scales_df_use[scales_df_use['Scale'] >= min_scale]
        if max_scale is not None:
            scales_df_use = scales_df_use[scales_df_use['Scale'] <= max_scale]

        self.adata = Clustering_transformed_signal(
            self.adata, 
            self.transformed_signals, 
            scales_df_use, 
            n_clusters_kmeans = n_clusters_kmeans, 
            resolutions = resolutions, 
            save_path = self.save_path, 
            runLabel = self.runLabel
            )
        self.scales_df_use = scales_df_use
        self.clustered = True
    
    def plot_simpson(self, resolution=None):
        if self.clustered is False:
            raise ValueError("Clustering has not been performed. Please run `clustering` first.")
        if resolution is None:
            resolution = self.resolutions[-1]
        
        self.adata = getKNN(
            self.adata, 
            use_rep = self.use_rep, 
            name = 'expr_knn', 
            n_neighbors = self.n_neighbors
            )
        typical_scales = self.scales_df_use['Scale'].values
        clusters_use = ['leiden_Raw_res' + str(resolution)] + ['leiden_scale' + str(scale) + '_res' + str(resolution) for scale in typical_scales]
        scales_plot = ['Raw'] + [str(scale) for scale in typical_scales]
        simpson_spatial_dict, simpson_expr_dict, ratios_dict = plot_simpson_indices(
            self.adata, 
            scales_plot, 
            clusters_use, 
            save_path = self.save_path)
        self.simpson_spatial_dict = simpson_spatial_dict
        self.simpson_expr_dict = simpson_expr_dict
        self.ratios_dict = ratios_dict
        # return simpson_spatial_dict, simpson_expr_dict, ratios_dict
    
    def perform_multiscale_pipeline(
            self, 
            use_rep, 
            n_clusters_kmeans=10000, 
            n_neighbors=10, 
            scales=None, 
            resolutions=None, 
            max_clusters=10, 
            min_clusters=3,
            min_scale=None,
            max_scale=None
            ):
        """
        Perform the entire multiscale analysis pipeline."
        """
        if self.transformed_signals is None:
            self.multiscale_transform(
                use_rep=use_rep, 
                scales=scales, 
                n_neighbors=n_neighbors
                )
        if self.typical_scales_df is None:
            self.identify_typical_scales(
                max_clusters=max_clusters, 
                min_clusters=min_clusters
                )
        if self.clustered is False:
            self.clustering(
                n_clusters_kmeans=n_clusters_kmeans, 
                resolutions=resolutions,
                min_scale=min_scale,
                max_scale=max_scale
                )
        self.plot_simpson()

    def plot_multiscale_clusters(self, resolution=0.4, pt_size_umap = 5, pt_size_scatter = 0.2, figsize = None, width_ratios=None):
        if self.clustered is False:
            raise ValueError("Clustering has not been performed. Please run `clustering` first.")
        
        typical_scales = self.scales_df_use['Scale'].values
        plot_clusters(
            self.adata, 
            typical_scales, 
            resolution=resolution,
            figsize=figsize,
            width_ratios=width_ratios,
            pt_size_umap= pt_size_umap,
            pt_size_scatter= pt_size_scatter
            )