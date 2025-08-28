import math
import time
import numpy as np
from collections import Counter
from matplotlib import gridspec
from sklearn.preprocessing import StandardScaler
import squidpy as sq
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.cluster import MiniBatchKMeans
import warnings
from scipy.sparse import issparse

warnings.filterwarnings("ignore")

def TwoStepClustering(adata, data, batch_size_kmeans = 10000, n_clusters_kmeans = None, scale_data = True, resolution = [0.2, 0.4, 0.6, 0.8, 1], sigma = 0, title = 'scale'+str(0.001)):
    if n_clusters_kmeans is None:
        n_clusters_kmeans = math.floor(data.shape[0]/10)
    #start_time = time.time()
    if scale_data:
        scaler = StandardScaler()
        signal = scaler.fit_transform(data)
    else:
        signal = data.copy()

    kmeans = MiniBatchKMeans(n_clusters = n_clusters_kmeans, init='k-means++', max_iter = 300, 
                             max_no_improvement = 30,
                             batch_size = batch_size_kmeans, random_state = 1)
    kmeans.fit(signal)
    kmeans_labels = kmeans.labels_
    kmeans_cluster_centers = kmeans.cluster_centers_
    adata.obs['Kmeans_'+ str(title)] = kmeans_labels
    
    cell_metadata = pd.DataFrame(index=[f'KMEANS_{i}' for i in range(kmeans_cluster_centers.shape[0])])
    gene_metadata = pd.DataFrame(index=[f'P_{i}' for i in range(kmeans_cluster_centers.shape[1])])
    adata_kmeans = sc.AnnData(X=kmeans_cluster_centers, obs=cell_metadata, var=gene_metadata)
    
    sc.pp.neighbors(adata_kmeans, use_rep = 'X')
    sc.tl.umap(adata_kmeans)

    cell_umap_coords = []
    for i, label in enumerate(kmeans_labels):
        center_coords = adata_kmeans.obsm['X_umap'][label]
        perturbed_coords = center_coords + np.random.normal(0, sigma, size=center_coords.shape)  # Adding small noise
        cell_umap_coords.append(perturbed_coords)
    cell_umap_coords = np.array(cell_umap_coords)
    adata.obsm['X_umap_'+str(title)] = cell_umap_coords

    for res in resolution:
        sc.tl.leiden(adata_kmeans, resolution = res)
        adata_kmeans.obs['leiden_'+ str(title)+ '_res' + str(res)] = adata_kmeans.obs['leiden']
        
        agg_labels = adata_kmeans.obs['leiden']
        assigned_labels = agg_labels[kmeans_labels]
        label_counts = Counter(assigned_labels)
        #assigned_labels = agg_labels[kmeans_labels]
        sorted_labels = sorted(label_counts, key=label_counts.get, reverse=True)
        label_names = {label: f'C{i}' for i, label in enumerate(sorted_labels)}
        named_labels = np.array([label_names[label] for label in assigned_labels])
        named_labels = pd.Categorical(named_labels, categories=list(label_names.values()))
        adata.obs['leiden_'+str(title)+'_res'+str(res)] = named_labels
        
    #end_time = time.time()
    #print(f"Time taken to perform clustering: {(end_time - start_time)} seconds")
    return adata, adata_kmeans




def plotClustering(adata, adata_kmeans, keys = ['leiden_scale'+str(0.01)+'_res0.6'], title = 'scale'+str(0.01), ref_keys = 'leiden_Raw_res1',
                   plot_number = 4, figsize = (22, 5.5), width_ratios=[8, 14], dpi = 300, show = False, save_path=None):
    #tmp_labels, counts = np.unique(adata.obs[keys], return_counts=True)
    #unique_labels = tmp_labels[np.argsort(-counts)]
    
    fig = plt.figure(figsize=figsize)  # 总宽度是 12+16=28，高度是 8
    
    if (plot_number == 2):
        gs = gridspec.GridSpec(1, 2, width_ratios=width_ratios)  # 12:16 的宽度比例

        # 左边绘制 UMAP 图，大小为 12x8
        ax1 = plt.subplot(gs[0, 0])
        sc.pl.embedding(
            adata_kmeans,
            basis='X_umap',
            color=keys[0],  # 使用 leinden_scale 对应的列
            ax=ax1,
            show=False  # 不立即显示图像
        )
        
        # 右上角：空间散点图，基于 leiden_scale
        ax2 = plt.subplot(gs[0, 1])
        sq.pl.spatial_scatter(
            adata,
            library_id="spatial",
            shape=None,
            color=keys[0],  # 使用 leinden_scale 对应的列
            wspace=0.4,
            ax=ax2
        )
    
    if (plot_number == 4):
        gs = gridspec.GridSpec(2, 2, width_ratios=width_ratios)  # 12:16 的宽度比例

        # 左边绘制 UMAP 图，大小为 12x8
        ax1 = plt.subplot(gs[0, 0])
        sc.pl.embedding(
            adata_kmeans,
            basis='X_umap',
            color=keys[0],  # 使用 leinden_scale 对应的列
            ax=ax1,
            show=False  # 不立即显示图像
        )
        
        # 右上角：空间散点图，基于 leiden_scale
        ax2 = plt.subplot(gs[0, 1])
        sq.pl.spatial_scatter(
            adata,
            library_id="spatial",
            shape=None,
            color=keys[0],  
            wspace=0.4,
            ax=ax2
        )

        # 左下角：UMAP 图，基于 ref_keys
        ax3 = plt.subplot(gs[1, 0])
        sc.pl.embedding(
            adata,
            basis='X_umap_' + str(title),
            color=ref_keys,
            size = 120000 /adata_kmeans.shape[0],
            ax=ax3,
            show=False  # 不立即显示图像
        )
        
        # 右下角：空间散点图，基于 ref_keys
        ax4 = plt.subplot(gs[1, 1])
        sq.pl.spatial_scatter(
            adata,
            library_id="spatial",
            shape=None,
            color=ref_keys,  
            wspace=0.4,
            ax=ax4
        )
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    if save_path:
        plt.savefig(save_path, dpi=dpi)
    
    # 关闭图像
    if show:
        plt.show()
    else:
        plt.close()




def Clustering_raw_signal(adata, 
                          use_rep = 'X_nmf', 
                          n_clusters_kmeans = 10000, 
                          resolutions = [0.4, 0.6, 0.8, 1, 1.2], 
                          scale_data = True, 
                          figsize = (22, 5.5),
                          title = 'Raw', save_path = '', runLabel = ''):
    """
    Perform clustering on the raw signal after NMF transformation
    """
    start_time = time.time()
    if use_rep == 'X':
        data = adata.X.copy()
        if issparse(data):
            print("Converting sparse matrix to dense format for 'X' representation. This may consume a lot of memory.")
            data = data.toarray()
    else:
        if use_rep not in adata.obsm:
            raise ValueError(f"Data representation '{use_rep}' not found in `adata.obsm`.")
        data = adata.obsm[use_rep].copy()

    adata, adata_kmeans =  TwoStepClustering(adata, data, n_clusters_kmeans = n_clusters_kmeans, scale_data = scale_data, resolution = resolutions, title = title)
    adata_kmeans.write(save_path + '/' + runLabel +'_object_kmeans_'+title+'.h5ad')

    # for res in resolutions:
    #     path  =  save_path +'/Clustering_signal_'+str(title) + '_res'+ str(res)+'.png'
    #     plotClustering(adata, adata_kmeans, keys = ['leiden_'+str(title)+'_res'+str(res)], title = title, save_path=path, dpi = 300, plot_number = 2, figsize = figsize)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken to perform UMAP and clustering for raw signal: {elapsed_time:.4f} seconds")
    return adata



def Clustering_transformed_signal(adata, transformed_signals, typical_scales_df, n_clusters_kmeans = 10000, 
                                  resolutions = [0.4, 0.6, 0.8, 1, 1.2], scale_data = True, figsize = (22, 11), save_path = '', runLabel = ''):
    """
    Perform clustering on the transformed signal after GWT transformation
    """
    start_time1 = time.time()
    typical_scales = [round(x, 2) for x in  typical_scales_df['Scale'].values]
    typical_scales_index = typical_scales_df['Index'].values
    typical_transformed_signals = transformed_signals[:, :, typical_scales_index]

    for i, scale in enumerate(typical_scales):
        start_time = time.time()
        data = typical_transformed_signals[:, :, i]
        title = 'scale'+str(scale)

        adata, adata_kmeans =  TwoStepClustering(adata, data, n_clusters_kmeans = n_clusters_kmeans, scale_data = scale_data,  resolution = resolutions, title = title)
        adata_kmeans.write(save_path + '/' + runLabel +'_object_kmeans_'+title+'.h5ad')

        # for res in resolutions:
        #     path  =  save_path +'/Clustering_signal_'+str(title) + '_res'+ str(res)+'.png'
        #     plotClustering(adata, adata_kmeans, keys = ['leiden_'+str(title)+'_res'+str(res)], title = title, save_path=path, dpi = 300, plot_number = 4, figsize = figsize, ref_keys = 'leiden_Raw_res1')

        end_time = time.time()
        print(f"Time taken to perform clustering for signal at scale {scale}: {end_time - start_time:.4f} seconds")

    elapsed_time = end_time - start_time1
    print(f"Total time: {elapsed_time:.4f} seconds")

    return adata