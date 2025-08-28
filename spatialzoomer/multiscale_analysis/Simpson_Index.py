import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
import matplotlib.pyplot as plt


def binarize_knn_matrix(knn_matrix):
    """
    Set all non-zero elements in the KNN matrix to 1.
    
    :param knn_matrix: KNN adjacency matrix (sparse matrix)
    :return: Binarized KNN matrix
    """
    binarized_knn = knn_matrix.copy()
    binarized_knn.data = np.ones_like(binarized_knn.data)
    
    return binarized_knn



def calculate_neighbor_composition(knn_matrix, clusters, num_clusters):
    """
    Calculate the proportion of each cluster type in the neighbors of each cell.
    
    :param knn_matrix: KNN adjacency matrix (sparse matrix)
    :param clusters: Array of cluster labels for cells
    :param num_clusters: Number of clusters
    :return: Proportion of each cluster type in each cell's neighbors (num_cells, num_clusters)
    """
    num_cells = knn_matrix.shape[0]

    label_matrix = np.zeros((num_cells, num_clusters))
    label_matrix[np.arange(num_cells), clusters] = 1

    neighbor_labels_count = knn_matrix.dot(label_matrix)

    row_sums = neighbor_labels_count.sum(axis=1)  
    row_sums[row_sums == 0] = 1  
    neighbor_labels_proportion = neighbor_labels_count / row_sums[:, np.newaxis]    
    return neighbor_labels_proportion


def calculate_simpson_index_from_composition(composition):
    """
    Calculate the Simpson index from the neighbor composition proportions.
    
    :Matrix of neighbor cluster composition proportions (num_cells, num_clusters)
    :return: Simpson index for each cell (num_cells,)
    """
    return  np.sum(composition ** 2, axis=1)



def calculate_simpson_indices(spatial_knn, transcriptomic_knn, adata, cluster_col):
    """
    Calculate Simpson indices for spatial and transcriptomic KNN graphs.
    
    :param spatial_knn: Spatial KNN graph (sparse matrix)
    :param transcriptomic_knn: Transcriptomic KNN graph (sparse matrix)
    :param adata: AnnData object containing cluster annotations
    :param cluster_col: Column name of cluster labels in adata.obs
    :return: Simpson indices for spatial and transcriptomic graphs

    """
    clusters = pd.Categorical(adata.obs[cluster_col]).codes
    num_clusters = len(np.unique(clusters))

    spatial_knn_binarized = binarize_knn_matrix(spatial_knn)
    transcriptomic_knn_binarized = binarize_knn_matrix(transcriptomic_knn)

    spatial_composition = calculate_neighbor_composition(spatial_knn_binarized, clusters, num_clusters)
    transcriptomic_composition = calculate_neighbor_composition(transcriptomic_knn_binarized, clusters, num_clusters)

    simpson_spatial = calculate_simpson_index_from_composition(spatial_composition)
    simpson_transcriptomic = calculate_simpson_index_from_composition(transcriptomic_composition)

    return simpson_spatial, simpson_transcriptomic



def plot_simpson_indices(adata, scales, clusters, save_path = None):
    """
    Plot the mean Simpson indices to compare spatial and transcriptomic influences.

    :param adata: AnnData object containing cluster labels
    :param spatial_knn_binarized: Spatial KNN matrix
    :param expr_knn_binarized: Transcriptomic KNN matrix
    :param scales: List of scales
    :param clusters: List of cluster column names at each scale
    :param save_path: Path to save the plots (optional)
    """

    spatial_knn = adata.obsp['spatial_knn'].copy()
    spatial_knn = spatial_knn.tocsr()  

    expr_knn = adata.obsp['expr_knn'].copy()
    expr_knn = expr_knn.tocsr()  
    spatial_knn_binarized = binarize_knn_matrix(spatial_knn)
    expr_knn_binarized = binarize_knn_matrix(expr_knn)

    simpson_spatial_dict = {}
    simpson_transcriptomic_dict = {}

    for scale in clusters:
        simpson_spatial, simpson_transcriptomic = calculate_simpson_indices(spatial_knn_binarized, expr_knn_binarized, adata, scale)
        simpson_spatial_dict[scale] = simpson_spatial
        simpson_transcriptomic_dict[scale] = simpson_transcriptomic

    mean_simpson_spatial = [np.mean(simpson_spatial_dict[scale]) for scale in clusters]
    mean_simpson_transcriptomic = [np.mean(simpson_transcriptomic_dict[scale]) for scale in clusters]

    plt.figure(figsize=(5, 4))

    plt.plot(scales, mean_simpson_spatial, marker='o', linestyle='-', color='b', label='Spatial information')
    plt.plot(scales, mean_simpson_transcriptomic, marker='s', linestyle='--', color='g', label='Transcriptomic information')

    plt.title('Average Simpson Indices across Scales', fontsize=8)
    plt.xlabel('Scales', fontsize=7)
    plt.ylabel('Mean Simpson Index', fontsize=7)
    plt.xticks(rotation=0)
    plt.legend(loc='upper left')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path + 'Average_Simpson_Indices_across_Scales.pdf', dpi=300)
    plt.show()


    ratios_dict = {}
    for scale in clusters:
        ratios_dict[scale] = simpson_spatial_dict[scale] / simpson_transcriptomic_dict[scale]

    return simpson_spatial_dict, simpson_transcriptomic_dict, ratios_dict