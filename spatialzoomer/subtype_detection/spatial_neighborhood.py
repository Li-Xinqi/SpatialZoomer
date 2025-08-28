from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
import numpy as np

def getKNN(adata, k):
    coords = adata.obsm['spatial'].copy()
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='auto', metric='euclidean').fit(coords)
    distances, indices = nbrs.kneighbors(coords)
    
    distances = distances[:, 1:]
    indices = indices[:, 1:]

    epsilon = 1e-10  
    similarities  = 100 / (distances + epsilon)

    n_cells = coords.shape[0]
    row_indices = np.repeat(np.arange(coords.shape[0]), k)
    col_indices = indices.flatten()
    data = similarities.flatten()
    
    knn_similarity_matrix = csr_matrix((data, (row_indices, col_indices)), shape=(n_cells, n_cells))
    knn_similarity_matrix = knn_similarity_matrix.maximum(knn_similarity_matrix.T)
    
    return knn_similarity_matrix

def calculate_neighbor_composition_fromKNN(knn_matrix, clusters, num_clusters):
    num_cells = knn_matrix.shape[0]
    label_matrix = np.zeros((num_cells, num_clusters))
    label_matrix[np.arange(num_cells), clusters] = 1

    neighbor_labels_count = knn_matrix.dot(label_matrix)
    
    row_sums = neighbor_labels_count.sum(axis=1)
    row_sums[row_sums == 0] = 1  # avoid / 0
    neighbor_labels_proportion = neighbor_labels_count / row_sums[:, np.newaxis]
    return neighbor_labels_proportion

def calculate_neighbor_composition(adata, k, celltype_key):
    celltype_labels = adata.obs[celltype_key]
    n_cluster = len(np.unique(celltype_labels))
    knn_similarity_matrix = getKNN(adata, k=k)
    
    # 建立标签到序号的映射
    cluster_labels = np.unique(celltype_labels)
    cluster_label_to_index = {label: i for i, label in enumerate(cluster_labels)}
    # 将标签转换为序号
    clusters = np.array([cluster_label_to_index[label] for label in celltype_labels])
    # 计算邻居组成比例
    component_matrix = calculate_neighbor_composition_fromKNN(
        knn_similarity_matrix, 
        clusters=clusters, 
        num_clusters=n_cluster)
    return component_matrix, cluster_labels