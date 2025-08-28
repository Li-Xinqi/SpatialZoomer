import numpy as np
import matplotlib.pyplot as plt
from kneed import KneeLocator
from scipy.spatial.distance import cdist
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
import seaborn as sns  
import pandas as pd
import time


def compute_within_variance(corr_matrix, start, end):
    """
    Compute the variance within a sub-matrix of the correlation matrix.

    Parameters:
    -----------
    corr_matrix : ndarray
        The full correlation matrix.
    
    start : int
        The starting index of the sub-matrix.
    
    end : int
        The ending index of the sub-matrix.

    Returns:
    --------
    float
        The variance of the sub-matrix from start to end.
    """
    sub_matrix = corr_matrix[start:end+1, start:end+1]
    return np.var(sub_matrix)



def dp_optimal_clustering_total_variance(corr_matrix, n_clusters):
    """
    Use dynamic programming to find the optimal clustering that minimizes the total within-cluster variance.

    Parameters:
    -----------
    corr_matrix : ndarray
        The correlation matrix where clustering is to be performed.
    
    n_clusters : int
        The desired number of clusters.

    Returns:
    --------
    float
        The total within-cluster variance for the optimal clustering.
    """

    n = corr_matrix.shape[0]
    dp = np.full((n + 1, n_clusters + 1), np.inf)
    dp[0][0] = 0
    div_points = np.zeros((n + 1, n_clusters + 1), dtype=int)

    for i in range(1, n+1):
        dp[i][1] = compute_within_variance(corr_matrix, 0, i-1)

    for k in range(2, n_clusters + 1):
        for i in range(1, n + 1):
            for j in range(1, i):
                cost = dp[j][k-1] + compute_within_variance(corr_matrix, j, i-1)
                if cost < dp[i][k]:
                    dp[i][k] = cost
                    div_points[i][k] = j

    return dp[n][n_clusters] 



def dp_optimal_clustering(corr_matrix, n_clusters):
    """
    Perform dynamic programming to find the optimal partitioning of the correlation matrix 
    into the specified number of clusters.

    Parameters:
    -----------
    corr_matrix : ndarray
        The correlation matrix where clustering is to be performed.
    
    n_clusters : int
        The desired number of clusters.

    Returns:
    --------
    list of tuple
        A list of tuples indicating the start and end indices of each cluster.
    """

    n = corr_matrix.shape[0]
    dp = np.full((n + 1, n_clusters + 1), np.inf)  # Dynamic programming array
    dp[0][0] = 0
    div_points = np.zeros((n + 1, n_clusters + 1), dtype=int)  # Store division points
    
    
    for i in range(1, n+1):
        dp[i][1] = compute_within_variance(corr_matrix, 0, i-1)

    # Recursive computation of optimal division
    for k in range(2, n_clusters + 1):  # Start from 2 clusters
        for i in range(1, n + 1):
            for j in range(1, i):
                cost = dp[j][k-1] + compute_within_variance(corr_matrix, j, i-1)
                if cost < dp[i][k]:
                    dp[i][k] = cost
                    div_points[i][k] = j  # Record optimal division point

    # Step 3: Trace back division points to get optimal partitions
    clusters = []
    end = n
    for k in range(n_clusters, 0, -1):
        start = div_points[end][k]
        clusters.append((start, end - 1))
        end = start

    return clusters[::-1]  # Return sorted partitions



def find_cluster_centers(corr_matrix, clusters):
    """
    Identify the cluster centers for each cluster. The center is defined as the point with 
    the minimum sum of distances to all other points in the cluster.

    Parameters:
    -----------
    corr_matrix : ndarray
        The correlation matrix where clustering was performed.
    
    clusters : list of tuple
        A list of tuples indicating the start and end indices of each cluster.

    Returns:
    --------
    list of int
        A list of indices representing the center point of each cluster.
    """

    centers = []
    for start, end in clusters:
        sub_matrix = corr_matrix[start:end+1, start:end+1]
        distances = np.sum(cdist(sub_matrix, sub_matrix, metric='euclidean'), axis=1)
        center_idx = np.argmin(distances)
        centers.append(start + center_idx)
    return centers



def Identify_Typical_Scales(transformed_signals, scales, min_clusters = 3, max_clusters = 10, figsize = (12, 6), save_path=None):
    """
    Identify the optimal number of frequency bands using dynamic programming and the elbow method 
    for clustering based on within-cluster variance.

    Parameters:
    -----------
    transformed_signals : ndarray
        Transformed signals where clustering is performed (typically the result of NMF and graph wavelet transform).
    
    scales : list of float
        The scales corresponding to the signal transformations.
    
    min_clusters : int, optional, default: 3
        Minimum number of clusters to consider.

    max_clusters : int, optional, default: 10
        Maximum number of clusters to consider for the elbow method.
    
    figsize : tuple, optional, default: (5, 10)
        Figure size for the elbow plot and heatmap.

    save_path : str, optional, default: None
        Path to save plots. If None, the plots will not be saved.

    Returns:
    --------
    mean_correlation_matrix : ndarray
        The mean correlation matrix across NMF programs.
    
    optimal_clusters : list of tuple
        The start and end indices of the optimal clusters.
    
    scale_centers : list of float
        The center scale of each cluster.
    """
    start_time = time.time()
    # Step 1: Compute the correlation matrix across scales for each NMF program
    correlation_matrices = []
    for nmf_idx in range(transformed_signals.shape[1]):  
        nmf_scales = transformed_signals[:, nmf_idx, :]        
        corr_matrix = np.corrcoef(nmf_scales.T)  
        correlation_matrices.append(corr_matrix)

    correlation_matrices = np.stack(correlation_matrices)
    mean_correlation_matrix = np.mean(correlation_matrices, axis=0)

    # Step 2: Compute total within-cluster variance for a range of cluster numbers
    variances = []
    for k in range(min_clusters, max_clusters + 1):
        total_variance = dp_optimal_clustering_total_variance(mean_correlation_matrix, k)
        variances.append(total_variance)
    
    # Step 3: Detect the elbow point
    kneedle = KneeLocator(range(3, max_clusters + 1), variances, curve='convex', direction='decreasing',  S=1)

    # fig, axs = plt.subplots(1, 2, figsize=figsize)
    plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 6])
    ax1 = plt.subplot(gs[0])
    ax1.plot(range(3, max_clusters + 1), variances, marker='o')
    ax1.set_xlabel('Number of Clusters')
    ax1.set_ylabel('Total Within-Cluster Variance')
    ax1.set_title('Elbow Method for Optimal Number of Clusters')
    optimal_k = kneedle.elbow
    if optimal_k is not None:
        ax1.axvline(x=optimal_k, color='red', linestyle='--', label=f'Elbow at k={optimal_k}')
        ax1.legend()

    ax2 = plt.subplot(gs[1])
    sns.heatmap(mean_correlation_matrix, cmap='coolwarm', annot=False, cbar=True, 
                xticklabels=scales, yticklabels=scales, ax=ax2)
    ax2.set_title('Optimal Clustering Using Dynamic Programming')

    optimal_clusters = dp_optimal_clustering(mean_correlation_matrix, optimal_k)
    cluster_centers = find_cluster_centers(mean_correlation_matrix, optimal_clusters)
    for i, (start, end) in enumerate(optimal_clusters):
        ax2.add_patch(Rectangle((start, start), end - start + 1, end - start + 1, 
                                fill=False, edgecolor='black', lw=2))
        ax2.scatter(start + 0.5, start + 0.5, color='blue', marker='o', s=50)
    for center in cluster_centers:
        ax2.scatter(center + 0.5, center + 0.5, color='blue', marker='o', s=50)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path + '/typical_scales.png', dpi=300)
        # plt.savefig(f'{save_path}_combined.pdf', dpi=300)
    plt.show()

    print(f"The optimal number of clusters is: {optimal_k}")
    
    # print("Optimal clusters (start, end):", optimal_clusters)
    scale_centers = [scales[scale] for scale in cluster_centers]
    # print("Cluster centers:", scale_centers)

    start_index = [optimal_clusters[i][0] for i in range(len(optimal_clusters))]
    typical_scales = scale_centers + [scales[i] for i in start_index]
    typical_scales = sorted(typical_scales)
    print("Typical scales:", typical_scales)

    typical_scales_indices = [scales.index(scale) for scale in typical_scales]
    
    typical_scales_df = pd.DataFrame({
    'Index': typical_scales_indices,
    'Scale': typical_scales
    })
    if save_path:
        typical_scales_df.to_csv(f'{save_path}_typical_scales.csv', index=False)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken to perform identification of typical scales: {elapsed_time:.4f} seconds")

    return mean_correlation_matrix, optimal_clusters, typical_scales_df