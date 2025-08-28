import numpy as np
from scipy import stats

def calculate_p_value(group1, group2):
    s_stat, p_value = stats.ttest_ind(group1, group2, equal_var=False)
    return p_value

def calculate_centroids(X, labels):
    unique_labels = np.unique(labels)
    centroids = np.vstack([X[labels == label].mean(axis=0) for label in unique_labels])
    return centroids, unique_labels

def ward_linkage(X, y):
    centroids, unique_labels = calculate_centroids(X, y)
    sizes = np.array([np.sum(y == label) for label in unique_labels])

    all_centroids_exp = centroids[:, np.newaxis, :] - centroids[np.newaxis, :, :]
    centroid_dists_sq = np.sum(all_centroids_exp ** 2, axis=-1)
    ess_diff = (sizes[:, np.newaxis] * sizes[np.newaxis, :]) / (sizes[:, np.newaxis] + sizes[np.newaxis, :]) * centroid_dists_sq
    
    # 确保只考虑上三角部分的有效合并差值
    mask = np.triu(np.ones_like(ess_diff, dtype=bool), k=1)
    ess_diff[~mask] = np.inf

    min_ward_value = np.min(ess_diff)
    return min_ward_value# / len(y)
