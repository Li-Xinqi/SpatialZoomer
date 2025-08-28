import numpy as np
from scipy.spatial.distance import pdist, squareform


class CFSFDP(object):
    """Cluster data by algorithm proposed in ``Clustering by fast search and find of density peaks''.
    Parameters
    ----------
    distance_cutoff: float, optional
        cutoff distance
    center_threshold: float, optional
        if delta_i * rho_i > center_thresold, data_i is regarded as a cluster center
    halos_threshold: float, optional
        if density rho_i < center_thresold, data_i is regarded as a halos

    Attributes
    ----------
    cluster_centers_ : array, [n_clusters, n_features]
        Coordinates of cluster centers.
    labels_ : array, [n_samples]
        Labels of each point. -1 for halos
    """
    def __init__(self, distance_cutoff=None, n_clusters=None, halos_threshold=0):
        self.distance_cutoff = distance_cutoff
        self.n_cluster = n_clusters
        self.halos_threshold = halos_threshold

        self.cluster_centers_ = None
        self.labels_ = None

    def get_distance_cutoff(self, dist):
        """Initialize parameter `distance_cutoff` from the rule described in the paper:
        Choose distance_cutoff so that the average number of neighbors is around 1% of the n_samples.

        Parameters
        ----------
        dist : array-like matrix, shape=(n_samples, n_samples)
            distance matrix
            Index of the cluster each sample belongs to.
        """
        if self.distance_cutoff is not None:
            return self.distance_cutoff

        kth = max(dist.shape[0] // 100, 1)
        partition = np.partition(dist, kth, axis=1)
        return partition[:, kth].mean()

    def get_n_center(self, confidence):
        if self.n_cluster is not None:
            return self.n_cluster

        # confidence_sorted = np.sort(confidence[confidence > 1e-3])[::-1]
        confidence_sorted = np.sort(confidence)[::-1]
        return (confidence_sorted[1:-1] / confidence_sorted[2:])[:confidence.shape[0] // 10].argmax() + 2

    def fit(self, X):
        """Compute clustering."""
        np.random.seed(0)
        n_samples = X.shape[0]

        dist = pdist(X, metric='euclidean')
        dist_square = squareform(dist)
        dist_square /= dist_square.max()

        distance_cutoff = self.get_distance_cutoff(dist_square)
        self.rho = (dist_square < distance_cutoff).sum(1) - 1 + np.random.normal(scale=1e-5, size=n_samples)

        self.delta = np.zeros(n_samples)
        parents = np.zeros(n_samples, dtype=int)
        idx_max_rho = self.rho.argmax()
        idx_all = np.arange(n_samples)

        for i in idx_all:
            if i == idx_max_rho:
                self.delta[i] = dist_square[i, :].max()
                parents[i] = i
            else:
                mask = self.rho > self.rho[i]
                dist_masked = dist_square[i, mask]
                self.delta[i] = dist_masked.min()
                parents[i] = idx_all[mask][dist_masked.argmin()]

        confidence = self.rho * self.delta
        n_center = self.get_n_center(confidence)
        center_threshold = np.sort(confidence)[-n_center]
        center_mask = confidence >= center_threshold
        self.cluster_centers_ = X[center_mask, :]

        parents[center_mask] = center_mask.nonzero()[0]

        self.labels_ = -1 * np.ones(n_samples, dtype=int)
        self.labels_[center_mask] = np.arange(self.cluster_centers_.shape[0])

        # 优化：减少排序次数，直接通过 rho 的索引来更新 labels
        sorted_idx = np.argsort(self.rho)[::-1]
        for i in sorted_idx:
            self.labels_[i] = self.labels_[parents[i]]

        self.labels_[self.rho <= 0] = -1

        return self

    def fit_predict(self, X):
        """Compute cluster centers and predict cluster index for each sample.
        Convenience method; equivalent to calling fit(X) followed by
        predict(X).
        
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            New data to transform.

        Returns
        -------
        labels : array, shape [n_samples,]
            Index of the cluster each sample belongs to.
        """
        return self.fit(X).labels_