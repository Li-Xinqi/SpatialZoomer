import time
import scanpy as sc
from sklearn.decomposition import NMF
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix, issparse
from pygsp import graphs, filters
import numpy as np
import h5py
import scanpy.external as sce


def performDR(adata, type = 'NMF', n_components=50, batch_col = None, time_report = True):
    """
    Perform dimensionality reduction on the input data using either NMF (Non-negative Matrix Factorization) 
    or PCA (Principal Component Analysis). The reduced data is stored in the adata.obsm, with 'X_nmf' for NMF and 'X_pca' for PCA.

    Parameters：
    ----------
    adata: AnnData
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    type: str
        The type of dimensionality reduction to perform. 
        Options: 'NMF' or 'PCA'.
    n_components: int
        The number of components to keep.
    time_report: bool
        Whether to report the time taken to perform the dimensionality reduction.
    
    Returns:
    --------
    adata : AnnData
        The input data with results of the dimensionality reduction. 
        For NMF, the following fields will be added:
        - adata.obsm['X_nmf']: Transformed matrix (W) representing the NMF scores for each cell.
        - adata.varm['nmf_loadings']: Loadings matrix (H.T) representing the NMF components for each feature.
        For PCA, the principal components will be added to adata as adata.obsm['X_pca'].

    Notes:
    ------
    - For NMF, the sklearn implementation of NMF is used, with initialization set to 'nndsvd' for improved convergence.
    - For PCA, the built-in scanpy function is used (sc.pp.pca).
    
    Example Usage:
    --------------
    # Perform NMF on the data with 50 components
    adata = performDR(adata, type='NMF', n_components=50)

    # Perform PCA on the data with 50 components and report the time
    adata = performDR(adata, type='PCA', n_components=50, time_report=True)
    """
    start_time = time.time()
    if type == 'NMF':
        model = NMF(n_components=n_components, init='nndsvd', random_state=0, l1_ratio=0.5)
        W = model.fit_transform(adata.X)
        H = model.components_

        adata.obsm['X_nmf'] = W
        adata.varm['nmf_loadings'] = H.T
        end_time = time.time()
        elapsed_time = end_time - start_time
        if time_report:
            print(f"Time taken to perform NMF: {elapsed_time:.4f} seconds")
    if type == 'PCA':
        sc.tl.pca(adata, n_comps=n_components)
        if batch_col is not None:
            sc.external.pp.harmony_integrate(adata, key=batch_col)
            adata.obsm['X_pca_raw'] = adata.obsm['X_pca'].copy()
            adata.obsm['X_pca'] = adata.obsm['X_pca_harmony'].copy()
        end_time = time.time()
        elapsed_time = end_time - start_time
        if time_report:
            print(f"Time taken to perform PCA: {elapsed_time:.4f} seconds")
    return adata 



def getKNN(adata, use_rep = 'spatial', name = 'spatial_knn', n_neighbors = 15, factor = 100, metric = 'euclidean', max_similarity = None, pattern = None, results_report = True, time_report= True):
    """
    Construct a KNN graph based on the input data and store the results in the adata.obsp.  
    The similarities between cells are calculated based on their distances, and the KNN graph is symmetrized.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix containing the spatial or other representation in adata.obsm to be used for KNN construction.
    
    use_rep : str, optional, default: 'spatial'
        The key in adata.obsm specifying which data representation to use for calculating KNN. 
        For example, 'spatial' will use the spatial coordinates stored in adata.obsm['spatial'].
    
    name : str, optional, default: 'spatial_knn'
        The name under which the resulting KNN similarity matrix will be stored in adata.obsp.
    
    n_neighbors : int, optional, default: 15
        The number of nearest neighbors to include when constructing the KNN graph (excluding the point itself).
    
    factor : int or float, optional, default: 100
        A factor used to compute similarities by dividing the factor by the distances (with small epsilon 1e-10 added to avoid division by zero).
    
    metric : str, optional, default: 'euclidean'
        The distance metric to use for finding nearest neighbors. Can be 'euclidean', 'manhattan', 'cosine', etc.
    
    results_report : bool, optional, default: True
        If True, prints out a summary of the computed similarities (maximum, minimum, and median values).
    
    time_report : bool, optional, default: True
        If True, reports the time taken to compute the KNN graph.
    
    Returns:
    --------
    adata : AnnData
        The input adata with the KNN similarity matrix added to adata.obsp under the specified name.

    Notes:
    ------
    - The function computes the nearest neighbors using the specified metric and transforms the distances into similarities.
    - The KNN graph is symmetrized so that if there is an edge from cell A to cell B, there will also be an edge from B to A.

    Example Usage:
    --------------
    # Construct a KNN graph using spatial coordinates and store it in adata.obsp['spatial_knn']
    adata = getKNN(adata, use_rep='spatial', n_neighbors=15, factor=100, metric='euclidean', results_report=True)

    """
    start_time = time.time()
    coords = adata.obsm[use_rep].copy()
    nbrs = NearestNeighbors(n_neighbors= n_neighbors + 1, algorithm='auto', metric=metric).fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    # remove the first column, which is the distance to the point itself
    distances = distances[:, 1:]
    indices = indices[:, 1:]

        
    if pattern == 'auto':
        factor = 5 * np.median(distances)
        if max_similarity is None and name == 'spatial_knn':
            max_similarity = 100.0

    similarities  = factor / (distances + 1e-10)
    
    if max_similarity is None:
        max_similarity = np.max(similarities[similarities < 1e10])  

    similarities = np.clip(similarities, None, max_similarity)  # 将所有值限制为不超过 max_similarity



    if results_report:
        # print('Maximum similarities: ', np.max(similarities))
        # print('Minimum similarities: ', np.min(similarities))
        print('Median similarities: ', np.median(similarities))

    n_cells = coords.shape[0]
    row_indices = np.repeat(np.arange(coords.shape[0]), n_neighbors)
    col_indices = indices.flatten()
    data = similarities.flatten()
    
    knn_similarity_matrix = csr_matrix((data, (row_indices, col_indices)), shape=(n_cells, n_cells))
    knn_similarity_matrix = knn_similarity_matrix.maximum(knn_similarity_matrix.T)
    adata.obsp[name] = knn_similarity_matrix.copy()

    if time_report:
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Time taken to construct {name}: {elapsed_time:.4f} seconds")

    return adata



def SGFB(G, data, scales=[2, 4, 8, 16], base_type = 'Heat', time_report = True):
    """
    Perform Spectral Graph Filter Bank (SGFB) on the input data based on a graph structure and given scales.

    Parameters:
    -----------
    G : Graph
        The input graph on which the SGFB will be applied. Typically, this is a graph object 
        (such as from the PyGSP library) that defines the nodes and edges of the graph, and it contains 
        the necessary adjacency or Laplacian matrix for processing.
    
    data : ndarray
        The input data to be transformed. This could be node features or signals aligned with the nodes of the graph.
    
    scales : list of int, optional, default: [2, 4, 8, 16]
        The scales at which the SGFB will be applied. Larger scales correspond to smoother signals, 
        while smaller scales focus on localized details.
    
    base_type : str, optional, default: 'Heat'
        The type of filter to be used as the base for the SGFB. Currently, the 'Heat' kernel is supported.
    
    time_report : bool, optional, default: True
        If True, the function will report the time taken to complete the SGFB.

    Returns:
    --------
    transformed_signals : ndarray
        The transformed signals corresponding to the input data, where each row represents the transformed signal 
        for a given node and each column corresponds to a different scale.

    Example Usage:
    --------------
    # Apply SGFB to a graph G with node data and scales [2, 4, 8, 16]
    transformed_signals = SGFB(G, data, scales=[2, 4, 8, 16])

    Notes:
    ------
    - The heat kernel is commonly used in graph signal processing to diffuse or smooth signals over a graph.
    - The filter scales control the degree of smoothing applied to the signal: smaller scales focus on local features, 
      while larger scales smooth the signal over a wider neighborhood.

    """

    start_time = time.time()
    if base_type == 'Heat':
        g = filters.Heat(G, scales)
    transformed_signals = g.filter(data)
    if time_report:
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Time taken to perform spectral graph filter bank: {elapsed_time:.4f} seconds")
    return transformed_signals



def check_graph(G_pygsp):
    print('is_connected: ', G_pygsp.is_connected())
    print('is_directed: ', G_pygsp.is_directed())
    print('check_weights: ', G_pygsp.check_weights())



def performSGFB(adata, use_knn = 'spatial_knn', use_rep = 'X_nmf', scales = None, runLabel = 'Xenium', save_path = '/'):
    """
    Perform Spectral Graph Filter Bank (SGFB) on the signals (like expressions) and save the transformed signals.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix. This function will use the KNN graph stored in adata.obsp and the specified representation after 
        dimensionality reduction in adata.obsm to perform the heat kernel-based spectral graph filter bank.

    use_knn : str, optional, default: 'spatial_knn'
        The key in adata.obsp that refers to the precomputed K-nearest neighbors (KNN) graph that will be used to construct 
        the graph for the spectral graph signal processing.

    use_rep : str, optional, default: 'X_nmf'
        The key in adata.obsm that refers to the data representation (e.g., NMF-transformed data) on which the heat kernel-based 
        spectral graph filter bank will be applied.

    scales : list of float, optional, default: None
        The list of scales at which to apply the SGFB. If no scales are provided, a default set of scales is generated, 
        which includes smaller and larger scales to cover a wide range of resolutions.

    runLabel : str, optional, default: 'Xenium'
        A label that will be used to name the output file containing the transformed signals.

    save_path : str, optional, default: '/'
        The directory path where the transformed signals will be saved. The output is saved as an HDF5 file.

    Returns:
    --------
    transformed_signals : ndarray
        The transformed signals for the input data at different scales.
    
    scales : list of float
        The list of scales used for the SGFB.

    Notes:
    ------
    - The function retrieves the KNN graph from adata.obsp and uses it to construct a graph object with PyGSP.
    - The graph Laplacian is computed before applying the SGFB, which uses the selected or default scales.
    - The resulting transformed signals are saved in an HDF5 file for later analysis or further processing.

    Example Usage:
    --------------
    # Perform spectral graph filter bank on adata using spatial KNN and NMF-transformed data with default scales
    transformed_signals, scales = performSGFB(adata, use_knn='spatial_knn', use_rep='X_nmf', runLabel='analysis', save_path='/output')

    """

    if scales is None:
        scales = [0.01] + np.arange(0.1, 2.1, 0.1).tolist() + np.arange(2.5, 15.5, 0.5).tolist() + np.arange(16, 21, 1).tolist() +  np.arange(25, 55, 5).tolist()
        scales = [round(x, 2) for x in scales]
    
    if use_knn not in adata.obsp:
        raise ValueError(f"KNN graph '{use_knn}' not found in adata.obsp.")

    spatial_knn = adata.obsp[use_knn]
    G_pygsp = graphs.Graph(spatial_knn)
    G_pygsp.compute_laplacian()
    G_pygsp.estimate_lmax() 
    
    
    if use_rep == 'X':
        data = adata.X.copy()
        if issparse(data):
            print("Converting sparse matrix to dense format for adata.X.")
            data = data.toarray()
    else:
        if use_rep not in adata.obsm:
            raise ValueError(f"Data representation '{use_rep}' not found in adata.obsm.")
        data = adata.obsm[use_rep].copy()

    transformed_signals = SGFB(G_pygsp, data = data, scales= scales, time_report = True)

    with h5py.File(save_path +'/'+ runLabel +'_transformed_signals.h5', 'w') as f:
        f.create_dataset('transformed_signals', data=transformed_signals)
    print(f"Transformed signals saved to {save_path +'/'+ runLabel +'_transformed_signals.h5'}")

    return transformed_signals, scales