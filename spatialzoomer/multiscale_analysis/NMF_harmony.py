import scanpy as sc
import numpy as np
from scipy.sparse import issparse

def NMF_harmony_correction(
    adata,
    harmony_vars="batch",
    librarysize_targetsum=1e4,
    quantile_thresh=0.9999,
    max_value=None,
    theta=1,
    max_iter_harmony=20,
    nonneg_strategy="shift",
    store_key_pca_harmony="X_pca_harmony"
):
    """
    Run Harmony on PCA, then apply Harmony's MoE ridge correction to expression matrix X,
    and optionally clip negative values to 0. Return corrected AnnData copy.

    Parameters
    ----------
    adata : AnnData
        Input AnnData.
    harmony_vars : str or list[str]
        Column(s) in adata.obs used as Harmony batch covariates (e.g., 'batch').
    librarysize_targetsum : float
        Target sum for library size normalization.
    quantile_thresh : float or None
        Quantile ceiling threshold (e.g. 0.9999). If None, no ceiling.
    max_value : float or None
        Passed to sc.pp.scale(max_value=...). Usually None if you use quantile ceiling.
    theta : float
        Harmony diversity clustering penalty.
    max_iter_harmony : int
        Harmony max iterations.
    nonneg_strategy : str
        Strategy to handle negative values after correction:
        - "clip": set negative values to 0.
        - "shift": shift all values by adding a small constant to make min value slightly > 0.
        - "none" or "keep": do not modify negative values.
    store_key_pca_harmony : str
        Key to store Harmony-corrected PCA embedding in .obsm.

    Returns
    -------
    _adata : AnnData
        Corrected copy of input adata.
    """
    try:
        import harmonypy
    except ImportError as exc:
        raise ImportError(
            "NMF_harmony_correction requires the optional dependency "
            "'harmonypy'. Please install it with:\n"
            "    pip install harmonypy"
        ) from exc
    
    def stdscale_quantile_celing(_adata, max_value=None, quantile_thresh=None):
        sc.pp.scale(_adata, zero_center=False, max_value=max_value)
        if quantile_thresh is not None:
            X = _adata.X
            if issparse(X):
                threshval = np.quantile(np.array(X.todense()).reshape(-1), quantile_thresh)
            else:
                threshval = np.quantile(X.reshape(-1), quantile_thresh)
            _adata.X[_adata.X > threshval] = threshval

    def moe_correct_ridge(Z_orig, R, K, Phi_moe, lamb):
        """
        Z_orig: (features, cells)  (i.e., X.T)
        Harmony objects:
          R: (n_clusters, cells) or (cells, n_clusters) depending on harmonypy version.
          Phi_moe: design matrix used in Harmony (dimensions depend on harmonypy).
        This function matches your current implementation (transpose choices).
        """
        Z_corr = Z_orig.copy()
        W = None
        Phi_Rk = None
        x = None

        for i in range(K):
            Phi_Rk = np.multiply(Phi_moe.T, R[:, i])
            x = np.dot(Phi_Rk, Phi_moe) + lamb
            W = np.dot(np.dot(np.linalg.inv(x), Phi_Rk), Z_orig.T)
            W[0, :] = 0  # do not remove the intercept
            Z_corr -= np.dot(W.T, Phi_Rk)

        Z_cos = Z_corr / np.linalg.norm(Z_corr, ord=2, axis=0)
        return Z_cos, Z_corr, W, Phi_Rk
    
    np.random.seed(14)
    _adata = adata.copy()

    # 1) normalization for PCA input
    anorm = _adata.copy()
    sc.pp.normalize_total(anorm, target_sum=librarysize_targetsum)
    stdscale_quantile_celing(anorm, max_value=max_value, quantile_thresh=quantile_thresh)

    # also apply same scaling+ceiling to _adata.X (as in your code)
    stdscale_quantile_celing(_adata, max_value=max_value, quantile_thresh=quantile_thresh)

    # 2) PCA on normalized data, then put into _adata
    sc.pp.pca(anorm,  zero_center=True)
    _adata.obsm["X_pca"] = anorm.obsm["X_pca"].copy()
    del anorm

    # 3) Harmony on PCA
    harmony_res = harmonypy.run_harmony(
        _adata.obsm["X_pca"],
        _adata.obs,
        harmony_vars,
        max_iter_harmony=max_iter_harmony,
        theta=theta,
        verbose=False
    )

    X_pca_harmony = harmony_res.Z_corr.copy()

    # 4) Apply MoE ridge correction to expression matrix
    X = _adata.X
    if issparse(X):
        XX = X.todense()
    else:
        XX = np.asarray(X)

    _, X_corr, _, _ = moe_correct_ridge(
        XX.T,
        harmony_res.R,
        harmony_res.K,
        harmony_res.Phi_moe,
        harmony_res.lamb,
    )
    X_corr = np.asarray(X_corr.T)

    if nonneg_strategy  == "clip":
        X_corr[X_corr < 0] = 0

    elif nonneg_strategy == "shift":
        qval = np.quantile(X_corr, 0.00001)
        if qval < 0:
            X_corr = X_corr - qval
        X_corr[X_corr < 0] = 0
    elif nonneg_strategy in ["none", "keep"]:
        pass
    else:
        raise ValueError(
            f"Unknown nonneg_strategy: {nonneg_strategy}. "
            "Choose from ['clip', 'shift', 'none']."
        )

    _adata.X = X_corr.copy()
    _adata.obsm[store_key_pca_harmony] = X_pca_harmony.copy()

    return _adata


# adata_corr = NMF_harmony_correction(
#     adata,
#     harmony_vars="batch",
#     librarysize_targetsum=1e4,
#     quantile_thresh=0.9999,
#     theta=1,
#     max_iter_harmony=20
# )
