import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scanpy as sc
import scanpy.external as sce
from scipy.sparse import issparse


def NMF_harmony_correction(
    adata,
    harmony_vars="batch",
    librarysize_targetsum=1e4,
    quantile_thresh=0.9999,
    theta=1,
    max_iter_harmony=20,
    nonneg_strategy="shift",
    n_pcs=50,
    random_state=0,
):
    """
    Apply Harmony-based batch correction to adata.X.

    Input adata.X should contain raw counts. The function:
      1. Normalizes to librarysize_targetsum
      2. Clips at quantile_thresh to reduce outliers
      3. Log1p transforms
      4. Scales (zero-centers) for PCA
      5. Runs PCA
      6. Runs Harmony on PCA embeddings
      7. Reconstructs corrected expression by inverse PCA (undo scaling + undo log1p)
      8. Applies nonneg_strategy to ensure non-negative values
      9. Returns adata copy with corrected X in normalized (non-log) space

    After this function returns, you may apply sc.pp.normalize_total and
    sc.pp.log1p before running NMF.

    Parameters
    ----------
    adata : AnnData
        Input AnnData with raw counts in adata.X.
    harmony_vars : str or list of str
        Column name(s) in adata.obs to use for Harmony batch correction.
    librarysize_targetsum : float
        Target total counts per cell for normalization before Harmony.
    quantile_thresh : float
        Quantile threshold (among non-zero values) for clipping extreme values
        after normalization. Set to 1.0 to disable clipping.
    theta : float
        Harmony diversity clustering penalty parameter.
    max_iter_harmony : int
        Maximum number of Harmony iterations.
    nonneg_strategy : str
        Strategy to ensure non-negative values after inverse PCA reconstruction.
        "shift" : shift by the global minimum so all values >= 0.
        "clip"  : clip values below 0 to 0.
    n_pcs : int
        Number of principal components for PCA + Harmony.
    random_state : int
        Random seed for PCA.

    Returns
    -------
    adata2 : AnnData
        Copy of input AnnData with batch-corrected expression in adata.X
        (normalized, non-log space, non-negative).
    """
    adata2 = adata.copy()

    # 1. Normalize by library size
    sc.pp.normalize_total(adata2, target_sum=librarysize_targetsum)

    # 2. Clip extreme values at the given quantile (computed over non-zero entries)
    if issparse(adata2.X):
        X_dense = adata2.X.toarray().astype(np.float64)
    else:
        X_dense = np.asarray(adata2.X, dtype=np.float64)

    non_zero_vals = X_dense[X_dense > 0]
    if non_zero_vals.size > 0:
        q_thresh = np.quantile(non_zero_vals, quantile_thresh)
    else:
        q_thresh = np.quantile(X_dense, quantile_thresh)
    X_dense = np.clip(X_dense, None, q_thresh)
    adata2.X = X_dense

    # 3. Log1p transform
    sc.pp.log1p(adata2)

    # 4. Scale for PCA (zero-center; stores mean and std in adata.var)
    sc.pp.scale(adata2, zero_center=True)

    # 5. PCA
    n_pcs_actual = min(n_pcs, min(adata2.n_obs, adata2.n_vars) - 1)
    sc.tl.pca(adata2, n_comps=n_pcs_actual, random_state=random_state)

    # 6. Harmony batch correction on PCA embeddings
    sce.pp.harmony_integrate(
        adata2,
        key=harmony_vars,
        theta=theta,
        max_iter_harmony=max_iter_harmony,
        random_state=random_state,
    )

    # 7. Reconstruct corrected expression in scaled space via inverse PCA
    # adata2.varm["PCs"]: (n_genes, n_pcs) - PCA loadings
    # adata2.obsm["X_pca_harmony"]: (n_cells, n_pcs) - Harmony-corrected scores
    loadings = adata2.varm["PCs"]               # (n_genes, n_pcs)
    W_harmony = adata2.obsm["X_pca_harmony"]    # (n_cells, n_pcs)
    X_scaled_corr = W_harmony @ loadings.T       # (n_cells, n_genes)

    # 8. Undo scaling: X_log_corr = X_scaled_corr * std + mean
    mean = adata2.var["mean"].values             # (n_genes,)
    std = adata2.var["std"].values               # (n_genes,)
    std = np.where(std == 0, 1.0, std)           # avoid division by zero
    X_log_corr = X_scaled_corr * std + mean      # (n_cells, n_genes)

    # 9. Undo log1p: back to normalized (non-log) space
    X_norm_corr = np.expm1(X_log_corr)           # (n_cells, n_genes)

    # 10. Apply non-negative strategy
    if nonneg_strategy == "shift":
        min_val = X_norm_corr.min()
        if min_val < 0:
            X_norm_corr -= min_val
    elif nonneg_strategy == "clip":
        X_norm_corr = np.clip(X_norm_corr, 0, None)
    else:
        raise ValueError(f"Unknown nonneg_strategy: '{nonneg_strategy}'. Use 'shift' or 'clip'.")

    adata2.X = X_norm_corr.astype(np.float32)

    return adata2
