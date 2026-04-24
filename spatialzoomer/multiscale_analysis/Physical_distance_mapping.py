import warnings
warnings.filterwarnings("ignore")

# Load packages
import numpy as np
import scanpy as sc
import squidpy as sq
import seaborn as sns
import os
import matplotlib.pyplot as plt
from pygsp import graphs, filters
from matplotlib.ticker import MultipleLocator, FuncFormatter


# Set figure parameters
plt.rcParams.update({
    "font.family":'sans-serif',
    "mathtext.fontset":'stix',
    "font.size": 7,
    'pdf.fonttype': 42,
    })

from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def calculate_density_matrices(adata, radius):
    print('Calculate density')
    coords = adata.obsm['spatial'].copy()
    tree = KDTree(coords)
    neighbor_indices = tree.query_ball_tree(tree, r=radius)
    density = np.array([len(neighbors) for neighbors in neighbor_indices])
    return density



import numpy as np
from typing import Optional

# Sample cells with minimum distance constraint
def sample_min_dist_grid(coords_um: np.ndarray,
                         r: float = 500.0,
                         seed: int = 0,
                         order: str = "random",
                         target_n: Optional[int] = None):
    """
    coords_um: (N,2) float, unit µm
    r: minimum distance (µm)
    order: "random" / "raster" / "given"
    target_n: if given, stop sampling after reaching target_n; None means sample maximally
    
    Returns:
      sel_idx: indices of selected points in coords_um (K,)
    """
    coords = np.asarray(coords_um, dtype=np.float32)
    N = coords.shape[0]
    if N == 0:
        return np.array([], dtype=np.int64)

    if target_n is not None:
        target_n = int(target_n)
        if target_n <= 0:
            return np.array([], dtype=np.int64)
        if target_n >= N:
            # Min-distance constraint still applies; target_n is just an upper limit
            pass

    # Traversal order
    if order == "random":
        rng = np.random.default_rng(seed)
        perm = rng.permutation(N)
    elif order == "raster":
        perm = np.lexsort((coords[:, 1], coords[:, 0]))
    elif order == "given":
        perm = np.arange(N)
    else:
        raise ValueError("order must be one of: random, raster, given")

    inv_r = 1.0 / r
    r2 = r * r

    selected = []
    grid = {}  # (gx,gy) -> list of selected indices (within coords array)

    for idx in perm:
        x, y = coords[idx]
        gx = int(np.floor(x * inv_r))
        gy = int(np.floor(y * inv_r))

        ok = True
        for nx in (gx - 1, gx, gx + 1):
            for ny in (gy - 1, gy, gy + 1):
                bucket = grid.get((nx, ny))
                if not bucket:
                    continue
                bx = coords[bucket, 0]
                by = coords[bucket, 1]
                dx = bx - x
                dy = by - y
                if np.any(dx * dx + dy * dy < r2):
                    ok = False
                    break
            if not ok:
                break

        if ok:
            selected.append(idx)
            grid.setdefault((gx, gy), []).append(idx)

            # Early stop when target count reached
            if target_n is not None and len(selected) >= target_n:
                break

    return np.asarray(selected, dtype=np.int64)


import numpy as np
import pandas as pd
from typing import Optional, Union, Sequence, Dict, List


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plot_spatial_sampling_rounds(
    adata,
    obs_key="sample_round_r1000",
    spatial_key="spatial",
    s_bg=2,                 # background point size
    s_fg=10,                # sampled point size
    alpha_bg=0.3,
    alpha_fg=0.9,
    cmap="tab10",           # colormap for rounds
    max_rounds=None,        # plot only the first K rounds (None = all rounds)
    figsize=(6, 6),
):
    coords = adata.obsm[spatial_key]
    rounds = adata.obs[obs_key].to_numpy()

    x = coords[:, 0]
    y = - coords[:, 1]

    fig, ax = plt.subplots(figsize=figsize)

    # -------------------------
    # Background: unsampled cells (grey)
    # -------------------------
    bg_mask = rounds == 0
    ax.scatter(
        x[bg_mask], y[bg_mask],
        s=s_bg,
        c="lightgrey",
        alpha=alpha_bg,
        linewidths=0,
        rasterized=True,
        zorder=1,
        label="Rounds"
    )

    # -------------------------
    # Foreground: sampled cells colored by round
    # -------------------------
    uniq_rounds = np.sort(np.unique(rounds[rounds > 0]))
    if max_rounds is not None:
        uniq_rounds = uniq_rounds[:max_rounds]

    cmap_obj = cm.get_cmap(cmap, len(uniq_rounds))

    for i, rd in enumerate(uniq_rounds):
        m = rounds == rd
        ax.scatter(
            x[m], y[m],
            s=s_fg,
            color=cmap_obj(i),
            alpha=alpha_fg,
            linewidths=0,
            zorder=2 + i,
            label=f"round {rd}"
        )

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_title("Spatial sampling by rounds")

    ax.legend(
        markerscale=1.5,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False
    )

    plt.tight_layout()
    return fig, ax




import numpy as np

def make_signal_matrix(
    adata,
    obs_key="sample_round_per200",
    dtype=np.float32
):
    """
    Returns:
      X: (N, R)  multi-round signal matrix
      rounds: list[int]  actual round numbers used (column order)
    """
    round_vec = adata.obs[obs_key].to_numpy()
    n_cells = round_vec.shape[0]

    # Valid rounds (exclude 0)
    rounds = np.sort(np.unique(round_vec[round_vec > 0]))
    n_rounds = len(rounds)

    # round_id -> column index
    rd2col = {rd: i for i, rd in enumerate(rounds)}

    # Build signal matrix
    X = np.zeros((n_cells, n_rounds), dtype=dtype)

    for rd in rounds:
        col = rd2col[rd]
        X[round_vec == rd, col] = 1.0

    return X, rounds


import numpy as np
import os
# 1) Get seed cells from sampling obs_key and prepare indices
def get_seed_indices_and_roundcols(adata, obs_key="sample_round_per200"):
    r = adata.obs[obs_key].to_numpy()
    seed_idx = np.flatnonzero(r > 0).astype(np.int32)
    round_col = (r[seed_idx] - 1).astype(np.int16)   # round -> column index in transformed dim 2
    return seed_idx, round_col

# 2) KDTree batch neighbor query
from scipy.spatial import cKDTree

def build_tree(coords_um):
    coords = np.asarray(coords_um, dtype=np.float32)
    tree = cKDTree(coords)
    return tree, coords

def query_neighbors_batch(tree, coords, seed_idx, radius_um=500.0):
    seed_coords = coords[seed_idx]
    # Returns list where each element is a list of neighbor indices
    neigh_lists = tree.query_ball_point(seed_coords, r=radius_um)
    return neigh_lists


# 3) Core: export distance + intensity matrix for each seed
import numpy as np

def export_one_seed_npz(
    out_dir,
    seed_global_idx: int,
    seed_cell_id: str,
    round_col: int,
    neighbors: np.ndarray,
    coords: np.ndarray,
    transformed: np.ndarray,  # (N, R, S)
    scales: np.ndarray,       # (S,)
    sort_by_distance: bool = True,
    dtype_intensity = np.float16,
):
    # Compute distances
    dxy = coords[neighbors] - coords[seed_global_idx]
    dist = np.sqrt((dxy * dxy).sum(axis=1)).astype(np.float32)

    if sort_by_distance:
        order = np.argsort(dist)
        neighbors = neighbors[order]
        dist = dist[order]

    # Extract intensity: (n_neighbors, S)
    intens = transformed[neighbors, round_col, :].astype(dtype_intensity, copy=False)

    # Save one file per seed (memory efficient and parallelizable)
    path = os.path.join(out_dir, f"seed_{seed_global_idx}.npz")
    np.savez_compressed(
        path,
        seed_index=np.int32(seed_global_idx),
        seed_cell_id=np.array(seed_cell_id),
        round_col=np.int16(round_col),
        neighbor_index=neighbors.astype(np.int32),
        distance_um=dist,
        intensity=intens,
        scales=scales.astype(np.float32),
    )
    return path, neighbors.size


# 4) Main export function: batch export (parallelizable)
import pandas as pd
from joblib import Parallel, delayed

def export_all_seeds(
    adata,
    transformed,                 # (N, R, S)
    scales,
    out_dir="seed_profiles_npz",
    obs_key="sample_round_per200",
    radius_um=500.0,
    n_jobs=8,                    # number of parallel jobs
    sort_by_distance=True,
    dtype_intensity=np.float16,
    obs_columns=None,            # obs columns to export: None=all columns, []=none
):
    os.makedirs(out_dir, exist_ok=True)

    # 1) seeds & round mapping
    seed_idx, round_col = get_seed_indices_and_roundcols(adata, obs_key=obs_key)
    cell_ids = adata.obs_names.to_numpy()

    # 2) build KDTree once
    tree, coords = build_tree(adata.obsm["spatial"])

    # 3) neighbor query (batch, fast)
    neigh_lists = query_neighbors_batch(tree, coords, seed_idx, radius_um=radius_um)

    scales = np.asarray(scales, dtype=np.float32)

    # 4) parallel export
    def _run(i):
        sidx = int(seed_idx[i])
        rcol = int(round_col[i])
        nbr = np.asarray(neigh_lists[i], dtype=np.int32)
        # Optional: exclude the seed cell itself
        # nbr = nbr[nbr != sidx]
        return export_one_seed_npz(
            out_dir=out_dir,
            seed_global_idx=sidx,
            seed_cell_id=str(cell_ids[sidx]),
            round_col=rcol,
            neighbors=nbr,
            coords=coords,
            transformed=transformed,
            scales=scales,
            sort_by_distance=sort_by_distance,
            dtype_intensity=dtype_intensity,
        )

    results = Parallel(n_jobs=n_jobs, prefer="processes", batch_size=16)(
        delayed(_run)(i) for i in range(len(seed_idx))
    )

    # 5) write manifest
    manifest = pd.DataFrame({
        "seed_index": seed_idx,
        "seed_cell_id": cell_ids[seed_idx],
        "round_col": round_col,
        "file": [p for p, k in results],
        "n_neighbors": [k for p, k in results],
        "radius_um": radius_um,
    })
    
    # Append adata.obs columns to manifest
    if obs_columns is None:
        # Export all obs columns
        obs_data = adata.obs.loc[cell_ids[seed_idx]].reset_index(drop=True)
        for col in obs_data.columns:
            if col not in manifest.columns:
                manifest[col] = obs_data[col].values
    elif isinstance(obs_columns, (list, tuple)) and len(obs_columns) > 0:
        # Export specified columns
        obs_data = adata.obs.loc[cell_ids[seed_idx], obs_columns].reset_index(drop=True)
        for col in obs_columns:
            if col in adata.obs.columns:
                manifest[col] = obs_data[col].values
    
    manifest_path = os.path.join(out_dir, "manifest.csv")
    manifest.to_csv(manifest_path, index=False)

    return manifest_path



# Classify cells into High / Medium / Low groups by density quantile
import pandas as pd

def classify_density_by_quantile(adata, density_key='Density', percentile=10, new_key='Density_group'):
    """
    Classify cells into three groups (High, Medium, Low) based on density quantiles.
    
    Parameters:
    -----------
    adata : AnnData
        Input AnnData object.
    density_key : str
        Column name for density values, default 'Density'.
    percentile : int
        Quantile parameter, default 10.
        - If 10, uses q10 and q90 as thresholds.
        - If 25, uses q25 and q75 as thresholds.
    new_key : str
        Name for the new grouping column, default 'Density_group'.
    
    Returns:
    --------
    q_low, q_high : tuple
        Lower and upper quantile values.
    """
    density_values = adata.obs[density_key].values
    
    # Compute quantile thresholds
    q_low = np.percentile(density_values, percentile)
    q_high = np.percentile(density_values, 100 - percentile)
    
    def classify_density(density_value):
        if density_value >= q_high:
            return "High"
        elif density_value <= q_low:
            return "Low"
        else:
            return "Medium"
    
    # Classify cells
    adata.obs[new_key] = adata.obs[density_key].apply(classify_density)
    adata.obs[new_key] = pd.Categorical(
        adata.obs[new_key], 
        categories=["Low", "Medium", "High"], 
        ordered=True
    )
    
    # Print statistics
    print(f"Using percentile: {percentile} (q{percentile}={q_low:.2f}, q{100-percentile}={q_high:.2f})")
    print(adata.obs[new_key].value_counts())
    
    return q_low, q_high

import numpy as np
import pandas as pd
from typing import Optional, Sequence, Dict, List, Set


def multi_round_sampling_ultrafast(
    adata,
    r: float = 500.0,
    target_total: int = 500,                 # final total sample count (strict)
    base_seed: int = 0,
    key_spatial: str = "spatial",
    obs_key_round: str = "sample_round",
    obs_key_flag: str = "sampled_any",

    # Round 1 acceleration: no strict upper limit
    first_round_target_n: int = 200,         # stop round 1 when this count is reached
    first_round_candidate_cap: Optional[int] = 200000,  # candidate pool size for round 1 (None=all)

    # Subsequent rounds: target per round (multiple of 25); None = auto-estimate from round 1
    per_round_25: Optional[int] = None,

    # Preferred cell sampling (optional)
    preferred_ids: Optional[Sequence[str]] = None, # cell names
    preferred_max_frac: float = 0.5,         # max fraction of preferred cells per round

    # Output
    save_csv_path: Optional[str] = None,
    verbose: bool = False,
) -> Dict[int, List[str]]:

    if key_spatial not in adata.obsm:
        raise KeyError(f"adata.obsm does not contain '{key_spatial}'")

    coords_all = np.asarray(adata.obsm[key_spatial], dtype=np.float32)
    ids_all = adata.obs_names.to_numpy()
    n_cells = len(ids_all)

    adata.obs[obs_key_round] = 0
    adata.obs[obs_key_flag] = False

    id2idx = {cid: i for i, cid in enumerate(ids_all)}

    preferred_idx: Optional[Set[int]] = None
    if preferred_ids is not None:
        preferred_idx = {id2idx[cid] for cid in preferred_ids if cid in id2idx}

    remaining_mask = np.ones(n_cells, dtype=bool)
    round2ids: Dict[int, List[str]] = {}

    rng = np.random.default_rng(base_seed)

    # -------- Round 1: candidate pool + early stop --------
    idx0 = np.arange(n_cells, dtype=np.int32)
    if first_round_candidate_cap is not None and first_round_candidate_cap < n_cells:
        idx0 = rng.choice(idx0, size=int(first_round_candidate_cap), replace=False)

    coords0 = coords_all[idx0]
    sel0_local = sample_min_dist_grid(
        coords0, r=r, seed=base_seed + 1, order="random", target_n=int(first_round_target_n)
    )
    sel0 = idx0[sel0_local]
    sel0_ids = ids_all[sel0].tolist()

    round2ids[1] = sel0_ids
    adata.obs.loc[sel0_ids, obs_key_round] = 1
    adata.obs.loc[sel0_ids, obs_key_flag] = True
    remaining_mask[sel0] = False

    total = len(sel0_ids)
    if verbose:
        print(f"[Round 01] sampled {len(sel0_ids)} cells | total={total}/{target_total}")

    if total >= target_total:
        # First round is sufficient (edge case): truncate to target_total
        round2ids[1] = round2ids[1][:target_total]
        if save_csv_path is not None:
            rows = [(1, cid) for cid in round2ids[1]]
            pd.DataFrame(rows, columns=["round", "cell_id"]).to_csv(save_csv_path, index=False)
        return round2ids

    # -------- decide per_round_25 --------
    if per_round_25 is None:
        # Estimate per_round_25 from round 1: <= half of round 1, multiple of 25, at least 25
        n0 = len(sel0_ids)
        per_round_25 = max((n0 // 25) * 25, 25)

    if per_round_25 < 1:
        per_round_25 = 25

    # -------- subsequent rounds: fill to target_total --------
    rd = 2
    while total < target_total:
        idx_remain = np.flatnonzero(remaining_mask)
        if idx_remain.size == 0:
            if verbose:
                print("[WARN] no remaining cells to sample.")
            break

        remaining_need = target_total - total

        # Key: last round fills exactly to target_total (no multiple-of-25 requirement)
        if remaining_need <= per_round_25:
            target_n = int(remaining_need)           # last round exact fill
        else:
            target_n = int(per_round_25)             # regular round 25-multiple

        # --- build ordered candidates (preferred first up to half) ---
        rng_rd = np.random.default_rng(base_seed + rd)
        rng_rd.shuffle(idx_remain)

        if preferred_idx is not None:
            is_pref = np.array([i in preferred_idx for i in idx_remain], dtype=bool)
            pref = idx_remain[is_pref]
            other = idx_remain[~is_pref]
            rng_rd.shuffle(pref)
            rng_rd.shuffle(other)

            n_pref_max = int(target_n * preferred_max_frac)
            ordered = np.concatenate([pref[:n_pref_max], other])
        else:
            ordered = idx_remain

        coords_use = coords_all[ordered]

        sel_local = sample_min_dist_grid(
            coords_use, r=r, seed=base_seed + rd, order="given", target_n=target_n
        )
        sel_global = ordered[sel_local]
        sel_ids = ids_all[sel_global].tolist()

        round2ids[rd] = sel_ids
        if len(sel_ids) > 0:
            adata.obs.loc[sel_ids, obs_key_round] = rd
            adata.obs.loc[sel_ids, obs_key_flag] = True
            remaining_mask[sel_global] = False

        total += len(sel_ids)
        if verbose:
            print(f"[Round {rd:02d}] sampled {len(sel_ids)} cells | total={total}/{target_total} "
                  f"(target_n={target_n})")

        # If 0 cells sampled (spatial constraint too strict or poor candidates), break to avoid infinite loop
        if len(sel_ids) == 0:
            if verbose:
                print("[WARN] 0 cells sampled this round. Consider lowering r or increasing candidate pool.")
            break
        rd += 1

    # -------- optional save --------
    if save_csv_path is not None:
        rows = [(rdi, cid) for rdi, cids in round2ids.items() for cid in cids]
        pd.DataFrame(rows, columns=["round", "cell_id"]).to_csv(save_csv_path, index=False)

    return round2ids



# heatmap
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def _bin_edges_centers(x_max_um, bin_um):
    edges = np.arange(0, x_max_um + bin_um, bin_um, dtype=np.float32)
    centers = (edges[:-1] + edges[1:]) * 0.5
    return edges, centers

def _scale_edges(scales):
    """Y-axis bin edges for pcolormesh (works for non-uniform scales too)"""
    s = np.asarray(scales, dtype=np.float32)
    if s.size == 1:
        return np.array([s[0]-0.5, s[0]+0.5], dtype=np.float32)
    mids = (s[:-1] + s[1:]) * 0.5
    first = s[0] - (mids[0] - s[0])
    last  = s[-1] + (s[-1] - mids[-1])
    return np.concatenate([[first], mids, [last]]).astype(np.float32)

def aggregate_distance_scale_heatmap(
    manifest_csv: str,
    x_max_um: float = 500.0,
    bin_um: float = 5.0,
    intensity_cut: float = 0.0,          # optional: treat intensities below this as invalid (excluded from interpolation)
    normalize_by_self: bool = False,      # normalize each seed's intensity per scale by its self-intensity
    normalize_rows: bool = True,         # row-normalize so each scale's distance bins sum to 1
    use_counts_weight: bool = True,      # weight by neighbor count in each distance bin
    smooth_sigma_bins: float = 0.0,      # optional: 1D Gaussian smoothing after interpolation (in bins)
    out_npy: str = None,                 # optional: save aggregated matrix as .npy
    energy_percentile: float = 0.95      # energy percentile threshold (0-1) for distance calculation
):
    """
    Returns:
      H: (S, B)  aggregated heatmap (scale x distance_bin_center)
      scales: (S,)
      dist_centers: (B,)
    """
    mf = pd.read_csv(manifest_csv)
    files = mf["file"].tolist()

    edges, centers = _bin_edges_centers(x_max_um, bin_um)
    B = centers.size

    # Load first file to get scales / S
    d0 = np.load(files[0], allow_pickle=True)
    scales = d0["scales"].astype(np.float32)
    S = scales.size
    d0.close()

    # Accumulators: sum_wI / sum_w
    sum_wI = np.zeros((S, B), dtype=np.float32)
    sum_w  = np.zeros((B,), dtype=np.float32) if use_counts_weight else np.zeros((S, B), dtype=np.float32)

    # Optional smoothing (applied to interpolated result only)
    if smooth_sigma_bins > 0:
        try:
            from scipy.ndimage import gaussian_filter1d
        except Exception:
            raise ImportError("smooth_sigma_bins>0 requires scipy (scipy.ndimage.gaussian_filter1d)")

    for fp in files:
        d = np.load(fp, allow_pickle=True)
        dist = d["distance_um"].astype(np.float32)          # (K,)
        nbr  = d["neighbor_index"].astype(np.int32)         # (K,)
        intens = d["intensity"].astype(np.float32)          # (K,S) or (K,S) after load
        seed_index = int(d["seed_index"])
        d.close()

        # Distance cutoff
        m = dist <= float(x_max_um)
        if not np.any(m):
            continue
        dist = dist[m]
        nbr = nbr[m]
        intens = intens[m, :]  # (K', S)

        # Sort by ascending distance for np.interp
        o = np.argsort(dist)
        dist = dist[o]
        nbr = nbr[o]
        intens = intens[o, :]

        # Optional: intensity threshold - treat values below this as invalid to reduce noise
        if intensity_cut > 0:
            intens = np.where(intens >= float(intensity_cut), intens, np.nan)

        # --- normalize by self (set seed's own intensity per scale = 1) ---
        if normalize_by_self:
            pos = np.where(nbr == seed_index)[0]
            if pos.size > 0:
                base = intens[pos[0], :]  # (S,)
            else:
                # If the seed itself was excluded from export, skip normalization
                # to avoid losing the seed
                base = None

            if base is not None:
                base = np.maximum(base, 1e-12)
                intens = intens / base

        # --- Weight: neighbor count per distance bin ---
        if use_counts_weight:
            counts, _ = np.histogram(dist, bins=edges)
            counts = counts.astype(np.float32)  # (B,)
            if counts.sum() == 0:
                continue
        else:
            counts = None

        # --- For each scale: interpolate discrete (dist, intensity) to bin centers ---
        # np.interp does not handle nan; remove nan per scale before interpolating
        Igrid = np.empty((S, B), dtype=np.float32)
        for si in range(S):
            y = intens[:, si]
            ok = np.isfinite(y)
            if ok.sum() < 2:
                Igrid[si, :] = np.nan
                continue
            Igrid[si, :] = np.interp(centers, dist[ok], y[ok], left=np.nan, right=np.nan)

        # Optional: smooth along distance axis after interpolation
        if smooth_sigma_bins > 0:
            # gaussian_filter1d does not handle nan: fill nan with 0, then divide by smoothed mask
            mask = np.isfinite(Igrid).astype(np.float32)
            I0 = np.nan_to_num(Igrid, nan=0.0)
            I_s = gaussian_filter1d(I0, sigma=smooth_sigma_bins, axis=1, mode="nearest")
            m_s = gaussian_filter1d(mask, sigma=smooth_sigma_bins, axis=1, mode="nearest")
            Igrid = np.where(m_s > 1e-6, I_s / m_s, np.nan)

        # --- Aggregate ---
        if use_counts_weight:
            # counts depend only on distance, same weight for all scales
            # accumulate only where Igrid is valid
            valid = np.isfinite(Igrid)
            w = counts[None, :]  # (1,B)
            sum_wI[valid] += (Igrid * w)[valid]
            sum_w += counts
        else:
            # Equal-weight average (each seed contributes weight=1 per bin)
            valid = np.isfinite(Igrid)
            sum_wI[valid] += Igrid[valid]
            sum_w[valid] += 1.0

    # Compute weighted mean
    if use_counts_weight:
        denom = np.maximum(sum_w[None, :], 1e-12)
        H = sum_wI / denom
    else:
        denom = np.maximum(sum_w, 1e-12)
        H = sum_wI / denom
    
    # Optional: row-normalize (each scale's distance bins sum to 1)
    if normalize_rows:
        for si in range(S):
            row_sum = np.nansum(H[si, :])  # ignore nan
            if row_sum > 1e-12:
                H[si, :] = H[si, :] / row_sum

    # ---- Compute energy percentile distance ----
    # sum_wI: (S, B) total unnormalized energy (some bins may be invalid/nan)
    E_percentile_dist = np.full(S, np.nan, dtype=np.float32)
    for si in range(S):
        E = sum_wI[si, :]  # (B,) total energy distribution for this scale
        valid_idx = np.isfinite(E)  # which distance bins have valid energy
        
        if np.any(valid_idx):
            E_valid = E[valid_idx]  # extract valid values
            total_energy = E_valid.sum()
            
            if total_energy > 0:
                # Compute cumulative energy fraction over valid bins
                E_norm = E_valid / total_energy  # normalize to [0,1]
                E_cum = np.cumsum(E_norm)  # cumulative energy
                
                # Find the first distance bin where cumulative energy >= energy_percentile
                hit = E_cum >= energy_percentile
                if np.any(hit):
                    idx_in_valid = np.argmax(hit)  # position within valid array
                    dist_idx = np.where(valid_idx)[0][idx_in_valid]  # map back to full distance array
                    E_percentile_dist[si] = centers[dist_idx]
                # else: keep nan

    # ---- Handle scale 0: intensity = 1 at distance 0, 0 elsewhere ----
    scale0_row = np.zeros((1, B), dtype=np.float32)
    # Find the bin closest to distance 0
    dist0_idx = np.argmin(np.abs(centers))
    scale0_row[0, dist0_idx] = 1.0
    
    # Check if scale 0 already exists
    scale0_exists = np.any(np.isclose(scales, 0.0, atol=1e-6))
    
    if scale0_exists:
        # If scale 0 exists, find its index and overwrite
        scale0_pos = np.where(np.isclose(scales, 0.0, atol=1e-6))[0][0]
        H[scale0_pos, :] = scale0_row[0, :]
        E_percentile_dist[scale0_pos] = 0.0
    else:
        # Scale 0 does not exist yet: prepend to H matrix
        H = np.vstack([scale0_row, H])  # (S+1, B)
        # Prepend 0 to scales array
        scales = np.concatenate([[0.0], scales])  # (S+1,)
        # Prepend 0 to E_percentile_dist (energy distance at scale 0 is 0)
        E_percentile_dist = np.concatenate([[0.0], E_percentile_dist])  # (S+1,)

    if out_npy is not None:
        np.save(out_npy, H)

    return H, scales, centers, E_percentile_dist

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, FixedLocator

def plot_distance_scale_heatmap(
    H=None, scales=None, dist_centers=None,  # data can be passed directly
    E95_dist=None,
    manifest_csv=None,                   # or compute from manifest
    bin_um=1.0,                          # parameter for manifest mode
    important_scales=None,               # optional: list of important scales to mark on the E95 curve
    intensity_cut=0.0,
    normalize_by_self=False,
    normalize_rows=True,                # row-normalize: each scale's distance bins sum to 1
    use_counts_weight=False,
    smooth_sigma_bins=0.0,
    energy_percentile=0.95,
    x_max_um=250.0,
    y_max_scale=None,                    # y-axis (scale) maximum cutoff
    cmap="turbo",
    vmin=None, vmax=None,
    figsize=(6, 4.0),
    show_colorbar=True,
    ylabel="Scale",
    xlabel="Distance (µm)",
    title="Scale-distance heatmap",
    E95_color="gold",                    # color for the 95% energy line
    colorbar_scale="linear",             # "linear" / "log" / "log2" / "symlog" / "power"
    colorbar_vmin_log=1e-5,              # minimum value for log scale (values below use darkest color)
    colorbar_vmax_symlog=0.8,            # maximum value for symlog scale (values above use brightest color)
    colorbar_power=0.5,                  # exponent for power scale (0.5=sqrt, 0.3=gentler, 0.7=steeper)
    save_pdf=True,                       # whether to save PDF
    dpi=300,
    save_data=True,                      # whether to save H matrix and E95 data
    out_dir=None,                        # output path; None uses the manifest directory
):
    """
    Plot scale-distance heatmap.

    Two usage modes:
    1. Pass pre-computed H, scales, dist_centers, E95_dist directly.
    2. Pass manifest_csv; the function calls aggregate_distance_scale_heatmap internally.

    Returns
    -------
    fig, ax : matplotlib objects
    results : dict (only in manifest_csv mode)
    """
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.ticker import MultipleLocator, FuncFormatter

    # If manifest_csv is provided, compute the H matrix first
    if manifest_csv is not None:
        mf = pd.read_csv(manifest_csv)
        print(f"Computing heatmap from manifest ({mf.shape[0]} seeds)...")
        H, scales, dist_centers, E95_dist = aggregate_distance_scale_heatmap(
            manifest_csv=manifest_csv,
            x_max_um=x_max_um,
            bin_um=bin_um,
            intensity_cut=intensity_cut,
            normalize_by_self=normalize_by_self,
            normalize_rows=normalize_rows,
            use_counts_weight=use_counts_weight,
            smooth_sigma_bins=smooth_sigma_bins,
            energy_percentile=energy_percentile
        )

        # Return results (numpy arrays only)
        results = {
            "H": np.asarray(H),
            "E95": np.asarray(E95_dist),
            "scales": np.asarray(scales),
            "dist_centers": np.asarray(dist_centers),
        }

        # Save output data
        if save_data:
            if out_dir is None:
                out_dir = os.path.dirname(manifest_csv)
            os.makedirs(out_dir, exist_ok=True)

            # 1) CSV (for human review / Excel)
            H_df = pd.DataFrame(
                H,
                index=[f"scale_{s:.2f}" for s in scales],
                columns=[f"dist_{d:.2f}" for d in dist_centers]
            )
            E95_df = pd.DataFrame({"scale": scales, "E95_distance_um": E95_dist})

            H_df.to_csv(os.path.join(out_dir, "H_matrix.csv"))
            E95_df.to_csv(os.path.join(out_dir, "E95_distance.csv"), index=False)

            # 2) NPZ (for direct loading in Python)
            npz_path = os.path.join(out_dir, "heatmap_arrays.npz")
            np.savez_compressed(
                npz_path,
                H=np.asarray(H, dtype=np.float32),
                E95=np.asarray(E95_dist, dtype=np.float32),
                scales=np.asarray(scales, dtype=np.float32),
                dist_centers=np.asarray(dist_centers, dtype=np.float32),
            )

            print(f"[Saved] NPZ -> {npz_path}")
            print(f"[Saved] CSV -> {out_dir}/H_matrix.csv, {out_dir}/E95_distance.csv")
    else:
        if H is None or scales is None or dist_centers is None:
            raise ValueError("Either provide (H, scales, dist_centers) or manifest_csv")
        results = None

    H = np.asarray(H, dtype=np.float32)       # (S,B)
    scales = np.asarray(scales, dtype=np.float32)
    dist_centers = np.asarray(dist_centers, dtype=np.float32)

    # edges for pcolormesh
    dx = dist_centers[1] - dist_centers[0] if dist_centers.size > 1 else 1.0
    x_edges = np.concatenate([[dist_centers[0] - dx / 2], dist_centers + dx / 2]).astype(np.float32)
    y_edges = _scale_edges(scales)

    fig, ax = plt.subplots(figsize=figsize)

    # Preprocessing: for log/log2, avoid NaN/near-zero values causing white patches
    H_plot = H.copy()
    if colorbar_scale in ("log", "log2"):
        vmin_val = float(colorbar_vmin_log)
        H_plot = np.where((H_plot < vmin_val) & np.isfinite(H_plot), vmin_val, H_plot)
        H_plot = np.nan_to_num(H_plot, nan=vmin_val)

    # ------------------ Select norm + plot ------------------
    if colorbar_scale == "log":
        from matplotlib.colors import LogNorm
        norm = LogNorm(vmin=float(colorbar_vmin_log),
                       vmax=vmax if vmax is not None else float(np.nanmax(H_plot)))
        mesh = ax.pcolormesh(x_edges, y_edges, H_plot, shading="auto", cmap=cmap, norm=norm)

    elif colorbar_scale == "log2":
        from matplotlib.ticker import FixedLocator

        vmin_val = float(colorbar_vmin_log)
        vmax_val = float(vmax) if vmax is not None else float(np.nanmax(H_plot))

        H_plot_log2 = np.log2(np.maximum(H_plot, vmin_val))
        vmin_log2 = np.log2(vmin_val)
        vmax_log2 = np.log2(vmax_val)

        mesh = ax.pcolormesh(
            x_edges, y_edges, H_plot_log2,
            shading="auto", cmap=cmap,
            vmin=vmin_log2, vmax=vmax_log2
        )

    elif colorbar_scale == "power":
        from matplotlib.colors import PowerNorm
        vmin_pow = float(vmin) if vmin is not None else float(colorbar_vmin_log)
        vmax_pow = float(vmax) if vmax is not None else float(np.nanmax(H_plot))
        norm = PowerNorm(gamma=float(colorbar_power), vmin=vmin_pow, vmax=vmax_pow)
        mesh = ax.pcolormesh(x_edges, y_edges, H_plot, shading="auto", cmap=cmap, norm=norm)

    elif colorbar_scale == "symlog":
        # ---- clip to min(vmax_data, colorbar_vmax_symlog); colorbar also only shown up to clip_val ----
        from matplotlib.colors import SymLogNorm

        vmin_global = 0.0
        vmax_data = float(np.nanmax(H_plot[np.isfinite(H_plot)])) if np.any(np.isfinite(H_plot)) else 0.0
        clip_val = min(vmax_data, float(colorbar_vmax_symlog))  # actual color mapping upper bound

        norm = SymLogNorm(
            linthresh=0.01,
            linscale=1.0,
            vmin=vmin_global,
            vmax=clip_val,
            base=10,
            clip=True
        )

        mesh = ax.pcolormesh(
            x_edges, y_edges, H_plot,
            shading="auto",
            cmap=cmap,
            norm=norm
        )

    else:  # linear
        mesh = ax.pcolormesh(
            x_edges, y_edges, H_plot,
            shading="auto",
            cmap=cmap,
            vmin=vmin, vmax=vmax
        )

    # ------------------ Axes style ------------------
    ax.set_xlim(0, float(x_max_um))
    if y_max_scale is not None:
        ax.set_ylim(0, float(y_max_scale))
    ax.set_xlabel(xlabel, fontsize=6)
    ax.set_ylabel(ylabel, fontsize=6)
    ax.set_title(title, fontsize=7)

    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.tick_params(axis="both", which="major", labelsize=6)

    # 95% energy line + important scale markers
    if E95_dist is not None:
        E95_dist = np.asarray(E95_dist, dtype=np.float32)

        # Compatibility: some older versions have E95 one shorter than scales (missing explicit scale=0 row)
        if E95_dist.shape[0] == scales.shape[0]:
            E95_x = E95_dist
            E95_y = scales
        elif E95_dist.shape[0] == scales.shape[0] - 1:
            E95_x = np.concatenate([[0.0], E95_dist]).astype(np.float32)
            E95_y = scales
        else:
            E95_x = None
            E95_y = None

        if E95_x is not None:
            ax.plot(E95_x, E95_y,
                    color=E95_color, linewidth=2.0, alpha=0.85, label="95% energy")

            # Mark important scale points (same color)
            if important_scales is not None:
                try:
                    imp = np.asarray(list(important_scales), dtype=np.float32)
                except Exception:
                    imp = np.asarray(important_scales, dtype=np.float32)

                # Automatically add scale 0
                if imp.size == 0:
                    imp = np.array([0.0], dtype=np.float32)
                else:
                    if not np.any(np.isclose(imp, 0.0, atol=float(1e-4))):
                        imp = np.concatenate([np.array([0.0], dtype=np.float32), imp])
                    # Deduplicate and sort (avoid duplicate markers)
                    imp = np.unique(imp)

                if imp.size > 0:
                    for s in imp:
                        # Find the closest scale
                        idx = int(np.nanargmin(np.abs(E95_y - float(s))))
                        sy = float(E95_y[idx])
                        sx = float(E95_x[idx])

                        # Only mark scales that exactly exist (strict match)
                        if np.abs(sy - float(s)) > float(1e-4):
                            continue

                        ax.scatter(
                            [sx], [sy],
                            s=14,
                            c=E95_color,
                            edgecolors=E95_color,
                            linewidths=0.5,
                            zorder=5
                        )
    # ------------------ colorbar ------------------
    if show_colorbar:
        # extend: symlog indicates truncation with 'max' arrow
        extend_opt = "max" if colorbar_scale == "symlog" else "neither"
        cbar = fig.colorbar(mesh, ax=ax, pad=0.02, fraction=0.05, extend=extend_opt)
        cbar.set_label("Mean intensity", fontsize=6)
        cbar.ax.tick_params(labelsize=6)

        if colorbar_scale == "log":
            cbar.ax.set_ylabel("Mean intensity (log₁₀)", rotation=270, labelpad=5, fontsize=6)

        elif colorbar_scale == "log2":
            from matplotlib.ticker import FixedLocator

            def log2_to_original(x, pos):
                return f"{2**x:.1e}"

            cbar.ax.yaxis.set_major_formatter(FuncFormatter(log2_to_original))
            tick_values = [1e-5, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0]
            tick_positions = [np.log2(v) for v in tick_values]
            cbar.ax.yaxis.set_major_locator(FixedLocator(tick_positions))
            cbar.ax.set_ylabel("Mean intensity (log₂)", rotation=270, labelpad=5, fontsize=6)

        elif colorbar_scale == "symlog":
            # Consistent with by_groups: only display up to clip_val (typically 0.5)
            vmin_global = 0.0
            vmax_data = float(np.nanmax(H_plot[np.isfinite(H_plot)])) if np.any(np.isfinite(H_plot)) else 0.0
            clip_val = min(vmax_data, float(colorbar_vmax_symlog))

            tick_values = [0.0, 1e-2, 0.02, 0.04, 0.06, 0.08, 1e-1, 0.2, 0.3, 0.4, 0.5]
            tick_values = [t for t in tick_values if t >= vmin_global and t <= clip_val]
            if len(tick_values) == 0 or not np.isclose(tick_values[-1], clip_val):
                tick_values.append(clip_val)
            tick_values = sorted(set([float(t) for t in tick_values]))
            cbar.set_ticks(tick_values)

            def symlog_label(v, pos):
                if np.isclose(v, 0.0):
                    return "0"
                return f"{v:.0e}" if v < 0.01 else f"{v:g}"

            cbar.ax.yaxis.set_major_formatter(FuncFormatter(symlog_label))
            cbar.ax.set_ylabel(
                f"Mean intensity (symlog)",
                rotation=270, labelpad=10, fontsize=6
            )

    # legend（E95）
    if E95_dist is not None:
        leg = ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.85), frameon=True, fontsize=6)
        for text in leg.get_texts():
            text.set_color("white")
        frame = leg.get_frame()
        frame.set_facecolor("white")
        frame.set_alpha(0.4)
        frame.set_edgecolor("none")

    plt.tight_layout()

    # Save PDF
    if save_pdf:
        if out_dir is None:
            out_dir = os.path.dirname(manifest_csv) if manifest_csv is not None else "."
        os.makedirs(out_dir, exist_ok=True)
        pdf_path = os.path.join(out_dir, f"scale_distance_heatmap_{colorbar_scale}.pdf")
        fig.savefig(pdf_path, dpi=dpi, bbox_inches="tight")
        print(f"[Saved] PDF -> {pdf_path}")

    if results is not None:
        return fig, ax, results
    else:
        return fig, ax




def plot_distance_scale_heatmap_by_groups(
    manifest_csv: str,
    group_column: str,                      # column name for grouping in manifest
    groups_order: Sequence[str] = ("Low", "Medium", "High"),  # group order
    x_max_um: float = 500.0,
    y_max_scale: Optional[float] = None,
    bin_um: float = 1.0,
    important_scales=None,                  # optional: list of important scales to mark on each E95 curve (scale 0 added automatically)
    intensity_cut: float = 0.0,
    normalize_by_self: bool = False,
    normalize_rows: bool = True,        # row-normalize: each scale's distance bins sum to 1
    use_counts_weight: bool = False,
    smooth_sigma_bins: float = 0.0,
    energy_percentile: float = 0.95,
    cmap: str = "turbo",
    colorbar_scale: str = "linear",
    colorbar_vmin_log: float = 1e-5,
    colorbar_vmax_symlog: float = 0.5,   # true vmax for color clipping in symlog mode
    colorbar_power: float = 0.5,
    figsize: tuple = (8, 6),
    E95_color: str = "gold",
    xlabel: str = "Distance (µm)",
    ylabel: str = "Scale",
    suptitle: Optional[str] = None,          # overall title; None = auto-generate
    save_data: bool = True,                 # whether to save H matrix and E95 data
    out_dir: Optional[str] = None,          # output path; None uses the manifest directory
    save_pdf: bool = True,                  # whether to save PDF
    dpi=300, 
):
    """
    Plot a 2x2 subplot heatmap by groups from manifest_csv.
    - First panel: overall heatmap (all seeds)
    - Remaining three panels: per-group heatmaps
    - Shared colorbar

    In symlog mode:
    - Values > colorbar_vmax_symlog are clipped to the brightest color
    - Colorbar is shown only up to colorbar_vmax_symlog (typically 0.5)
    - A marker line is drawn at 0.5 on the colorbar
    - Value 1 is not displayed and no arrows are drawn
    """
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.ticker import MultipleLocator, FuncFormatter

    # ------------------ Load manifest ------------------
    mf = pd.read_csv(manifest_csv)
    if group_column not in mf.columns:
        raise ValueError(f"Column '{group_column}' not found in manifest_csv")

    # ------------------ 1) Overall heatmap ------------------
    print(f"Computing overall heatmap (all {len(mf)} seeds)...")
    H_all, scales, dist_centers, E95_all = aggregate_distance_scale_heatmap(
        manifest_csv=manifest_csv,
        x_max_um=x_max_um,
        bin_um=bin_um,
        intensity_cut=intensity_cut,
        normalize_by_self=normalize_by_self,
        normalize_rows=normalize_rows,
        use_counts_weight=use_counts_weight,
        smooth_sigma_bins=smooth_sigma_bins,
        energy_percentile=energy_percentile
    )

    # ------------------ 2) Per-group heatmaps ------------------
    H_groups = {}
    E95_groups = {}
    n_seeds_groups = {}

    for grp in groups_order:
        sub_mf = mf[mf[group_column] == grp]
        n_seeds_groups[grp] = len(sub_mf)

        if len(sub_mf) == 0:
            print(f"[WARN] No seeds in group {grp}")
            H_groups[grp] = None
            E95_groups[grp] = None
            continue

        print(f"Computing heatmap for group {grp} ({len(sub_mf)} seeds)...")
        tmp_manifest = manifest_csv.replace(".csv", f"_tmp_{grp}.csv")
        sub_mf.to_csv(tmp_manifest, index=False)

        H_g, _, _, E95_g = aggregate_distance_scale_heatmap(
            manifest_csv=tmp_manifest,
            x_max_um=x_max_um,
            bin_um=bin_um,
            intensity_cut=intensity_cut,
            normalize_by_self=normalize_by_self,
            normalize_rows=normalize_rows,
            use_counts_weight=use_counts_weight,
            smooth_sigma_bins=smooth_sigma_bins,
            energy_percentile=energy_percentile
        )

        H_groups[grp] = H_g
        E95_groups[grp] = E95_g

        if os.path.exists(tmp_manifest):
            os.remove(tmp_manifest)

    # ------------------ 3) Global color range ------------------
    all_H = [H_all] + [H_groups[g] for g in groups_order if H_groups[g] is not None]
    vmin_global = 0.0
    vmax_data = max([np.nanmax(H[np.isfinite(H)]) for H in all_H if H is not None])

    # symlog: clip true color mapping maximum to colorbar_vmax_symlog
    if colorbar_scale == "symlog":
        vmax_global = min(vmax_data, float(colorbar_vmax_symlog))
    else:
        vmax_global = vmax_data

    # ------------------ 4) Grid ------------------
    scales = np.asarray(scales, dtype=np.float32)
    dist_centers = np.asarray(dist_centers, dtype=np.float32)

    dx = dist_centers[1] - dist_centers[0] if dist_centers.size > 1 else 1.0
    x_edges = np.concatenate([[dist_centers[0] - dx / 2], dist_centers + dx / 2]).astype(np.float32)
    y_edges = _scale_edges(scales)

    # ------------------ 5) Subplots ------------------
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    # ------------------ 6) norm ------------------
    if colorbar_scale == "log":
        from matplotlib.colors import LogNorm
        norm = LogNorm(vmin=colorbar_vmin_log, vmax=vmax_global)
    elif colorbar_scale == "symlog":
        from matplotlib.colors import SymLogNorm
        norm = SymLogNorm(
            linthresh=0.01,
            linscale=1.0,
            vmin=vmin_global,
            vmax=vmax_global,
            base=10,
            clip=True
        )
    elif colorbar_scale == "power":
        from matplotlib.colors import PowerNorm
        norm = PowerNorm(gamma=colorbar_power, vmin=colorbar_vmin_log, vmax=vmax_global)
    else:
        norm = Normalize(vmin=vmin_global, vmax=vmax_global)

    # ------------------ 7) Plot heatmaps ------------------
    titles = [f"All"] + [f"{g}" for g in groups_order]
    H_list = [H_all] + [H_groups[g] for g in groups_order]
    E95_list = [E95_all] + [E95_groups[g] for g in groups_order]

    mesh_for_cbar = None

    for idx, (ax, H, E95, title) in enumerate(zip(axes, H_list, E95_list, titles)):
        if H is None:
            ax.text(0.5, 0.5, f"No data\nfor {title}",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        H_plot = H.copy()

        if colorbar_scale in ("log", "log2"):
            vmin_val = float(colorbar_vmin_log)
            H_plot = np.where((H_plot < vmin_val) & np.isfinite(H_plot), vmin_val, H_plot)
            H_plot = np.nan_to_num(H_plot, nan=vmin_val)

        mesh = ax.pcolormesh(x_edges, y_edges, H_plot, shading="auto", cmap=cmap, norm=norm)
        if mesh_for_cbar is None:
            mesh_for_cbar = mesh

        ax.set_xlim(0, float(x_max_um))
        if y_max_scale is not None:
            ax.set_ylim(0, float(y_max_scale))

        if idx >= 2:
            ax.set_xlabel(xlabel, fontsize=6)
        if idx % 2 == 0:
            ax.set_ylabel(ylabel, fontsize=6)

        ax.set_title(title, fontsize=7)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.tick_params(axis="both", which="major", labelsize=6)

        # 95% energy line + important scale markers
        if E95 is not None:
            E95 = np.asarray(E95, dtype=np.float32)

                # Compatibility: check whether E95 already includes scale=0
            if E95.shape[0] == scales.shape[0]:
                E95_x = E95
                E95_y = scales
            elif E95.shape[0] == scales.shape[0] - 1:
                E95_x = np.concatenate([[0.0], E95]).astype(np.float32)
                E95_y = scales
            else:
                E95_x = None
                E95_y = None

            if E95_x is not None:
                ax.plot(E95_x, E95_y, color=E95_color, linewidth=1.5, alpha=0.85, label="95% energy")

                if important_scales is not None:
                    try:
                        imp = np.asarray(list(important_scales), dtype=np.float32)
                    except Exception:
                        imp = np.asarray(important_scales, dtype=np.float32)

                        # Automatically add scale 0
                    if imp.size == 0:
                        imp = np.array([0.0], dtype=np.float32)
                    else:
                        if not np.any(np.isclose(imp, 0.0, atol=float(1e-4))):
                            imp = np.concatenate([np.array([0.0], dtype=np.float32), imp])
                        imp = np.unique(imp)

                    for s in imp:
                        idx_imp = int(np.nanargmin(np.abs(E95_y - float(s))))
                        sy = float(E95_y[idx_imp])
                        sx = float(E95_x[idx_imp])

                            # Only mark scales that exactly exist (strict match)
                        if np.abs(sy - float(s)) > float(1e-4):
                            continue

                        ax.scatter(
                            [sx], [sy],
                            s=12,
                            c=E95_color,
                            edgecolors=E95_color,
                            linewidths=0.5,
                            zorder=5
                        )

            leg = ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.85),
                            frameon=True, fontsize=6)
            for text in leg.get_texts():
                text.set_color("white")
            frame = leg.get_frame()
            frame.set_facecolor("white")
            frame.set_alpha(0.4)
            frame.set_edgecolor("none")

    # ------------------ 8) colorbar ------------------
    if mesh_for_cbar is None:
        raise RuntimeError("No valid heatmap mesh was created; cannot draw colorbar.")

    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])

    extend_opt = "max" if colorbar_scale == "symlog" else "neither"
    cbar = fig.colorbar(mesh_for_cbar, cax=cbar_ax, extend=extend_opt)
    cbar.set_label("Mean intensity", rotation=270, labelpad=5, fontsize=6)
    cbar.ax.tick_params(labelsize=6)

    # ------------------ 9) Colorbar format ------------------
    if colorbar_scale == "log":
        cbar.ax.set_ylabel("Mean intensity (log₁₀)", rotation=270, labelpad=5, fontsize=6)

    elif colorbar_scale == "log2":
        from matplotlib.ticker import FixedLocator
        def log2_to_original(x, pos):
            return f"{2**x:.1e}"
        cbar.ax.yaxis.set_major_formatter(FuncFormatter(log2_to_original))
        tick_values = [1e-5, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0]
        tick_positions = [np.log2(v) for v in tick_values]
        cbar.ax.yaxis.set_major_locator(FixedLocator(tick_positions))
        cbar.ax.set_ylabel("Mean intensity (log₂)", rotation=270, labelpad=5, fontsize=6)

    elif colorbar_scale == "symlog":
        clip_val = float(vmax_global)  # = min(vmax_data, colorbar_vmax_symlog)

        # Only display up to clip_val (typically 0.5)
        tick_values = [0.0, 1e-2, 0.02, 0.04, 0.06, 0.08, 1e-1, 0.2, 0.3, 0.4, 0.5]
        tick_values = [t for t in tick_values if t >= vmin_global and t <= clip_val]
        if len(tick_values) == 0 or not np.isclose(tick_values[-1], clip_val):
            tick_values.append(clip_val)
        tick_values = sorted(set([float(t) for t in tick_values]))
        cbar.set_ticks(tick_values)

        def symlog_label(v, pos):
            if np.isclose(v, 0.0):
                return "0"
            return f"{v:.0e}" if v < 0.01 else f"{v:g}"

        cbar.ax.yaxis.set_major_formatter(FuncFormatter(symlog_label))
        cbar.ax.set_ylabel(
            f"Mean intensity (symlog)",
            rotation=270, labelpad=10, fontsize=6
        )

    # ------------------ 10) Overall title ------------------
    if suptitle is None:
        suptitle = f"Scale-distance heatmaps by {group_column}"
    fig.suptitle(suptitle, fontsize=7, y=0.98)

    plt.tight_layout(rect=[0, 0, 0.92, 0.96])

    # ------------------ 11) Collect results and save data ------------------
    results = {}

    base_out = out_dir if out_dir is not None else os.path.dirname(manifest_csv)
    data_dir = os.path.join(base_out, "heatmap_by_groups")
    if save_data:
        os.makedirs(data_dir, exist_ok=True)

    results["All"] = {
        "H": H_all,
        "E95": E95_all,
        "scales": scales,
        "dist_centers": dist_centers,
        "n_seeds": len(mf)
    }

    H_all_df = pd.DataFrame(
        H_all,
        index=[f"scale_{s:.2f}" for s in scales],
        columns=[f"dist_{d:.2f}" for d in dist_centers]
    )
    results["All"]["H_df"] = H_all_df

    E95_all_df = pd.DataFrame({"scale": scales, "E95_distance_um": E95_all})
    results["All"]["E95_df"] = E95_all_df

    if save_data:
        H_all_df.to_csv(os.path.join(data_dir, "H_matrix_All.csv"))
        E95_all_df.to_csv(os.path.join(data_dir, "E95_distance_All.csv"), index=False)
        print(f"[Saved] All group data -> {data_dir}/H_matrix_All.csv, E95_distance_All.csv")

    for grp in groups_order:
        H_g = H_groups.get(grp, None)
        E95_g = E95_groups.get(grp, None)

        results[grp] = {
            "H": H_g,
            "E95": E95_g,
            "scales": scales,
            "dist_centers": dist_centers,
            "n_seeds": n_seeds_groups.get(grp, 0)
        }

        if H_g is not None:
            H_g_df = pd.DataFrame(
                H_g,
                index=[f"scale_{s:.2f}" for s in scales],
                columns=[f"dist_{d:.2f}" for d in dist_centers]
            )
            results[grp]["H_df"] = H_g_df

            E95_g_df = pd.DataFrame({"scale": scales, "E95_distance_um": E95_g})
            results[grp]["E95_df"] = E95_g_df

            if save_data:
                H_g_df.to_csv(os.path.join(data_dir, f"H_matrix_{grp}.csv"))
                E95_g_df.to_csv(os.path.join(data_dir, f"E95_distance_{grp}.csv"), index=False)
                print(f"[Saved] {grp} group data -> {data_dir}/H_matrix_{grp}.csv, E95_distance_{grp}.csv")
        else:
            results[grp]["H_df"] = None
            results[grp]["E95_df"] = None

    # ------------------ 12) Save PDF ------------------
    if save_pdf:
        os.makedirs(data_dir, exist_ok=True)
        pdf_path = os.path.join(out_dir, f"Heatmap_by_{group_column}.pdf")
        fig.savefig(pdf_path, dpi=dpi, bbox_inches="tight")
        print(f"[Saved] PDF -> {pdf_path}")

    return fig, axes, results



import os
import re
import time
import numpy as np
import matplotlib.pyplot as plt
from pygsp import graphs, filters

def extract_typical_scales_from_obsm(adata, prefix="X_umap_scale"):
    """Extract scales from adata.obsm key names, e.g. X_umap_scale8.0 -> 8.0"""
    scales = []
    for k in adata.obsm.keys():
        if k.startswith(prefix):
            m = re.search(rf"{re.escape(prefix)}(.+)", k)
            if m is not None:
                scales.append(float(m.group(1)))
    return sorted(set(scales))


def mapping_scale_to_distance(
    adata,
    out_dir,
    run_label="sample",

    # ---- typical_scales ----
    typical_scales=None,                # None -> auto-extract from obsm; or provide a manual list
    plot_typical_scales=True,           # whether to mark typical_scales on the heatmap

    # ---- Density ----
    density_radius=20,
    density_key="Density",
    density_percentile=10,              # 10 -> q10/q90
    density_group_key="Density_group",

    # ---- Sampling ----
    sample_r=500,
    target_total=2000,
    first_round_target_n=400,
    first_round_candidate_cap=10000,
    preferred_max_frac=0.6,
    sampling_obs_key="Sample_round",
    sampled_flag_key="Sampled",
    plot_sampling_rounds=True,

    # ---- Graph heat kernel scales (filter scales) ----
    candidate_scales=None,                 # None -> use the default scale grid

    # ---- Seed export ----
    export_radius_um=250.0,
    n_jobs=8,

    # ---- Heatmap ----
    group_column="Density_group",
    groups_order=("Low", "Medium", "High"),
    x_max_um=125.0,
    bin_um=1.0,
    normalize_rows=True,
    colorbar_scale="symlog",
    cmap="turbo",
    heatmap_figsize=(8, 6),
    colorbar_vmax_symlog=0.5,

    # ---- Other ----
    show=True,
    dpi=300,
):
    """
    Run the complete scale-to-distance mapping pipeline:
      1) Compute cell density + group by percentile
      2) Use q_low/q_high extreme-density cells as preferred sampling candidates
      3) Multi-round spatial sampling (writes to adata.obs[sampling_obs_key])
      4) Build graph heat-kernel filter bank, produce transformed signals
      5) export_all_seeds -> manifest.csv
      6) plot_distance_scale_heatmap_by_groups (optionally mark typical_scales)

    Results are saved to a 'scale_to_distance_mapping' subfolder inside out_dir.

    Returns
    -------
    adata : AnnData
        Input AnnData with density and sampling columns added.
    """

    out_dir = os.path.join(out_dir, "scale_to_distance_mapping")
    os.makedirs(out_dir, exist_ok=True)
    start_time = time.time()

    # --------- typical_scales ---------
    if typical_scales is None:
        typical_scales = extract_typical_scales_from_obsm(adata, prefix="X_umap_scale")
    if (not plot_typical_scales) or (typical_scales is None) or (len(typical_scales) == 0):
        important_scales = None
    else:
        important_scales = list(typical_scales)

    # --------- candidate_scales ---------
    if candidate_scales is None:
        scales = ([0.001]
                  + np.arange(0.1, 2.1, 0.1).tolist()
                  + np.arange(2.5, 15.5, 0.5).tolist()
                  + np.arange(16, 21, 1).tolist()
                  + np.arange(25, 55, 5).tolist())
        candidate_scales = [round(x, 2) for x in scales]
    else:
        candidate_scales = [float(x) for x in candidate_scales]

    # --------- 1) density ---------
    density = calculate_density_matrices(adata, radius=density_radius)
    adata.obs[density_key] = density.copy()

    q_low, q_high = classify_density_by_quantile(
        adata,
        density_key=density_key,
        percentile=density_percentile,
        new_key=density_group_key
    )

    # preferred cells: high-density + low-density extremes
    high_density_cells = adata.obs_names[adata.obs[density_key] >= q_high]
    low_density_cells = adata.obs_names[adata.obs[density_key] <= q_low]
    preferred_cells = list(high_density_cells) + list(low_density_cells)

    # --------- 2) graph + filter ---------
    if "spatial_knn" not in adata.obsp:
        raise KeyError("'spatial_knn' not found in adata.obsp (required for graph construction)")

    spatial_knn = adata.obsp["spatial_knn"].copy()
    G_pygsp = graphs.Graph(spatial_knn)
    G_pygsp.estimate_lmax()

    g = filters.Heat(G_pygsp, candidate_scales)

    # --------- 3) sampling ---------
    round2ids = multi_round_sampling_ultrafast(
        adata,
        r=float(sample_r),
        target_total=int(target_total),
        base_seed=0,  # fixed seed for reproducibility
        key_spatial="spatial",
        obs_key_round=sampling_obs_key,
        obs_key_flag=sampled_flag_key,
        first_round_target_n=int(first_round_target_n),
        first_round_candidate_cap=int(first_round_candidate_cap) if first_round_candidate_cap is not None else None,
        preferred_ids=preferred_cells,
        preferred_max_frac=float(preferred_max_frac),
        verbose=False
    )

    # Visualize sampling rounds
    sampling_fig_path = None
    if plot_sampling_rounds:
        fig, ax = plot_spatial_sampling_rounds(
            adata,
            obs_key=sampling_obs_key,
            s_bg=1,
            s_fg=5,
            figsize=(10, 4)
        )
        sampling_fig_path = os.path.join(out_dir, f"{run_label}_sampling_rounds.pdf")
        fig.savefig(sampling_fig_path, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()
        plt.close(fig)

    # --------- 4) make signal + filter ---------
    raw_signals, rounds = make_signal_matrix(adata, obs_key=sampling_obs_key)   # (N, R)
    transformed = g.filter(raw_signals)                                         # (N, R, S)

    # --------- 5) export seeds ---------
    seed_dir = os.path.join(out_dir, "seed_profiles")
    manifest_csv = export_all_seeds(
        adata=adata,
        transformed=transformed,
        scales=candidate_scales,
        out_dir=seed_dir,
        obs_key=sampling_obs_key,
        radius_um=float(export_radius_um),
        n_jobs=int(n_jobs)
    )

    # --------- 6) heatmap by groups ---------
    fig_hm, axes_hm, results_hm = plot_distance_scale_heatmap_by_groups(
        manifest_csv=manifest_csv,
        group_column=group_column,
        groups_order=groups_order,
        x_max_um=float(x_max_um),
        bin_um=float(bin_um),
        colorbar_scale=colorbar_scale,
        cmap=cmap,
        figsize=heatmap_figsize,
        important_scales=important_scales,
        normalize_rows=bool(normalize_rows),
        colorbar_vmax_symlog = colorbar_vmax_symlog, 
        out_dir = out_dir,
        dpi=dpi
    )
    if show:
        plt.show()
    plt.close(fig_hm)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken to perform physical distance analysis: {elapsed_time:.4f} seconds")

    return adata


# Backward-compatible alias
run_physical_distance_pipeline = mapping_scale_to_distance
