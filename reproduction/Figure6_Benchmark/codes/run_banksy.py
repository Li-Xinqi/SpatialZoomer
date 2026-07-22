from memory_profiler import memory_usage
import time
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import pandas as pd
import numpy as np
import os

import sys
sys.path.insert(0, '/home/glab/lxq/SpatialZoomer/Final_data/benchmark/Banksy_py/')
from banksy_utils.filter_utils import normalize_total, filter_hvg, print_max_min
from banksy.main import median_dist_to_nearest_neighbour
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.main import concatenate_all
from banksy_utils.umap_pca import pca_umap
from banksy.cluster_methods import run_Leiden_partition
from banksy.plot_banksy import plot_results

import random
# Note that BANKSY itself is deterministic, here the seeds affect the umap clusters and leiden partition
seed = 1234
np.random.seed(seed)
random.seed(seed)

def run_banksy(
        adata, 
        data_name,
        save_path="./saves/"
        ):
    print("Running Banksy on "+data_name+"...")
    
    start_time = time.time()

    # clustering method -- Banksy
    coord_keys = ('X', 'Y', 'spatial')
    adata.obs['X'] = adata.obsm['spatial'][:, 0]
    adata.obs['Y'] = adata.obsm['spatial'][:, 1]

    adata = normalize_total(adata)
    adata, adata_allgenes = filter_hvg(adata,
                n_top_genes = 2000,
                flavor="seurat")
    plot_graph_weights = True
    k_geom = 15 # only for fixed type
    max_m = 1 # azumithal transform up to kth order
    nbr_weight_decay = "scaled_gaussian" # can also be "reciprocal", "uniform" or "ranked"

    # Find median distance to closest neighbours, the median distance will be `sigma`
    nbrs = median_dist_to_nearest_neighbour(adata, key = coord_keys[2])
    banksy_dict = initialize_banksy(
        adata,
        coord_keys,
        k_geom,
        nbr_weight_decay=nbr_weight_decay,
        max_m=max_m,
        plt_edge_hist=False,
        plt_nbr_weights=False,
        plt_agf_angles=False, # takes long time to plot
        plt_theta=False,
    )

    # The following are the main hyperparameters for BANKSY
    resolutions = [0.2, 0.5, 0.7, 1.0, 1.2, 2.2, 2.5, 3.0] # clustering resolution for UMAP
    pca_dims = [20] # Dimensionality in which PCA reduces to
    lambda_list = [0.2] # list of lambda parameters

    banksy_dict, banksy_matrix = generate_banksy_matrix(adata, banksy_dict, lambda_list, max_m)

    banksy_dict["nonspatial"] = {
        # Here we simply append the nonspatial matrix (adata.X) to obtain the nonspatial clustering results
        0.0: {"adata": concatenate_all([adata.X], 0, adata=adata), }
    }

    pca_umap(banksy_dict,
            pca_dims = pca_dims,
            add_umap = True,
            plt_remaining_var = False,
            )
    results_df, max_num_labels = run_Leiden_partition(
        banksy_dict,
        resolutions,
        num_nn = 50,
        num_iterations = -1,
        partition_seed = seed,
        match_labels = True,
    )

    c_map =  'tab20' # specify color map
    weights_graph =  banksy_dict['scaled_gaussian']['weights'][0]

    plot_results(
        results_df,
        weights_graph,
        c_map,
        match_labels = True,
        coord_keys = coord_keys,
        max_num_labels  =  max_num_labels, 
        save_path = './banksy_xbp_brain_tmp_png',
        save_fig = False, # save the spatial map of all clusters
        save_seperate_fig = False, # save the figure of all clusters plotted seperately
    )
    cluster_res_df_list = []
    for idx in results_df.index:
        cluster_res_df_list.append(results_df.loc[idx, 'adata'].obs['labels_'+idx])
    combined_df = pd.concat(cluster_res_df_list, axis=1)
    # get total time consumption
    end_time = time.time()


    # save the results
    combined_df.to_csv(save_path+data_name+"_banksy.csv")
    pd.DataFrame(
        {'time_consumption(s)': [end_time - start_time]}
    ).to_csv(save_path+data_name+"_time_banksy.csv", index=False)

    print("time consumption {:.2f} seconds".format(end_time - start_time))
    return

def read_data_and_run(
        data_path, 
        data_name, 
        save_path="./saves/", 
        n_subsample=None
        ):
    ad = sc.read_h5ad(data_path + data_name + ".h5ad")
    if n_subsample is not None:
        sc.pp.subsample(ad, n_obs=n_subsample)
    run_banksy(
        ad, 
        data_name,
        save_path=save_path
    )
    # Force garbage collection
    import gc
    gc.collect()
    return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run on a single spatial dataset.')
    parser.add_argument('--data_path', type=str, default='./datasets/standard/', help='Path to dataset directory')
    parser.add_argument('--save_path', type=str, default='./saves/run_new_0714/', help='Path to save results')
    parser.add_argument('--data_name', type=str, required=True, help='Name of the dataset (without extension)')
    parser.add_argument('--n_subsample', type=int, default=None, help='Number of cells to subsample (optional)')
    
    args = parser.parse_args()
    
    # Ensure the save path exists
    os.makedirs(args.save_path, exist_ok=True)
    
    mem = memory_usage(
        (read_data_and_run,
        (args.data_path, args.data_name), 
        {'save_path': args.save_path,
        'n_subsample': None}),
        interval=1
        )
    max_mem = max(mem)
    pd.DataFrame(
        {'memory_usage(mb)':[max_mem]}
    ).to_csv(args.save_path + args.data_name+"_memory_banksy.csv", index=False)
    print(f"Peak memory usage: {max_mem:.2f} MiB")