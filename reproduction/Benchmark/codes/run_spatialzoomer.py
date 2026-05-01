# import tracemalloc
from memory_profiler import memory_usage
import time
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import os
import spatialzoomer as sz

def run_spatialzoomer(
        adata, 
        data_name,
        save_path="./saves/"
        ):
    print("Running SpatialZoomer on "+data_name+"...")
    
    # tracemalloc.start()
    start_time = time.time()

    # clustering method -- SpatialZoomer
    adata = sz.Preprocess(adata)
    adata = sz.performDR(adata, type = 'NMF', n_components=50)

    runLabel = data_name
    sz_analyzer = sz.MultiscaleAnalysis(adata=adata, runLabel=runLabel, save_path='/home/glab/lxq/SpatialZoomer/Final_data/benchmark/saves/run_singleres_all/#./saves/spatialzoomer/trash/')
    sz_analyzer.multiscale_transform(use_rep='X_nmf', n_neighbors=20)    # scales: use default
    sz_analyzer.identify_typical_scales(max_clusters=10, min_clusters=3, show = False)
    resolutions = [0.1, 0.4, 0.6, 0.8, 1.0, 2.0]
    sz_analyzer.clustering(
        n_clusters_kmeans=min(10000, adata.n_obs), 
        resolutions=resolutions, 
        max_scale=None)

    scales = sz_analyzer.scales_df_use['Scale']

    cluster_res_df = pd.DataFrame(index=adata.obs_names)
    for res in resolutions:
        for scale in scales:
            cluster_key = 'leiden_scale'+str(scale)+'_res'+str(res)
            cluster_res_df[cluster_key] = sz_analyzer.adata.obs[cluster_key]

    # get total time consumption
    end_time = time.time()

    # get current and peak memory usage
    # current, peak = tracemalloc.get_traced_memory()
    # tracemalloc.stop()

    # save the results
    cluster_res_df.to_csv(save_path+data_name+"_spatialzoomer.csv")
    pd.DataFrame(
        {'time_consumption(s)': [end_time - start_time]}
    ).to_csv(save_path+data_name+"_time_spatialzoomer.csv", index=False)

    print("time consumption {:.2f} seconds".format(end_time - start_time))
    # print(f"peak memory usage: {peak / 1024 / 1024:.2f} MB")
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
    run_spatialzoomer(
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
    ).to_csv(args.save_path + args.data_name+"_memory_spatialzoomer.csv", index=False)
    print(f"Peak memory usage: {max_mem:.2f} MiB")