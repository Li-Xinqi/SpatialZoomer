from memory_profiler import memory_usage
import time
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import pandas as pd
import numpy as np
import os
import STAGATE_pyG
import torch

def run_stagate(
        adata, 
        data_name,
        resolution=0.5,
        save_path="./saves/"
        ):
    print("Running STAGATE on "+data_name+"...")
    
    start_time = time.time()

    # clustering method -- STAGATE
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=50)
    STAGATE_pyG.Stats_Spatial_Net(adata)
    adata = STAGATE_pyG.train_STAGATE(adata, device = torch.device('cpu'))
    sc.pp.neighbors(adata, use_rep='STAGATE')
    for res in [0.2, 0.5, 1.0, 1.5, 2.8, 3.0, 3.5, 4.0]:
        sc.tl.louvain(adata, resolution=res)
        adata.obs['louvain_'+str(res)] = adata.obs['louvain']
    
    cluster_res_df = adata.obs

    # get total time consumption
    end_time = time.time()

    # save the results
    cluster_res_df.to_csv(save_path+data_name+"_stagate.csv")
    pd.DataFrame(
        {'time_consumption(s)': [end_time - start_time]}
    ).to_csv(save_path+data_name+"_time_stagate.csv", index=False)

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
    run_stagate( 
        ad, 
        data_name,
        resolution=0.5,
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
    ).to_csv(args.save_path + args.data_name+"_memory_stagate.csv", index=False)
    print(f"Peak memory usage: {max_mem:.2f} MiB")