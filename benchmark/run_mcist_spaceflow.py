# -*- coding: utf-8 -*-

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
import sys
from pathlib import Path

#from MCIST_GATE import MCIST_GATE
#from MCIST_GraphST import MCIST_GraphST
from MCIST_SpaceFlow import MCIST_SpaceFlow
import matplotlib.pyplot as plt
from mclustpy import mclustpy
from sklearn.metrics.cluster import normalized_mutual_info_score


def run_mcistspaceflow(
        adata, 
        data_name,
        save_path="./saves/",
        ):
    print("Running SpatialZoomer on "+data_name+"...")
    
    # tracemalloc.start()
    start_time = time.time()

    # clustering method -- MCIST_spaceflow
    n = 10
        
    adata = MCIST_SpaceFlow(adata = adata, n_clusters = n, clustering_algo='Mclust')
    cluster_res_df = adata.obs['MCIST_spatial_domains']
    

    # get total time consumption
    end_time = time.time()


    # save the results
    cluster_res_df.to_csv(save_path+data_name+f"_mcistspaceflow_{n}.csv")
    pd.DataFrame(
        {'time_consumption(s)': [end_time - start_time]}
    ).to_csv(save_path+data_name+"_time_mcistspaceflow.csv", index=False)

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
    run_mcistspaceflow(
        ad, 
        data_name,
        save_path=save_path
    )
    # ǿ�ƴ�����������
    import gc
    gc.collect()
    return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run on a single spatial dataset.')
    parser.add_argument('--data_path', type=str, default='./datasets/standard/', help='Path to dataset directory')
    parser.add_argument('--save_path', type=str, default='./saves/run_new_2602/', help='Path to save results')
    parser.add_argument('--data_name', type=str, required=True, help='Name of the dataset (without extension)')
    parser.add_argument('--n_subsample', type=int, default=None, help='Number of cells to subsample (optional)')
    
    args = parser.parse_args()
    
    # ȷ������·������
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
    ).to_csv(args.save_path + args.data_name+"_memory_mcistspaceflow.csv", index=False)
    print(f"Peak memory usage: {max_mem:.2f} MiB")