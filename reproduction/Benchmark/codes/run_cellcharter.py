from memory_profiler import memory_usage
import time
import warnings
warnings.filterwarnings("ignore")

import cellcharter as cc
import scanpy as sc
import squidpy as sq
import scvi
import pandas as pd
import numpy as np
import os

def run_cellcharter(
        adata, 
        data_name,
        n_clusters=(20,50),
        save_path="./saves/"
        ):
    print("Running CellCharter on "+data_name+"...")
    
    start_time = time.time()

    # clustering method -- CellCharter
    adata.obs['sample'] = 'current_sample'
    adata.obs['sample'] = adata.obs['sample'].astype('category')
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    scvi.model.SCVI.setup_anndata(
        adata, 
        layer="counts",  
        batch_key='sample',
    )
    model = scvi.model.SCVI(adata)
    model.train(early_stopping=True, enable_progress_bar=True)
    adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
    sq.gr.spatial_neighbors(adata, library_key='sample', coord_type='generic', delaunay=True, spatial_key='spatial', percentile=99)
    cc.gr.aggregate_neighbors(adata, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter', sample_key='sample')
    for n_cluster in [10, 20, 30, 40, 50]:
        autok = cc.tl.ClusterAutoK(
            n_clusters=(n_cluster, n_cluster+10), 
            max_runs=10,
            convergence_tol=0.001
        )
        autok.fit(adata, use_rep='X_cellcharter')
        adata.obs['cluster_cellcharter'+str(n_cluster)] = autok.predict(adata, use_rep='X_cellcharter')

    cluster_res_df = adata.obs

    # get total time consumption
    end_time = time.time()

    # save the results
    cluster_res_df.to_csv(save_path+data_name+"_cellcharter.csv")
    pd.DataFrame(
        {'time_consumption(s)': [end_time - start_time]}
    ).to_csv(save_path+data_name+"_time_cellcharter.csv", index=False)

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
    run_cellcharter(
        ad, 
        data_name,
        n_clusters=(20, 50),
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
    ).to_csv(args.save_path + args.data_name+"_memory_cellcharter.csv", index=False)
    print(f"Peak memory usage: {max_mem:.2f} MiB")