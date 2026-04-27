from memory_profiler import memory_usage
import time
import warnings
warnings.filterwarnings("ignore")

import MENDER
import scanpy as sc
import pandas as pd
import os

def run_mender(
        ad, 
        data_name,
        spatial_key='spatial_physical',
        cluster_resolution=1.5,
        save_path="./saves/"
        ):
    ad.obsm['spatial'] = ad.obsm[spatial_key]
    print("Running MENDER on "+data_name+"...")
    
    start_time = time.time()

    # clustering method -- MENDER
    ad.var_names_make_unique()
    ad.layers["counts"] = ad.X.copy()
    sc.pp.highly_variable_genes(ad, flavor="seurat_v3", n_top_genes=4000)
    sc.pp.normalize_total(ad)
    sc.pp.pca(ad)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)
    sc.tl.leiden(ad,key_added='ct',resolution=2)

    # input parameters of MENDER
    scale = 6

    # the default radius is 15 um
    radius = 15

    # main body of MENDER
    msm = MENDER.MENDER_single(
        ad,
        # determine which cell state to use
        # we use the cell state got by Leiden
        ct_obs='ct'
    )
    msm.estimate_radius()

    msm.set_MENDER_para(
        # default of n_scales is 6
        n_scales=scale,

        # for single cell data, nn_mode is set to 'radius'
        nn_mode='radius',

        # default of n_scales is 15 um (see the manuscript for why).
        # MENDER also provide a function 'estimate_radius' for estimating the radius
        # this para in 3D might be smaller than 2D
        nn_para=radius,

    )
    # construct the context representation
    msm.run_representation(
        # the number of processings
    )

    # set the spatial clustering parameter
    # positive values for the expected number of domains
    # negative values for the clustering resolution
    for res in [0.2, 0.5, 1.0, 1.5, 1.8, 2.0]:
        msm.run_clustering_normal(-res)
        msm.output_cluster('MENDER')
        msm.adata_MENDER.obs['MENDER'+str(res)] = msm.adata_MENDER.obs['MENDER']
    cluster_res_df = msm.adata_MENDER.obs

    # get total time consumption
    end_time = time.time()


    # save the results
    cluster_res_df.to_csv(save_path+data_name+"_mender.csv")
    pd.DataFrame(
        {'time_consumption(s)': [end_time - start_time]}
    ).to_csv(save_path+data_name+"_time_mender.csv", index=False)

    print("time consumption {:.2f} seconds".format(end_time - start_time))
    return

def read_data_and_run(
        data_path, 
        data_name, 
        save_path="./saves/", 
        n_subsample=None
        ):
    try:
        ad = sc.read_h5ad(data_path + data_name + ".h5ad")
        if n_subsample is not None:
            sc.pp.subsample(ad, n_obs=n_subsample)
        run_mender(
            ad, 
            data_name,
            spatial_key='spatial',
            cluster_resolution=1.0,
            save_path=save_path
        )
    except Exception as e:
        print(f"Error processing {data_name}: {e}")
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
    ).to_csv(args.save_path + args.data_name+"_memory_mender.csv", index=False)
    print(f"Peak memory usage: {max_mem:.2f} MiB")