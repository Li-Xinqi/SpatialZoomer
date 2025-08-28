import matplotlib.pyplot as plt
import squidpy as sq
import scanpy as sc
import numpy as np

def filter_similar_colors(colors, threshold=0.1):
    np.random.seed(0)  
    filtered_colors = []
    for color in colors:
        if all(np.linalg.norm(color - fc) >= threshold for fc in filtered_colors):
            filtered_colors.append(color)
    return np.array(filtered_colors)


def generate_colors(num_colors):
    np.random.seed(0) 
    cmap_list = [plt.cm.tab20, plt.cm.tab20b, plt.cm.tab20c, plt.cm.tab10, plt.cm.Paired, plt.cm.Set3]
    colors = np.vstack([cmap(np.linspace(0, 1, cmap.N))[:, :3] for cmap in cmap_list])
    if len(colors) < num_colors:
        additional_colors = np.random.rand(num_colors - len(colors), 3)
        colors = np.vstack([colors, additional_colors])
        new_colors = filter_similar_colors(colors, threshold=0.1)
    return new_colors[:num_colors]

def plot_clusters(ad, typical_scales, resolution, figsize = None, width_ratios = None, pt_size_umap = 5, pt_size_scatter = 0.2):
    scales_plot = ['Raw'] + ['scale' + str(scale) for scale in typical_scales]
    new_colors = generate_colors(200)

    if figsize is None:
        figsize = (12, (len(typical_scales)+1) * 5)

    if width_ratios is None:
        width_ratios = [1, 2]

    fig, axes = plt.subplots(
        len(typical_scales)+1, 2,
        figsize=figsize,
        dpi=300,
        gridspec_kw={'width_ratios': width_ratios}
    )
    
    for i, scale in enumerate(scales_plot):
        cluster_key = 'leiden_' + str(scale) + '_res' + str(resolution)
        ad.obsm['X_umap'] = ad.obsm['X_umap_' + str(scale)]

        cell_type_counts = ad.obs[cluster_key].value_counts()
        sorted_cell_types = cell_type_counts.index.tolist()
        color_mapping = {cell_type: new_colors[jj] for jj, cell_type in enumerate(sorted_cell_types)}
        sorted_colors = [color_mapping[cell_type] for cell_type in ad.obs[cluster_key].cat.categories]
        ad.uns[cluster_key + '_colors'] = sorted_colors

        scale_clean = scale.replace("scale", "")
        sc.pl.umap(
            ad,  
            color=cluster_key, 
            size=pt_size_umap, 
            ax=axes[i, 0], 
            show=False, 
            title='UMAP of ' + str(scale)
        )
        for collection in axes[i, 0].collections:
            collection.set_rasterized(True)
        # hide legend
        axes[i, 0].legend_.remove()

        sq.pl.spatial_scatter(
            ad, 
            library_id="spatial", 
            shape=None, 
            color= cluster_key,
            wspace=0.1,
            frameon=False,
            size=pt_size_scatter,
            figsize=(18, 18),
            dpi=100,
            outline=False,
            img=False,
            marker='.',
            ax=axes[i, 1],
            title=f'Scale={scale_clean}, Resolution={resolution}',
        )
        for collection in axes[i, 1].collections:
            collection.set_rasterized(True)
    plt.tight_layout()
    plt.show()
    plt.close()
