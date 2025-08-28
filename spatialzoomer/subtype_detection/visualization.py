import numpy as np
import matplotlib.pyplot as plt

def generate_colors(num_colors):
    cmap_list = [plt.cm.tab20, plt.cm.tab20b, plt.cm.tab20c, plt.cm.tab10, plt.cm.Paired, plt.cm.Set3]
    colors = np.vstack([cmap(np.linspace(0, 1, cmap.N))[:, :3] for cmap in cmap_list])
    if len(colors) < num_colors:
        additional_colors = np.random.rand(num_colors - len(colors), 3)
        colors = np.vstack([colors, additional_colors])
    return colors[:num_colors]