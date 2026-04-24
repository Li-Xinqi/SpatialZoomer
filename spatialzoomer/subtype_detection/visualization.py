import matplotlib.pyplot as plt
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