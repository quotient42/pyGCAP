import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.colors import LinearSegmentedColormap

def visualize_by_genus(output_dir):
    main_df = pd.read_csv(f"{output_dir}/tsv/contents_genus.tsv", sep='\t', comment='#')

    main_df = main_df.fillna(0)

    heatmap_data = main_df.iloc[:, 1:]
    heatmap_data[heatmap_data == 0] = np.nan

    cmap = LinearSegmentedColormap.from_list('custom_cmap', ['#FFFFFF', '#999EBC', '#666E9B', '#4D568A', '#333D79'], N=256)

    plt.figure(figsize=(6, 18))
    plt.imshow(heatmap_data.values, cmap=cmap, interpolation='nearest', aspect='auto', vmin=0, vmax=4)
    plt.xticks(np.arange(len(heatmap_data.columns)), heatmap_data.columns, rotation='vertical', ha='center')
    plt.tick_params(axis='x', bottom=False, top=True, labelbottom=False, labeltop=True)

    italic_font = FontProperties()
    italic_font.set_style('italic')
    plt.yticks(np.arange(len(heatmap_data)), main_df['genus'], fontproperties=italic_font)
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=True, labelright=False)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/img/target_content_genus.png')
    plt.close()