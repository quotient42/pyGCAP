import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.colors import LinearSegmentedColormap

def visualize_by_species(output_dir):
  main_df = pd.read_csv(f"{output_dir}/tsv/contents_species.tsv", sep='\t', comment='#')

  main_df = main_df.fillna(0)
  main_df = main_df.sort_values(by='name')

  heatmap_data = main_df.iloc[:, 2:]
  heatmap_data[heatmap_data == 0] = np.nan

  cmap = LinearSegmentedColormap.from_list('custom_cmap', ['#FFFFFF', '#92A2A0', '#5C7371', '#415C59', '#254441'], N=256)

  num_species = len(heatmap_data)
  num_features = len(heatmap_data.columns)

  plt.figure(figsize=(7, 120))
  plt.imshow(heatmap_data.values, cmap=cmap, interpolation='nearest', aspect='auto', vmin=0, vmax=4)

  plt.xticks(np.arange(num_features), heatmap_data.columns, rotation='vertical', ha='center')
  plt.tick_params(axis='x', bottom=False, top=True, labelbottom=False, labeltop=True)

  italic_font = FontProperties()
  italic_font.set_style('italic')
  plt.yticks(np.arange(num_species), [name.split()[0] + ' ' + name.split()[1] for name in main_df['name']],
            fontproperties=italic_font)
  plt.tick_params(axis='y', which='both', left=False, right=False,
                  labelleft=True, labelright=False)
  plt.tight_layout()

  plt.savefig(f'{output_dir}/img/contents_species.png')
  plt.close()
