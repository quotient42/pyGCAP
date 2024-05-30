import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties
from matplotlib.colors import LinearSegmentedColormap

#------------------------------------------------------------
def visualize_cluster_type(output_dir, df):
  x = np.arange(len(df))
  width = 0.3

  fig, ax = plt.subplots(figsize=(42, 6))

  bar1 = ax.bar(x - width, df['c_ratio'], width, label='Conserved', color='#3a7ca5')
  bar2 = ax.bar(x, df['f_ratio'], width, label='Fragmented', color='#F1A208')
  bar3 = ax.bar(x + width, df['d_ratio'], width, label='Disrupted', color='#588157')

  ax.set_xticks(x)
  ax.grid(axis='y', linestyle='-', alpha=0.3)
  ax.set_xticklabels(df['genus'], rotation=90, ha='right')
  ax.set_title('cluster type by genus', pad=20)

  ax.legend()

  plt.tight_layout()
  plt.savefig(f"{output_dir}/img/cluster_type.png")
  plt.close(fig)

#------------------------------------------------------------
def collect_gene_data(main_df):
    bar_data = {}

    for index, row in main_df.iterrows():
        bar_id = row['accession']

        if bar_id not in bar_data:
            bar_data[bar_id] = []

        gene_info = {
            'Prediction': row['Prediction'],
            'block': row['block'],
            'start': row['start'],
            'end': row['end'],
        }
        bar_data[bar_id].append(gene_info)

    return bar_data

def plot_subplot(ax, gene_data, bar_id, color_dict):
    handles = []

    for gene in gene_data:
        if not gene['Prediction'] in color_dict:
            gene['Prediction'] = 'etc'

        start = int(gene['start'])
        end = int(gene['end'])
        block = gene['block']

        bar = ax.barh(0, width=end - start, left=start, height=0.2, color=color_dict[gene['Prediction']], edgecolor='none')

    ax.set_yticks([0])
    ax.set_yticklabels([bar_id])

    ax.set_xlim(0, 500)
    ax.set_xticks([])

    return handles

def draw_cluster(cluster_data, genus, color_dict):
  fig, axes = plt.subplots(len(cluster_data) + 1, 1,
               figsize=(10, 0.5 * (len(cluster_data) + 1)),
               gridspec_kw={'height_ratios': [0.5] + [0.5] * len(cluster_data), 'hspace': 1})

  ax_legend = axes[0]
  ax_legend.axis('off')

  all_handles = []
  for idx, (bar_id, gene_data) in enumerate(cluster_data.items()):
    ax1 = axes[idx + 1]
    handles = plot_subplot(ax1, gene_data, bar_id, color_dict)
    all_handles.extend(handles)

    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)

  ax_legend.set_frame_on(False)

  legend_patches = [mpatches.Patch(color=color, label=gene) for gene, color in color_dict.items()]
  legend = ax_legend.legend(handles=legend_patches,
                               loc='center',
                ncol=7,
                edgecolor='white',
                handlelength=0.7,
                bbox_to_anchor=(0.5, 0.9),
                bbox_transform=plt.gcf().transFigure)

  axes[0].set_ylim(0, 1)
  legend.set_title(f"{genus}")

  return fig

def visualize_cluster(input_path, output_path, genus, color_dict):
  collected_df = pd.DataFrame()

  files = glob.glob(os.path.join(input_path, 'GCF_*.tsv'))

  for file_path in files:
    main_df = pd.read_csv(file_path, sep='\t')
    accession = os.path.basename(file_path).split('.tsv')[0]
    main_df['accession'] = accession
    collected_df = collected_df._append(main_df, ignore_index=True)

  collected_df = collected_df[['accession', 'block', 'Prediction', 'start', 'end']]

  bar_data = collect_gene_data(collected_df)
  fig = draw_cluster(bar_data, genus, color_dict)

  bar_img = f"{output_path}/{genus}.png"
  fig.savefig(bar_img, bbox_inches='tight')
  plt.close(fig)

#------------------------------------------------------------
def visualize_freq(output_dir, TAXON):
  df = pd.read_csv(f"{output_dir}/tsv/target_frequency.tsv", sep='\t', comment='#')
  df.set_index('index', inplace=True)

  bar_width = 0.8

  bar_colors = ['#7b9acc', '#9CC3D5']

  fig = plt.figure(figsize=(40, 20))
  ax = df.plot(kind='bar', legend=False, width=bar_width)

  for i, bar in enumerate(ax.patches):
      bar.set_color(bar_colors[i % len(bar_colors)])

  plt.title(f'Total frequency of target gene in {TAXON}', pad=20)

  ax.grid(axis='y', linestyle='-', alpha=0.3)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)

  ax.set_xlabel('')
  ax.set_ylabel('')
  ax.set_ylim()

  plt.savefig(f"{output_dir}/img/target_frequency.png")
  plt.close(fig)

#------------------------------------------------------------
def visualize_heatmap(output_dir, filename, mode):
  main_df = pd.read_csv(f"{output_dir}/tsv/{filename}.tsv", sep='\t', comment='#')

  main_df = main_df.fillna(0)
  
  if mode == 0 or mode == 2:
    main_df = main_df.sort_values(by='name')
    heatmap_data = main_df.iloc[:, 2:]
  elif mode == 1:
    main_df = main_df.sort_values(by='genus')
    heatmap_data = main_df.iloc[:, 1:]

  heatmap_data = heatmap_data.copy()
  heatmap_data[heatmap_data == 0] = np.nan

  cmap = LinearSegmentedColormap.from_list('custom_cmap', ['#FFFFFF', '#999EBC', '#666E9B', '#4D568A', '#333D79'], N=256)

  num_species = len(heatmap_data)
  num_features = len(heatmap_data.columns)

  if mode == 0 or mode == 1:
    figsize_width = min(15, num_features * 0.4)
    figsize_height = max(20, num_species * 0.2)
  elif mode == 2:
    figsize_width = max(15, num_features * 0.3)
    figsize_height = max(20, num_species * 0.2)
  
  plt.figure(figsize=(figsize_width, figsize_height))
  plt.imshow(heatmap_data.values, cmap=cmap, interpolation='nearest', aspect='auto', vmin=0, vmax=1)

  plt.xticks(np.arange(num_features), heatmap_data.columns, rotation='vertical', ha='center')
  plt.tick_params(axis='x', bottom=False, top=True, labelbottom=False, labeltop=True)

  italic_font = FontProperties()
  italic_font.set_style('italic')

  if mode == 0 or mode == 2:
    plt.yticks(np.arange(num_species), [name.split()[0] + ' ' + name.split()[1] for name in main_df['name']],
              fontproperties=italic_font)
  elif mode == 1:
    plt.yticks(np.arange(len(heatmap_data)), main_df['genus'], fontproperties=italic_font)

  plt.tick_params(axis='y', which='both', left=False, right=False,
                  labelleft=True, labelright=False)
  plt.tight_layout()

  plt.savefig(f'{output_dir}/img/{filename}.png')
  plt.close()
