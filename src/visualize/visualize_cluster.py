import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import sys
sys.path.append('../info/')
from target_info import target_color

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.max_colwidth', None)

def collect_gene_data(main_df):
    bar_data = {}

    for index, row in main_df.iterrows():
        bar_id = row['accession']

        if bar_id not in bar_data:
            bar_data[bar_id] = []

        gene_info = {
            'prediction': row['prediction'],
            'block': row['block'],
            'start': row['start'],
            'end': row['end'],
        }
        bar_data[bar_id].append(gene_info)

    return bar_data

def plot_subplot(ax, gene_data, bar_id):
    color_dict = {}
    handles = []

    for gene in gene_data:
        if gene['prediction'] in target_color:
            color_dict[gene['prediction']] = target_color[gene['prediction']]
        else:
            color_dict[gene['prediction']] = target_color['etc']

        start = int(gene['start'])
        end = int(gene['end'])
        block = gene['block']

        bar = ax.barh(0, width=end - start, left=start, height=0.2, color=color_dict[gene['prediction']], edgecolor='none')

    ax.set_yticks([0])
    ax.set_yticklabels([bar_id])

    ax.set_xlim(0, 300)
    ax.set_xticks([])

    return handles


def draw_cluster(cluster_data, genus):
	fig, axes = plt.subplots(len(cluster_data) + 1, 1,
							 figsize=(10, 0.5 * (len(cluster_data) + 1)),
							 gridspec_kw={'height_ratios': [0.5] + [0.5] * len(cluster_data), 'hspace': 1})

	ax_legend = axes[0]
	ax_legend.axis('off')

	all_handles = []
	for idx, (bar_id, gene_data) in enumerate(cluster_data.items()):
		ax1 = axes[idx + 1]
		handles = plot_subplot(ax1, gene_data, bar_id)
		all_handles.extend(handles)

		ax1.spines['right'].set_visible(False)
		ax1.spines['left'].set_visible(False)
		ax1.spines['top'].set_visible(False)
		ax1.spines['bottom'].set_visible(False)

	ax_legend.set_frame_on(False)

	legend_patches = [mpatches.Patch(color=color, label=gene) for gene, color in target_color.items()]
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

def visualize_cluster(input_path, output_path, genus):
	collected_df = pd.DataFrame()

	files = glob.glob(os.path.join(input_path, 'GCF_*.tsv'))

	for file_path in files:
		main_df = pd.read_csv(file_path, sep='\t')
		accession = os.path.basename(file_path).split('.tsv')[0]
		main_df['accession'] = accession
		collected_df = collected_df._append(main_df, ignore_index=True)

	collected_df = collected_df[['accession', 'block', 'prediction', 'start', 'end']]

	bar_data = collect_gene_data(collected_df)
	fig = draw_cluster(bar_data, genus)

	bar_img = f"{output_path}/{genus}.png"
	fig.savefig(bar_img, bbox_inches='tight')
	plt.close(fig)
