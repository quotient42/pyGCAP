import matplotlib.pyplot as plt
import numpy as np

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