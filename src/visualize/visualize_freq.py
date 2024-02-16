import pandas as pd
import matplotlib.pyplot as plt

def visualize_freq(output_dir):
	df = pd.read_csv(f"{output_dir}/tsv/target_frequency.tsv", sep='\t', comment='#')
	df.set_index('index', inplace=True)

	bar_width = 0.8

	bar_colors = ['#7b9acc', '#9CC3D5']

	plt.figure(figsize=(40, 20))
	ax = df.plot(kind='bar', legend=False, width=bar_width)

	for i, bar in enumerate(ax.patches):
			bar.set_color(bar_colors[i % len(bar_colors)])

	plt.title('Total frequency of target gene in Lactobacillales', pad=20)

	ax.grid(axis='y', linestyle='-', alpha=0.3)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	ax.set_xlabel('')
	ax.set_ylabel('')
	ax.set_ylim()

	plt.savefig(f"{output_dir}/img/target_frequency.png")