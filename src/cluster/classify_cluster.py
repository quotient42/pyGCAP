import pandas as pd
import os
import time
import sys

sys.path.append("../visualize")
from visualize_cluster_type import visualize_cluster_type

def classify_target_cluster(project_info):
	print("<< classifying target gene cluster...")
	input_dir = project_info['input']
	output_dir = project_info['output']
	all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]

	start_time = time.time()

	cluster_df = pd.read_csv(f"{output_dir}/tsv/contents_species.tsv", sep='\t', comment='#').copy()
	cluster_df.insert(1, 'cluster', "-")

	for genus in all_directories:
			for accession in os.listdir(os.path.join(input_dir, genus)):
					if not "GCF" in accession:
							continue
					
					filename = f"{input_dir}/{genus}/{accession}/target_cluster.tsv"
					df = pd.read_csv(filename, sep='\t', comment='#')

					protein_cnt = len(df[(df['Prediction'] != '----') & (df['cluster'] == True)])
					df = df[df['block'] > 0]
					block_sizes = []
					for block_num in df['block'].unique():
							Predictions_in_block = df[df['block'] == block_num]['Prediction'].tolist()
							non_dash_count = sum(1 for pred in Predictions_in_block if pred != '----')
							block_sizes.extend([non_dash_count] * len(Predictions_in_block))

					df['block_size'] = block_sizes
					block_cnt = df['block_size'].max()

					if protein_cnt < 9:
						result = "D"
					elif block_cnt >= 9:
						result = "C"
					else:
						result = "F"
					
					cluster_df.loc[cluster_df['accession'] == accession, 'cluster'] = result

	split_names = cluster_df['name'].str.split(' ', n=1, expand=True)
	cluster_df['genus'] = split_names[0]
	cluster_df['species'] = split_names[1]

	cluster_df = cluster_df[['genus', 'species', 'cluster']]
	cluster_df.to_csv(f"{output_dir}/tsv/cluster_type_species.tsv", sep='\t', index=False)

	conserved_count = cluster_df[cluster_df['cluster'] == 'C'].groupby('genus').size()
	fragmented_count = cluster_df[cluster_df['cluster'] == 'F'].groupby('genus').size()
	disrupted_count = cluster_df[cluster_df['cluster'] == 'D'].groupby('genus').size()

	cluster_df['conserved'] = cluster_df['genus'].map(conserved_count).fillna(0).astype(int)
	cluster_df['fragmented'] = cluster_df['genus'].map(fragmented_count).fillna(0).astype(int)
	cluster_df['disrupted'] = cluster_df['genus'].map(disrupted_count).fillna(0).astype(int)
	cluster_df = cluster_df[['genus', 'conserved', 'fragmented', 'disrupted']]
	unique_genus_df = cluster_df.groupby('genus').first().reset_index()
	total_counts = unique_genus_df[['conserved', 'fragmented', 'disrupted']].sum(axis=1)

	unique_genus_df['c_ratio'] = round(unique_genus_df['conserved'] / total_counts, 2)
	unique_genus_df['f_ratio'] = round(unique_genus_df['fragmented'] / total_counts, 2)
	unique_genus_df['d_ratio'] = round(unique_genus_df['disrupted'] / total_counts, 2)

	unique_genus_df.fillna(0, inplace=True)

	unique_genus_df.to_csv(f"{output_dir}/tsv/cluster_type_genus.tsv", sep='\t', index=False)

	visualize_cluster_type(output_dir, unique_genus_df)

	end_time = time.time()
	total = end_time - start_time
	print(f"   └── classification done (elapsed time: {round(total / 60, 3)} min)")
	print("----------------------------------------------------------")