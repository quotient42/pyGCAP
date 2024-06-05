import pandas as pd
import os
from ._visualize import visualize_heatmap

def final_profile(project_info):
	input_dir = project_info['input']
	seqlib_dir = project_info['seqlib']
	output_dir = project_info['output']

	seqlib = pd.read_csv(f"{seqlib_dir}/seqlib.tsv", sep='\t', comment='#')
	acclib = pd.read_csv(f"{output_dir}/tsv/acclib.tsv", sep='\t', comment='#')

	seqlib_repnames = seqlib['RepName'].tolist()
	acclib_repnames = acclib['RepName'].tolist()

	all_repnames = ['accession', 'name'] + seqlib_repnames + acclib_repnames

	profile = pd.DataFrame(columns=all_repnames)

	all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
	all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

	for genus in all_directories:
			for accession in os.listdir(os.path.join(input_dir, genus)):
					if not "GC" in accession:
							continue
					
					cur_dir = f"{input_dir}/{genus}/{accession}"

					with open(f"{cur_dir}/genome_summary.tsv", 'r', encoding='utf-8') as file:
							lines = file.readlines()
							species = lines[1][1:].strip()

					probe_final = pd.read_csv(f"{cur_dir}/probe_final.tsv", sep='\t', comment='#')

					new_row = pd.Series(0, index=profile.columns)
            
					new_row['accession'] = new_row['accession'].astype(str)
					new_row['accession'] = accession
					new_row['name'] = species
					
					for repname in probe_final['RepName']:
							if repname in profile.columns:
									new_row[repname] += 1
					
					profile = profile._append(new_row, ignore_index=True)

	profile = profile.sort_values('name')
	profile.to_csv(f"{output_dir}/tsv/final_profile_1.tsv", sep='\t', index=False)
	visualize_heatmap(output_dir, 'final_profile_1', 2)

	acc_columns = [col for col in profile.columns if col.startswith('Acc.')]
	simplified_profile = profile[['accession', 'name'] + acc_columns]
	simplified_profile.to_csv(f"{output_dir}/tsv/contents_acc_species.tsv", sep='\t', index=False)
	profile = profile.drop(columns=['accession', 'name'] + acc_columns)

	for col_name, col_data in profile.items():
			split_name = col_name.split('.')[0]
			if split_name not in simplified_profile.columns:
					simplified_profile[split_name] = 0
			simplified_profile[split_name] += col_data

	acc_columns = [col for col in simplified_profile.columns if col.startswith('Acc.')]
	other_columns = [col for col in simplified_profile.columns if not col.startswith('Acc.')]
	rearranged_columns = other_columns + acc_columns
	simplified_profile = simplified_profile[rearranged_columns]

	simplified_profile.to_csv(f"{output_dir}/tsv/final_profile_2.tsv", sep='\t', index=False)
	visualize_heatmap(output_dir, 'final_profile_2', 2)