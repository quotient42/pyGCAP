import pandas as pd

def make_target_fasta(project_info):
	data_dir = project_info['data']

	df = pd.read_csv(f"{data_dir}/target.tsv", sep="\t")
	output_file = f"{data_dir}/target.fasta"
	with open(output_file, 'w') as fasta_out:
		for index, row in df.iterrows():
			header = f">{row['name']} {row['protein id']} [{row['organism']}]\n"
			translation = row['translation']
			translation_lines = [translation[i:i+50] for i in range(0, len(translation), 50)]
			sequence = '\n'.join(translation_lines) + '\n'
			fasta_out.write(header + sequence)
