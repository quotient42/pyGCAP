import pandas as pd

def search_repseq(project_info):
	seqlib_dir = project_info['seqlib']
	data_dir = project_info['data']

	seq_df = pd.read_csv(f"{data_dir}/result_cluster.tsv", sep="\t").copy()

	seq_df['subseq'] = seq_df['subseq'].apply(lambda x: float(x.replace('WP_', '')))
	seq_df = seq_df.sort_values('subseq')

	target_df = pd.read_csv(f"{seqlib_dir}/seqlib_data2.tsv", sep='\t')
	target_df['Accession_tmp'] = target_df['Accession'].apply(lambda x: float(x.replace('WP_', '')))

	target_df['repseq'] = "-"

	for index, row in target_df.iterrows():
			accession = row['Accession_tmp']
			low = 0
			high = len(seq_df) - 1
			while low <= high:
					mid = (low + high) // 2
					if seq_df.iloc[mid]['subseq'] == accession:
							target_df.at[index, 'repseq'] = seq_df.iloc[mid]['repseq']
							break
					elif seq_df.iloc[mid]['subseq'] < accession:
							low = mid + 1
					else:
							high = mid - 1

	target_df = target_df[target_df['repseq'] != '-']
	target_df = target_df[['prediction', 'Accession', 'repseq']]
	target_df.to_csv(f"{seqlib_dir}/seqlib_rep.tsv", sep='\t', index=False)