'''
	create_seqlib
	├── search_repseq  ˑˑˑˑˑˑˑˑ  to search representative sequence of each protein
	│  													 by using result cluster from mmseq2 and blast output
	└── construct_seqlib  ˑˑˑˑˑ  to create sequence library using representative sequences

'''

import time
import pandas as pd
from _utils import remove_non_numeric

#===============================================================================
def create_seqlib(project_info):
	print("<< creating seqlib...")
	start_time = time.time()

	search_repseq(project_info)
	construct_seqlib(project_info)

	end_time = time.time()
	total = end_time - start_time
	print(f"   └── seqlib created (elapsed time: {round(total / 60, 3)} min)")

#===============================================================================
def search_repseq(project_info):
	seqlib_dir = project_info['seqlib']
	data_dir = project_info['data']

	file_path = f"{data_dir}/result_cluster.tsv"

	with open(file_path, 'r') as file:
					lines = file.readlines()

	if not lines[0].startswith("repseq\tsubseq"):
			lines.insert(0, "repseq\tsubseq\n")

	with open(file_path, 'w') as file:
			file.writelines(lines)

	seq_df = pd.read_csv(file_path, sep="\t").copy()

	seq_df['subseq_tmp'] = seq_df['subseq'].apply(lambda x: float(remove_non_numeric(x)))
	seq_df = seq_df.sort_values('subseq_tmp')

	target_df = pd.read_csv(f"{seqlib_dir}/blast_output2.tsv", sep='\t')
	target_df['Accession_tmp'] = target_df['Accession'].apply(lambda x: float(remove_non_numeric(x)))

	target_df['repseq'] = "-"

	for index, row in target_df.iterrows():
			accession = row['Accession_tmp']
			low = 0
			high = len(seq_df) - 1
			while low <= high:
					mid = (low + high) // 2
					if seq_df.iloc[mid]['subseq_tmp'] == accession:
							target_df.at[index, 'repseq'] = seq_df.iloc[mid]['repseq']
							break
					elif seq_df.iloc[mid]['subseq_tmp'] < accession:
							low = mid + 1
					else:
							high = mid - 1

	target_df = target_df[target_df['repseq'] != '-']
	target_df = target_df[['Prediction', 'Accession', 'repseq', 'Query']]
	target_df.to_csv(f"{seqlib_dir}/seqlib_rep.tsv", sep='\t', index=False)

#===============================================================================
def construct_seqlib(project_info):
		seqlib_dir = project_info['seqlib']
		data_dir = project_info['data']

		seq_df = pd.read_csv(f"{data_dir}/result_cluster.tsv", sep='\t').copy()
		seq_df = seq_df.groupby('repseq').agg({'subseq': list}).reset_index()
		seq_df.rename(columns={'subseq': 'SeqList'}, inplace=True)
		seq_df['Total'] = seq_df['SeqList'].apply(len)
		seq_df.drop_duplicates(subset=['repseq'], keep='first', inplace=True)
		seq_df.rename(columns={'repseq': 'RepSeq'}, inplace=True)
		seq_df = seq_df[['RepSeq', 'Total', 'SeqList']]

		seq_df.to_csv(f"{seqlib_dir}/mmseq_rep.tsv", sep='\t', index=False)

		seq_df['RepSeq_tmp'] = seq_df['RepSeq'].apply(lambda x: float(remove_non_numeric(x)))
		seq_df = seq_df.sort_values('RepSeq_tmp')
		seq_df['Prediction'] = pd.Series([], dtype=object)

		target_text = pd.read_csv(f"{seqlib_dir}/seqlib_rep.tsv", sep='\t', comment='#')
		target_text['repseq_tmp'] = target_text['repseq'].apply(lambda x: float(remove_non_numeric(x)))
		
		for index, row in target_text.iterrows():
				low = 0
				high = len(seq_df) - 1
				while low <= high:
						mid = (low + high) // 2
						if seq_df.iloc[mid]['RepSeq_tmp'] == row['repseq_tmp']:
								seq_df.at[mid, 'Prediction'] = row['Prediction']
								break
						elif seq_df.iloc[mid]['RepSeq_tmp'] < row['repseq_tmp']:
								low = mid + 1
						else:
								high = mid - 1

		seq_df = seq_df.dropna(subset=['Prediction'])
		seq_df['RepName'] = ""
		seq_df = seq_df[['RepName', 'Prediction', 'RepSeq', 'Total', 'SeqList']]
		seq_df = seq_df.sort_values('Prediction')

		name_count = {}
		for index, row in seq_df.iterrows():
				prediction = row['Prediction']
				if prediction not in name_count:
						name_count[prediction] = 1
				else:
						name_count[prediction] += 1
				name = f"{prediction}.{name_count[prediction]}"
				seq_df.at[index, 'RepName'] = name

		seq_df.to_csv(f"{seqlib_dir}/seqlib.tsv", sep='\t', index=False)