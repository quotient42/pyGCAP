import os
import pandas as pd
import time

def construct_lib(project_info):
		seqlib_dir = project_info['seqlib']
		data_dir = project_info['data']

		seq_df = pd.read_csv(f"{data_dir}/result_cluster.tsv", sep='\t').copy()
		seq_df = seq_df.groupby('repseq').agg({'subseq': list}).reset_index()
		seq_df.rename(columns={'subseq': 'Seq_List'}, inplace=True)
		seq_df['Total'] = seq_df['Seq_List'].apply(len)
		seq_df.drop_duplicates(subset=['repseq'], keep='first', inplace=True)
		seq_df.rename(columns={'repseq': 'RepSeq'}, inplace=True)
		seq_df = seq_df[['RepSeq', 'Total', 'Seq_List']]

		seq_df['RepSeq_tmp'] = seq_df['RepSeq'].apply(lambda x: float(x.replace('WP_', '')))
		seq_df = seq_df.sort_values('RepSeq_tmp')
		seq_df['Prediction'] = pd.Series([], dtype=object)

		start_time = time.time()

		target_text = pd.read_csv(f"{seqlib_dir}/seqlib_rep.tsv", sep='\t', comment='#')
		target_text['repseq_tmp'] = target_text['repseq'].apply(lambda x: float(x.replace('WP_', '')))
		
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
		seq_df = seq_df[['Prediction', 'RepSeq', 'Total', 'Seq_List']]
		seq_df = seq_df.sort_values('Prediction')

		seq_df.to_csv(f"{seqlib_dir}/seqlib.tsv", sep='\t', index=False)

		end_time = time.time()
		total = end_time - start_time
		print(f"   seqlib created (elapsed time: {round(total / 60, 3)} min)")
