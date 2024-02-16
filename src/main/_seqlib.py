import os
import pandas as pd
from Bio.Blast import NCBIXML
import time
import sys

sys.path.append("../seqlib/")
from collect_data import collect_target_data
from search_repseq import search_repseq
from construct_lib import construct_lib

def process_raw_data(seqlib_dir):
	xml_path = os.path.join(seqlib_dir, "blastp_result.xml")
	tsv_path = os.path.join(seqlib_dir, "blastp_result.tsv")

	with open(xml_path, "r") as result_file:
			blast_records = NCBIXML.parse(result_file)
			data = {'Query': [], 'Accession': [], 'Length': [], 'E-value': [], 'Score': [], 'Identities': [], 'Coverage': [], 'Organism': []}

			for blast_record in blast_records:
					query_id = blast_record.query_id
					query_len = blast_record.query_length
					for alignment in blast_record.alignments:
							for hsp in alignment.hsps:
									percent_identity = round((hsp.identities / hsp.align_length) * 100, 2)
									coverage = round((hsp.align_length / query_len) * 100, 2)
									data['Query'].append(query_id)
									data['Accession'].append(alignment.title)
									data['Length'].append(alignment.length)
									data['E-value'].append(hsp.expect)
									data['Score'].append(hsp.score)
									data['Identities'].append(percent_identity)
									data['Coverage'].append(coverage)

									hit_def = alignment.hit_def
									organism_start = hit_def.find('[')
									organism_end = hit_def.find(']')
									organism_name = hit_def[organism_start + 1:organism_end]
									data['Organism'].append(organism_name)

			for i, name in enumerate(data['Accession']):
					split_name = name.split('|')
					ref_index = split_name.index('ref')
					if ref_index != -1 and ref_index + 1 < len(split_name):
							data['Accession'][i] = split_name[ref_index + 1]

			df = pd.DataFrame(data)

			df.to_csv(tsv_path, sep='\t', index=False)

def create_seqlib(project_info):
	print("<< creating seqlib...")
	data_dir = project_info['data']
	seqlib_dir = project_info['seqlib']

	# collect_target_data(project_info)
	process_raw_data(seqlib_dir)

	df = pd.read_csv(f"{seqlib_dir}/blastp_result.tsv", sep='\t')
	ref_df = pd.read_csv(f"{data_dir}/target.tsv", sep='\t')

	prediction_list = []
	for index, row in df.iterrows():
			for _, ref_row in ref_df.iterrows():
					if row['Query'] in ref_row['Accession']:
							prediction_list.append(ref_row['prediction'])
							break

	df['prediction'] = prediction_list
	col1 = df.columns[:-2].to_list()
	col2 = df.columns[-1:].to_list()
	new_col = col2 + col1
	df = df[new_col]

	output_tsv_file = "seqlib_data1.tsv"
	output_tsv_path = os.path.join(f"{seqlib_dir}", output_tsv_file)
	df.to_csv(output_tsv_path, sep='\t', index=False)

	combined_df = df[['prediction', 'Accession']]
	combined_df = combined_df.drop_duplicates()

	output_tsv_file = "seqlib_data2.tsv"
	output_tsv_path = os.path.join(f"{seqlib_dir}", output_tsv_file)
	combined_df.to_csv(output_tsv_path, sep='\t', index=False)

	search_repseq(project_info)
	construct_lib(project_info)