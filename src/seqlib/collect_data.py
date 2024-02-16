from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
import time

def blast_query(base_dir_path, accessions):
		for query_protein_id in accessions:
				print(f"> processing {query_protein_id}...")
				xml_filename = f"{base_dir_path}/seqlib/tmp/{query_protein_id}.xml"
				result_handle = NCBIWWW.qblast("blastp", "refseq_protein", query_protein_id, entrez_query="txid186826[ORGN]")

				with open(xml_filename, "w") as result_file:
						result_file.write(result_handle.read())
				result_handle.close()

def parse_blast_results(base_dir_path, accessions):
		for query_protein_id in accessions:
				xml_filename = f"{base_dir_path}/seqlib/tmp/{query_protein_id}.xml"
				with open(xml_filename, "r") as result_file:
						blast_records = NCBIXML.parse(result_file)
						data = {'Name': [], 'Length': [], 'E-value': [], 'Score': [], 'Identities': [], 'Coverage': []}

						for blast_record in blast_records:
								query_len = blast_record.query_length
								for alignment in blast_record.alignments:
										for hsp in alignment.hsps:
												percent_identity = round((hsp.identities / hsp.align_length) * 100, 2)
												coverage = round((hsp.align_length / query_len) * 100, 2)
												data['Name'].append(alignment.title)
												data['Length'].append(alignment.length)
												data['E-value'].append(hsp.expect)
												data['Score'].append(hsp.score)
												data['Identities'].append(percent_identity)
												data['Coverage'].append(coverage)
						df = pd.DataFrame(data)

						tsv_filename = f"{base_dir_path}/seqlib/tmp/{query_protein_id}.tsv"
						df.to_csv(tsv_filename, sep='\t', index=False)

def filter_and_combine_dfs(accessions, prediction, all_filtered_dfs, base_dir_path):
    for query_protein_id in accessions:
        tsv_filename = f"{base_dir_path}/seqlib/tmp/{query_protein_id}.tsv"
        filtered_df = pd.read_csv(tsv_filename, sep='\t')
        filtered_df = filtered_df[filtered_df['E-value'] <= 0.001]
        filtered_df = filtered_df[['Name']].copy()
        filtered_df['Accession'] = "-"
        filtered_df['prediction'] = prediction

        for index, name in filtered_df.iterrows():
            parts = name['Name'].split('|')
            if 'ref' in parts:
                ref_index = parts.index('ref')
                if ref_index + 1 < len(parts):
                    filtered_df.at[index, 'Accession'] = parts[ref_index + 1]
            else:
                filtered_df.at[index, 'Accession'] = "-"

        all_filtered_dfs.append(filtered_df)

    combined_df = pd.concat(all_filtered_dfs, ignore_index=True)
    combined_df = combined_df[combined_df['Accession'] != "-"]
    combined_df = combined_df[['prediction', 'Accession']]
    combined_df = combined_df.drop_duplicates()

    return combined_df

def collect_target_data(base_dir_path):
    df = pd.read_csv(f"{base_dir_path}/target.tsv", sep='\t')
    df['Accession'] = df['Accession'].apply(lambda x: [item.strip(" []'") for item in x.split(',')])
    all_filtered_dfs = []

    for index, row in df.iterrows():
        prediction = row['prediction']
        accessions = row['Accession']

        start_time = time.time()
        blast_query(base_dir_path, accessions)
        parse_blast_results(base_dir_path, accessions)
        end_time = time.time()
        total = end_time - start_time
        print(f"<< search for {prediction} done (elapsed time: {round(total / 60, 3)} min)")

        combined_df = filter_and_combine_dfs(accessions, prediction, all_filtered_dfs, base_dir_path)
        all_filtered_dfs.append(combined_df)

    combined_result = pd.concat(all_filtered_dfs, ignore_index=True)
    combined_result = combined_result[['prediction', 'Accession']]
    combined_result = combined_result.drop_duplicates()

    combined_result.to_csv(f"{base_dir_path}/seqlib/target_expanded.tsv", sep='\t', index=False)
		