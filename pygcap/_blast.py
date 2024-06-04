'''
	run_blastp
	├── create_probe_fasta  ˑˑˑˑ  to create FASTA file from probe data.
	├── create_input_fasta  ˑˑˑˑˑ  to create FASTA file from input genome data
	└── blast_data  ˑˑˑˑˑˑˑˑˑˑˑˑˑ  to get blast result from 'blastmakedb' and 'blastp'

'''

import time
import pandas as pd
import os
import subprocess

#===============================================================================
def run_blastp(project_info):
	input_dir = project_info['input']
	start_time = time.time()

	# create_probe_fasta(project_info)
	create_input_fasta(project_info)
	blast_data(project_info)      
	split_blast_result(project_info)  

	end_time = time.time()
	total = end_time - start_time
	print(f"   └── blastp complete (elapsed time: {round(total / 60, 3)} min)")

#===============================================================================
def create_probe_fasta(project_info):
	print("<< creating probe FASTA...")
	data_dir = project_info['data']

	df = pd.read_csv(f"{data_dir}/probe.tsv", sep="\t")
	output_file = f"{data_dir}/probe.fasta"
	with open(output_file, 'w') as fasta_out:
		for index, row in df.iterrows():
			header = f">{row['name']} {row['protein id']} [{row['organism']}]\n"
			translation = row['translation']
			translation_lines = [translation[i:i+50] for i in range(0, len(translation), 50)]
			sequence = '\n'.join(translation_lines) + '\n'
			fasta_out.write(header + sequence)

#===============================================================================
def create_input_fasta(project_info):
    print("<< creating input FASTA...")
    input_dir = project_info['input']
    seqlib_dir = project_info['seqlib']
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
                binomial1 = species.split(" ")[0]
                binomial2 = species.split(" ")[1]

            genome_summary_df = pd.read_csv(f"{cur_dir}/genome_summary.tsv", sep='\t', comment='#')
            input_info = [accession, cur_dir, species]

            output_file = f"{cur_dir}/{accession}.fasta"
            with open(output_file, 'w') as fasta_file:
                for index, row in genome_summary_df.iterrows():
                    protein_id = row['protein_id']
                    gene = row['gene']
                    product = row['product']
                    translation = row['translation']

                    fasta_file.write(f">{protein_id}|{accession}|{binomial1}_{binomial2} {gene} {product} [{species}]\n")
                    
                    for i in range(0, len(translation), 50):
                        fasta_file.write(f"{translation[i:i+50]}\n")

            with open(output_file, 'r') as fasta_file:
                fasta_content = fasta_file.read()
                with open(f'{seqlib_dir}/all.fasta', 'a') as all_fasta_file:
                    all_fasta_file.write(fasta_content)

#===============================================================================
def blast_data(project_info):
	print("<< makeing blastdb...")
	db_path = f"{project_info['seqlib']}/lacto.aa"
	data_path = f"{project_info['seqlib']}/all.fasta"
	probe_path = f"{project_info['data']}/probe.fasta"
	output1_path = f"{project_info['seqlib']}/blast_output1.tsv"
	output2_path = f"{project_info['seqlib']}/blast_output2.tsv"

	makeblastdb_command = [
    'makeblastdb',
    '-in', data_path,
    '-dbtype', 'prot',
    '-hash_index',
    '-out', db_path,
    '-title', '"LAB ncbi protein DB"'
	]

	subprocess.run(makeblastdb_command)

	print("<< running blastp...")
	blastp_command = [
			'blastp',
			'-query', probe_path,
			'-db', db_path,
			'-evalue', "0.001",
			'-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
			'-out', output1_path,
	]

	subprocess.run(blastp_command)

	header_string = "qseqid\tsseqid\torgid\torgname\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"

	with open(output1_path, 'r') as file:
		content = file.read()

	content_with_header = header_string + content.replace('|', '\t')

	with open(output1_path, 'w') as file:
		file.write(content_with_header)

	df = pd.read_csv(output1_path, sep='\t')

	df['coverage'] = round((df['qend'] - df['qstart'] + 1) / df['length'] * 100, 3)
	df.to_csv(output1_path, sep='\t', index=False)

	df = df[['qseqid', 'sseqid']]
	df.columns = ["Prediction", "Accession"]
	df['Query'] = df["Prediction"]
	df['Prediction'] = df['Prediction'].apply(lambda x: x.split('.')[0])
	df.to_csv(output2_path, sep='\t', index=False)

def split_blast_result(project_info):
	output1_path = f"{project_info['seqlib']}/blast_output1.tsv"

	input_dir = project_info['input']
	
	all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
	all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

	df = pd.read_csv(output1_path, sep='\t')

	for genus in all_directories:
		for accession in os.listdir(os.path.join(input_dir, genus)):
				if not "GC" in accession:
						continue

				cur_dir = f"{input_dir}/{genus}/{accession}"
				with open(f"{cur_dir}/genome_summary.tsv", 'r', encoding='utf-8') as file:
						lines = file.readlines()
						species = lines[1][1:].strip()
				
				matching_rows = df[df['orgid'] == accession]
				matching_rows = matching_rows.drop(columns=['orgid', 'orgname'])
				matching_rows = matching_rows.sort_values(by='qseqid')

				# with open(f"{cur_dir}/probe_blast.tsv", 'w', encoding='utf-8') as blast_file:
				# 	blast_file.write(f"# {accession}\n")
				# 	blast_file.write(f"# {species}\n")
				# 	matching_rows.to_csv(blast_file, sep='\t', index=False)

				genome_summary = pd.read_csv(f"{cur_dir}/genome_summary.tsv", sep='\t', comment='#')
				merged_df = pd.merge(matching_rows, genome_summary, left_on='sseqid', right_on='protein_id', how='inner')
				merged_df = merged_df.rename(columns={'qseqid': 'TarName'})
				merged_df['Prediction'] = merged_df['TarName'].apply(lambda x: x.split('.')[0])
				merged_df = merged_df[['TarName', 'Prediction', 'contig', 'protein_id', 'strand', 'start', 'end']]

				with open(f"{cur_dir}/probe_blast.tsv", 'w', encoding='utf-8') as blast_file:
					blast_file.write(f"# {accession}\n# {species}\n")
					merged_df.to_csv(blast_file, sep='\t', index=False)
				
				merged_df.drop_duplicates(subset='protein_id', keep='first', inplace=True)
				with open(f"{cur_dir}/probe_blast2.tsv", 'w', encoding='utf-8') as blast_file:
					blast_file.write(f"# {accession}\n# {species}\n")
					merged_df.to_csv(blast_file, sep='\t', index=False)
				