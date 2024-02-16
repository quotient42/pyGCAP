import pandas as pd
import subprocess

def blast_data(project_info):
	print("<< makeing blastdb...")
	db_path = f"{project_info['seqlib']}/lacto.aa"
	data_path = f"{project_info['seqlib']}/all.fasta"
	target_path = f"{project_info['data']}/target.fasta"
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
			'-query', target_path,
			'-db', db_path,
			'-evalue', "0.001",
			'-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
			'-out', output1_path,
	]

	subprocess.run(blastp_command)

	header_string = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"

	with open(output1_path, 'r') as file:
		content = file.read()

	content_with_header = header_string + content

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
