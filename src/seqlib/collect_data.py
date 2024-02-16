import pandas as pd
import os

def collect_target_data(project_info):
    print("<< collecting fasta data...")
    input_dir = project_info['input']
    seqlib_dir = project_info['seqlib']
    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    for genus in all_directories:
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GCF" in accession:
                continue
            
            cur_dir = f"{input_dir}/{genus}/{accession}"

            with open(f"{cur_dir}/genome_summary.tsv", 'r', encoding='utf-8') as file:
                lines = file.readlines()
                species = lines[1][1:].strip()

            genome_summary_df = pd.read_csv(f"{cur_dir}/genome_summary.tsv", sep='\t', comment='#')
            input_info = [accession, cur_dir, species]

            output_file = f"{cur_dir}/{accession}.fasta"
            with open(output_file, 'w') as fasta_file:
                for index, row in genome_summary_df.iterrows():
                    protein_id = row['protein_id']
                    gene = row['gene']
                    product = row['product']
                    translation = row['translation']

                    fasta_file.write(f">{protein_id} {gene} {product} [{species}]\n")
                    
                    for i in range(0, len(translation), 50):
                        fasta_file.write(f"{translation[i:i+50]}\n")

            with open(output_file, 'r') as fasta_file:
                fasta_content = fasta_file.read()
                with open(f'{seqlib_dir}/all.fasta', 'a') as all_fasta_file:
                    all_fasta_file.write(fasta_content)
            