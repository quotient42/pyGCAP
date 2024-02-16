import pandas as pd
import sys

sys.path.append("../cluster/")
sys.path.append("../parse/")
sys.path.append("../utils/")
from parse_assembly_report import parse_assembly_report
from organize_input_dir import organize_input_dir
from _parse import parse_genome_data
from _seqlib import create_seqlib
from _search import search_target_seq
from _count import count_target
from _cluster import cluster_target
from classify_cluster import classify_target_cluster

import os

def test(project_info):
    input_dir = project_info['input']
    input_len = project_info['input_len']
    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    progress = 0

    for genus in all_directories:
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GCF" in accession:
                continue
            progress += 1
            print(f"   ├── processing... ({progress}/{input_len}): {accession}") 
            
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
                    
                    # Split translation into lines of 50 characters each
                    for i in range(0, len(translation), 50):
                        fasta_file.write(f"{translation[i:i+50]}\n")


with open("./project_info.tsv", 'r') as file:
    lines = file.readlines()
    if len(lines) < 2:
        print("<< Initialize your project first.")
        exit()
    project_info_data = lines[1].strip().split('\t')
    project_info = {
        'project_name': project_info_data[0],
        'root': project_info_data[1],
        'input': project_info_data[2],
        'output': project_info_data[3],
        'data': project_info_data[4],
        'seqlib': project_info_data[5]
		}

input_len = parse_assembly_report(project_info)
project_info['input_len'] = input_len
organize_input_dir(project_info)

project_info['input_len'] = 696
test(project_info)
# parse_genome_data(project_info)
# create_seqlib(project_info)
# search_target_seq(project_info)
# count_target(project_info)
# cluster_target(project_info)
# classify_target_cluster(project_info)
