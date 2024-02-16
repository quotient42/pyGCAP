import os
import pandas as pd
import sys
import time

sys.path.append('../search/')
sys.path.append('../utils/')
from search_by_seq import search_target_by_seq_main

def search_target_seq(project_info):
    print("<< searching target based on seq data...")
    input_dir = project_info['input']
    input_len = project_info['input_len']
    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    seq_lib = pd.read_csv(f"{project_info['seqlib']}/seqlib.tsv", sep='\t')
    seq_lib['prediction'] = seq_lib['prediction'].astype(str)
    seq_lib['RepSeq'] = seq_lib['RepSeq'].astype(str)
    seq_lib['Total'] = seq_lib['Total'].astype(int)
    seq_lib['Seq_List'] = seq_lib['Seq_List'].astype(str)

    progress = 0
    start_time = time.time()

    for genus in all_directories:
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GCF" in accession:
                continue
            progress += 1
            print(f"   ├── searching... ({progress}/{input_len}): {accession}") 
            
            cur_dir = f"{input_dir}/{genus}/{accession}"

            with open(f"{cur_dir}/genome_summary.tsv", 'r', encoding='utf-8') as file:
                lines = file.readlines()
                species = lines[1][1:].strip()

            genome_summary = pd.read_csv(f"{cur_dir}/genome_summary.tsv", sep='\t', comment='#')
            input_info = [accession, cur_dir, species]
            target_seq = search_target_by_seq_main(input_info, genome_summary, seq_lib)  # target_seq.tsv

    end_time = time.time()
    total = end_time - start_time
    print(f"   └── searching done (elapsed time: {round(total / 60, 3)} min)")
    print("----------------------------------------------------------")
