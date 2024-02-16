import os
import pandas as pd
import sys
import time

sys.path.append('../cluster/')
from search_cluster import search_cluster_main
from adjust_cluster import adjust_cluster_main

sys.path.append('../visualize/')
from visualize_cluster import visualize_cluster

def cluster_target(project_info):
    print("<< clustering target genes...")
    input_dir = project_info['input']
    input_len = project_info['input_len']
    output_dir = project_info['output']

    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    progress = 0
    start_time = time.time()

    for genus in all_directories:
        if ("output" in genus) or ("seqlib" in genus):
            continue
        laps_start = time.time()
        path = os.path.join(input_dir, genus)
        if not os.path.exists(f"{path}/tmp"):
            os.makedirs(f"{path}/tmp")

        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GCF" in accession:
                continue
            progress += 1
            print(f"   ├── processing... ({progress}/{input_len}): {accession}") 
            
            cur_dir = f"{input_dir}/{genus}/{accession}"

            with open(f"{cur_dir}/genome_summary.tsv", 'r') as file:
                lines = file.readlines()
                species = lines[1][1:].strip()

            genome_summary = pd.read_csv(f"{cur_dir}/genome_summary.tsv", sep='\t', comment='#')
            target_list = pd.read_csv(f"{cur_dir}/target_seq.tsv", sep='\t', comment='#')
            input_info = [accession, cur_dir, species]

            cluster_summary = search_cluster_main(input_info, genome_summary, target_list)
            adjust_cluster_main(input_info, cluster_summary, f"{path}/tmp")

        print(f"   └── drawing cluster of {genus}...")
        visualize_cluster(f"{path}/tmp", f"{output_dir}/genus", genus)

    end_time = time.time()
    total = end_time - start_time
    print(f"   └── clustring done (elapsed time: {round(total / 60, 3)} min)")
    print("----------------------------------------------------------")