import os
import pandas as pd
import sys
import time

sys.path.append('../utils/')
from count_target_content import count_target_by_species
from count_target_content import count_target_by_genus

sys.path.append('../visualize/')
from visualize_genus import visualize_by_genus
from visualize_species import visualize_by_species
from visualize_freq import visualize_freq

sys.path.append('../info/')
from target_info import target_list

def count_target(project_info):
    print("<< counting target by species and genus level...")
    input_dir = project_info['input']
    output_dir = project_info['output']

    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    start_time = time.time()

    target_names = []
    for target in target_list:
        target_names.append(target['protein'])
    target_cnt_list = []

    for genus in all_directories:
        genus_df_list = []
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GCF" in accession:
                continue
            file_path = f"{input_dir}/{genus}/{accession}/target_seq.tsv"
            tmp = count_target_by_species(file_path, target_names)
            genus_df_list.append(tmp)
        
        if genus_df_list:
            genus_df = pd.concat(genus_df_list, ignore_index=True)
            target_cnt_list.append(genus_df)

    target_cnt = pd.concat(target_cnt_list, ignore_index=True)

    target_cnt = target_cnt.sort_values("name")
    target_cnt.to_csv(f"{output_dir}/tsv/contents_species.tsv", sep='\t', index=False)

    count_target_by_genus(output_dir, target_cnt)
    visualize_by_species(output_dir)
    visualize_by_genus(output_dir)

    target_sum = target_cnt[target_names]
    target_sum = target_sum.sum(axis=0).reset_index()
    target_sum.to_csv(f"{output_dir}/tsv/target_frequency.tsv", sep='\t', index=False)

    visualize_freq(output_dir)

    end_time = time.time()
    total = end_time - start_time
    print(f"   └── counting done (elapsed time: {round(total / 60, 3)} min)")
    print("----------------------------------------------------------")

