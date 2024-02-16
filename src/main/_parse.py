import os
import sys
import time

sys.path.append('../parse/')
sys.path.append('../utils/')
from parse_gbff import parse_gbff_main
from parse_gff import parse_gff_main
from merge_files import merge_gbff_gff

def parse_genome_data(base_info):
    print("<< parsing genome data...")
    base_dir_path = base_info['input']
    base_dir_len = base_info['input_len']
    all_directories = [d for d in os.listdir(base_dir_path) if os.path.isdir(os.path.join(base_dir_path, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    progress = 0
    start_time = time.time()

    for genus in all_directories:
        # print(f"* {genus}")
        laps_start = time.time()
        for accession in os.listdir(os.path.join(base_dir_path, genus)):
            if accession == "tmp":
                continue
            progress += 1
            print(f"   └── processing... ({progress}/{base_dir_len}): {accession}") 
            
            cur_dir = f"{base_dir_path}/{genus}/{accession}"
            gbff_input = f"{cur_dir}/genomic.gbff"
            gff_input = f"{cur_dir}/genomic.gff"
            
            input_info = [accession, cur_dir]

            species, gbff_summary = parse_gbff_main(gbff_input)
            input_info.append(species)
            gff_summary = parse_gff_main(gff_input)
            genome_summary = merge_gbff_gff(input_info, gbff_summary, gff_summary)

        # end_time = time.time()
        # total = end_time - start_time
        # print(f"└── done (elapsed time: {round(total / 60, 3)} min)\n")

    end_time = time.time()
    total = end_time - start_time
    print(f"   parsing done (elapsed time: {round(total / 60, 3)} min)")
    print(f"   check each GCF* directory to see genome summary")
    print("----------------------------------------------------------")