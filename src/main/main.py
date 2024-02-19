import time
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
from _accessory import collect_accessory
from classify_cluster import classify_target_cluster

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

start_time = time.time()

input_len = parse_assembly_report(project_info)
project_info['input_len'] = input_len
organize_input_dir(project_info)

project_info['input_len'] = 696
# parse_genome_data(project_info)
# create_seqlib(project_info)
# search_target_seq(project_info)
# count_target(project_info)
# cluster_target(project_info)
# classify_target_cluster(project_info)

collect_accessory(project_info)

end_time = time.time()
total = end_time - start_time
print(f"ALL DONE (elapsed time: {round(total / 60, 3)} min)")
