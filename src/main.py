import time

from _utils import sum_assembly_report
from _utils import organize_input_dir
from _parse import parse_genome
from _blast import run_blastp
from _count import count_blastp_result
from _seqlib import create_seqlib
from _search import search_target
from _count import count_mmseq_result
from _cluster import cluster_target
from _cluster import classify_cluster_type
from _accessory import collect_accessory
from _profile import final_profile

def find_gene_cluster():
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

    input_len = sum_assembly_report(project_info)
    project_info['input_len'] = input_len
    organize_input_dir(project_info)
    # parse_genome(project_info)

    run_blastp(project_info)
    count_blastp_result(project_info)

    create_seqlib(project_info)
    search_target(project_info)
    count_mmseq_result(project_info)

    cluster_target(project_info)
    classify_cluster_type(project_info)

    collect_accessory(project_info)
    final_profile(project_info)

    end_time = time.time()
    total = end_time - start_time
    print("------------------------------------------")
    print(f"<< ALL DONE (elapsed time: {round(total / 60, 3)} min)")
    print("------------------------------------------")

find_gene_cluster()