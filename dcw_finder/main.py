import time
import subprocess
import os
import glob

from ._utils import sum_assembly_report
from ._utils import organize_input_dir
from ._parse import parse_genome
from ._blast import run_blastp
from ._count import count_blastp_result
from ._seqlib import create_seqlib
from ._search import search_target
from ._count import count_mmseq_result
from ._cluster import cluster_target
from ._cluster import classify_cluster_type
from ._accessory import collect_accessory
from ._profile import final_profile
from ._target import process_target_data

def process_genome_data(working_dir, TAXON):
    ROOT_DIR = f"{working_dir}/{TAXON}"
    DATA_DIR = f"{ROOT_DIR}/data"
    INPUT_DIR = f"{ROOT_DIR}/input"

    answer = input(">> Need to download new genome data? (y/n): ")
    if answer.lower() == "y":
        print("<< Downloading genome data...")
        subprocess.run(["datasets", "download", "genome", "taxon", TAXON, "--reference", "--include", "protein,gff3,gbff", "--filename", f"{working_dir}/{TAXON}.zip"])
    else:
        print("<< Skip.")

    answer = input(">> Need to unzip genome data? (y/n): ")
    if answer.lower() == "y":
        print("<< Unzipping files...")
        subprocess.run(["unzip", f"{working_dir}/{TAXON}.zip", "-d", f"{working_dir}/{TAXON}"])
        
        os.rename(f"{ROOT_DIR}/ncbi_dataset/data", INPUT_DIR)
        os.rmdir(f"{ROOT_DIR}/ncbi_dataset")
        if not os.path.exists(DATA_DIR):
            os.mkdir(DATA_DIR)
    else:
        print("<< Skip.")

    answer = input(">> Need to format genome and protein data? (y/n): ")
    if answer.lower() == "y":
        print("<< Formatting data...")
        subprocess.run(f"dataformat tsv genome --inputfile {INPUT_DIR}/assembly_data_report.jsonl --fields accession,organism-name,organism-tax-id,source_database,checkm-completeness,checkm-contamination,assmstats-total-sequence-len,annotinfo-featcount-gene-total,annotinfo-featcount-gene-protein-coding,assminfo-bioproject,assminfo-biosample-accession,wgs-project-accession,wgs-url > {DATA_DIR}/assembly_report.tsv", shell=True)

        print("<< Concatenating protein files...")
        protein_files = glob.glob(f"{INPUT_DIR}/GCF*/protein.faa")
        with open(f"{DATA_DIR}/{TAXON}.faa", "wb") as outfile:
            for protein_file in protein_files:
                with open(protein_file, "rb") as infile:
                    outfile.write(infile.read())
    else:
        print("<< Skip.")

    answer = input(">> Need to run mmseqs2? (y/n): ")
    if answer.lower() == "y":
        print("<< Clustering protein sequences...")
        subprocess.run(["mmseqs", "easy-cluster", f"{DATA_DIR}/{TAXON}.faa", f"{DATA_DIR}/result", f"{DATA_DIR}/tmp", "--min-seq-id", "0.5", "--threads", "55"])
    else:
        print("<< Skip.")

    print("<< Organizing directories...")
    output_dir = f"{ROOT_DIR}/output"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    directories = ["genus", "img", "tsv"]
    for directory in directories:
        subdir = os.path.join(output_dir, directory)
        if not os.path.exists(subdir):
            os.mkdir(subdir)

    seqlib_dir = f"{ROOT_DIR}/seqlib"
    if not os.path.exists(seqlib_dir):
        os.mkdir(seqlib_dir)

    print("<< Data organization complete.")


def find_gene_cluster(working_dir, TAXON, target_path):
    if not os.path.exists(working_dir):
        print("<< Working directory does not exist. Please provide a valid file path.")
        return
    if not os.path.exists(target_path):
        print("<< Target file does not exist. Please provide a valid file path.")
        return

    start_time = time.time()

    process_genome_data(working_dir, TAXON)

    project_info = {
            'project_name': TAXON + " gene cluster",
            'root': f"{working_dir}/{TAXON}",
            'input': f"{working_dir}/{TAXON}/input",
            'output': f"{working_dir}/{TAXON}/output",
            'data': f"{working_dir}/{TAXON}/data",
            'seqlib': f"{working_dir}/{TAXON}/seqlib"
            }

    input_len = sum_assembly_report(project_info)
    project_info['input_len'] = input_len
    organize_input_dir(project_info)
    parse_genome(project_info)

    process_target_data(project_info, target_path)
    run_blastp(project_info)
    count_blastp_result(project_info)

    create_seqlib(project_info)
    search_target(project_info)
    count_mmseq_result(project_info, TAXON)

    cluster_target(project_info)
    classify_cluster_type(project_info)

    collect_accessory(project_info)
    final_profile(project_info)

    end_time = time.time()
    total = end_time - start_time
    print("------------------------------------------")
    print(f"<< ALL DONE (elapsed time: {round(total / 60, 3)} min)")
    print("------------------------------------------")
