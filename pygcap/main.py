import time
import subprocess
import os
import glob


from ._utils import sum_assembly_report
from ._utils import organize_input_dir
from ._parse import parse_genome
from ._blast import make_blastdb, run_blast
from ._count import count_blastp_result
from ._seqlib import create_seqlib
from ._search import search_probe
from ._count import count_mmseq_result
from ._cluster import cluster_probe
from ._cluster import classify_cluster_type
from ._accessory import collect_accessory
from ._profile import final_profile
from ._probe import process_probe_data

def process_genome_data(working_dir, TAXON, dir_info):
    ROOT_DIR = dir_info['root']
    DATA_DIR = dir_info['data']
    INPUT_DIR = dir_info['input']

    print("<< Downloading genome data...")
    subprocess.run(["datasets", 
                    "download", 
                    "genome", 
                    "taxon", 
                    TAXON, 
                    "--reference", 
                    "--include", "protein,gff3,gbff",
                    "--filename", f"{working_dir}/{TAXON}.zip"])
    
    print("<< Unzipping files...")
    subprocess.run(["unzip", 
                    f"{working_dir}/{TAXON}.zip", 
                    "-d", f"{working_dir}/{TAXON}"])
    os.rename(f"{ROOT_DIR}/ncbi_dataset/data", INPUT_DIR)
    os.rmdir(f"{ROOT_DIR}/ncbi_dataset")
    if not os.path.exists(DATA_DIR):
        os.mkdir(DATA_DIR)
    
    print("<< Formatting data...")
    subprocess.run(f"dataformat tsv genome --inputfile {INPUT_DIR}/assembly_data_report.jsonl --fields accession,organism-name,organism-tax-id,source_database,checkm-completeness,checkm-contamination,assmstats-total-sequence-len,annotinfo-featcount-gene-total,annotinfo-featcount-gene-protein-coding,assminfo-bioproject,assminfo-biosample-accession,wgs-project-accession,wgs-url > {DATA_DIR}/assembly_report.tsv", shell=True)

    print("<< Concatenating protein files...")
    protein_files = glob.glob(f"{INPUT_DIR}/GCF*/protein.faa")
    with open(f"{DATA_DIR}/{TAXON}.faa", "wb") as outfile:
        for protein_file in protein_files:
            with open(protein_file, "rb") as infile:
                outfile.write(infile.read())

def process_protein_data(TAXON, dir_info, thread, identity):
    DATA_DIR = dir_info['data']
    
    print("<< Clustering protein sequences...")
    subprocess.run(["mmseqs", 
                    "easy-cluster", 
                    f"{DATA_DIR}/{TAXON}.faa", 
                    f"{DATA_DIR}/result", 
                    f"{DATA_DIR}/tmp", 
                    "--min-seq-id", f"{identity}", 
                    "--threads", f"{thread}"])

def organize_directory(dir_info):
    ROOT_DIR = dir_info['root']

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

def find_gene_cluster(TAXON, taxid, probe_path,
                      working_dir=".", thread=55, identity=0.5, max_target_seqs=500,
                      skip_ncbi=False, skip_mmseqs2=False, 
                      skip_parsing=False, skip_uniprot=False, skip_blastdb=False):
    print(f"<< TAXON info: {TAXON} [taxid={taxid}]")
    start_time = time.time()

    dir_info = {
            'project_name': TAXON + " gene cluster",
            'root': f"{working_dir}/{TAXON}",
            'input': f"{working_dir}/{TAXON}/input",
            'output': f"{working_dir}/{TAXON}/output",
            'data': f"{working_dir}/{TAXON}/data",
            'seqlib': f"{working_dir}/{TAXON}/seqlib"
            }
    
    if not skip_ncbi:
        process_genome_data(working_dir, TAXON, dir_info)
    if not skip_mmseqs2:
        process_protein_data(TAXON, dir_info, thread, identity)
    organize_directory(dir_info)
    print("<< Data organization complete.")

    input_len = sum_assembly_report(dir_info)
    dir_info['input_len'] = input_len
    organize_input_dir(dir_info)

    if not skip_parsing:
        parse_genome(dir_info)
    if not skip_uniprot:
        process_probe_data(dir_info, probe_path)
    if not skip_blastdb:
        make_blastdb(dir_info, thread)

    run_blast(dir_info, thread, max_target_seqs)
    count_blastp_result(dir_info)

    create_seqlib(dir_info)
    search_probe(dir_info)
    count_mmseq_result(dir_info, TAXON)

    cluster_probe(dir_info)
    classify_cluster_type(dir_info)

    collect_accessory(dir_info)
    final_profile(dir_info)

    end_time = time.time()
    total = end_time - start_time
    print("------------------------------------------")
    print(f"<< ALL DONE (elapsed time: {round(total / 60, 3)} min)")
    print("------------------------------------------")
