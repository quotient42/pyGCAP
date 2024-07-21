import time
import subprocess
import os
import glob
import requests
from xml.etree import ElementTree

from ._utils import sum_assembly_report
from ._utils import organize_input_dir
from ._parse import parse_genome
from ._blast import run_blast
from ._count import count_blastp_result
from ._seqlib import create_seqlib
from ._search import search_probe
from ._count import count_mmseq_result
from ._cluster import cluster_probe
from ._cluster import classify_cluster_type
from ._accessory import collect_accessory
from ._profile import final_profile
from ._probe import process_probe_data

def get_taxon_name_from_taxid(taxid):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    params = {
        "db": "taxonomy",
        "id": taxid,
        "retmode": "xml"
    }
    response = requests.get(base_url + "efetch.fcgi", params=params)
    
    if response.status_code == 200:
        root = ElementTree.fromstring(response.content)
        scientific_name = root.find(".//ScientificName")
        if scientific_name is not None:
            return scientific_name.text
        else:
            print("[Error] Taxon not found. Check the taxid and try again.")
            exit()
    else:
        print ("[Error] something went wrong while searching taxon.")
        exit()

def search_taxon_by_name(name):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    params = {
        "db": "taxonomy",
        "term": name,
        "retmode": "xml"
    }
    response = requests.get(base_url + "esearch.fcgi", params=params)

    if response.status_code == 200:
        root = ElementTree.fromstring(response.content)
        id_list = root.findall(".//IdList/Id")
        if not id_list:
            print("[Error] No taxon found for the given name. Check the name and try again.")
            exit()
        else:
            taxid = id_list[0].text
            return taxid
    else:
        print("[Error] something went wrong while searching taxon by name.")
        exit()

def process_genome_data(working_dir, TAXON, skip_ncbi, skip_mmseqs2):
    if TAXON.isdigit():
        taxid = TAXON
        TAXON = get_taxon_name_from_taxid(TAXON)
    else:
        taxid = search_taxon_by_name(TAXON)

    print(f"<< TAXON info: {TAXON} [taxid={taxid}]")
    ROOT_DIR = f"{working_dir}/{TAXON}"
    DATA_DIR = f"{ROOT_DIR}/data"
    INPUT_DIR = f"{ROOT_DIR}/input"

    if not skip_ncbi:
        print("<< Downloading genome data...")
        subprocess.run(["datasets", "download", "genome", "taxon", TAXON, "--reference", "--include", "protein,gff3,gbff", "--filename", f"{working_dir}/{TAXON}.zip"])
        print("<< Unzipping files...")
        subprocess.run(["unzip", f"{working_dir}/{TAXON}.zip", "-d", f"{working_dir}/{TAXON}"])
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

    if not skip_mmseqs2:
        print("<< Clustering protein sequences...")
        subprocess.run(["mmseqs", "easy-cluster", f"{DATA_DIR}/{TAXON}.faa", f"{DATA_DIR}/result", f"{DATA_DIR}/tmp", "--min-seq-id", "0.5", "--threads", "55"])

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
    return TAXON


def find_gene_cluster(working_dir, TAXON, probe_path, 
                      skip_ncbi=False, skip_mmseqs2=False, 
                      skip_parsing=False, skip_uniprot=False, skip_blastdb=False):
    if not os.path.exists(working_dir):
        print("<< Working directory does not exist. Please provide a valid file path.")
        return
    if not os.path.exists(probe_path):
        print("<< Probe file does not exist. Please provide a valid file path.")
        return

    start_time = time.time()

    TAXON = process_genome_data(working_dir, TAXON, skip_ncbi, skip_mmseqs2)

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

    if not skip_parsing:
        parse_genome(project_info)
    if not skip_uniprot:
        process_probe_data(project_info, probe_path)
    run_blast(project_info, skip_blastdb)
    count_blastp_result(project_info)

    create_seqlib(project_info)
    search_probe(project_info)
    count_mmseq_result(project_info, TAXON)

    cluster_probe(project_info)
    classify_cluster_type(project_info)

    collect_accessory(project_info)
    final_profile(project_info)

    end_time = time.time()
    total = end_time - start_time
    print("------------------------------------------")
    print(f"<< ALL DONE (elapsed time: {round(total / 60, 3)} min)")
    print("------------------------------------------")
