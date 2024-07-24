import argparse
import os
import requests
from xml.etree import ElementTree
from .main import find_gene_cluster

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

def check_argument(working_dir, TAXON, probe_path):
    if not os.path.exists(working_dir):
        print("<< Working directory does not exist. Please provide a valid file path.")
        exit(1)
    if not os.path.exists(probe_path):
        print("<< Probe file does not exist. Please provide a valid file path.")
        exit(1)
    if TAXON.isdigit():
        taxid = TAXON
        TAXON = get_taxon_name_from_taxid(TAXON)
    else:
        taxid = search_taxon_by_name(TAXON)
    
    return taxid

def main():
    parser = argparse.ArgumentParser(description="Find Gene Cluster")
    parser.add_argument("TAXON", type=str, help="Taxon identifier")
    parser.add_argument("probe_path", type=str, help="Probe file path")
    
    parser.add_argument("-w", "--working_dir", type=str, default=".", help="Working directory path")
    parser.add_argument("-t", "--thread", type=int, default=50, help="")
    parser.add_argument("-i", "--identity", type=int, default=0.5, help="")
    parser.add_argument("-m", "--max_target_seqs", type=int, default=500, help="")
    parser.add_argument("-s", "--skip", action='append', default=[], choices=['ncbi', 'mmseqs2', 'parsing', 'uniprot', 'blastdb', 'all'], help="Options to skip steps, e.g., ncbi, mmseqs2, parsing, uniprot, blastdb or all")
    
    args = parser.parse_args()
    
    taxid = check_argument(args.working_dir, args.TAXON, args.probe_path)
    
    if args.thread < 1:
        print("Thread number must be an integer greater than 0.")
        exit(1)
    if args.identity < 0.4:
        print("Identity must be a value between 0 and 1.")
        exit(1)
    if args.max_target_seqs < 1:
        print("Max target seqs must be an integer greater than 0.")
        exit(1)

    skip_options = {
        'ncbi': False,
        'mmseqs2': False,
        'parsing': False,
        'uniprot': False,
        'blastdb': False,
    }
    
    if 'all' in args.skip:
        for key in skip_options.keys():
            skip_options[key] = True
    else:
        for opt in args.skip:
            skip_options[opt] = True
    
    find_gene_cluster(
        args.TAXON,
        taxid,
        args.probe_path,
        args.working_dir, 
        args.thread,
        args.identity,
        args.max_target_seqs,
        skip_ncbi=skip_options['ncbi'],
        skip_mmseqs2=skip_options['mmseqs2'],
        skip_uniprot=skip_options['uniprot'],
        skip_parsing=skip_options['parsing'],
        skip_blastdb=skip_options['blastdb']
    )

if __name__ == "__main__":
    main()
