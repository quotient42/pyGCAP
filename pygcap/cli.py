import argparse
from .main import find_gene_cluster

def main():
    parser = argparse.ArgumentParser(description="Find Gene Cluster")
    parser.add_argument("working_dir", type=str, help="Working directory path")
    parser.add_argument("TAXON", type=str, help="Taxon identifier")
    parser.add_argument("probe_path", type=str, help="Probe file path")
    
    parser.add_argument("-s", "--skip", action='append', default=[], choices=['ncbi', 'mmseqs2', 'parsing', 'uniprot', 'blastdb', 'all'], help="Options to skip steps, e.g., ncbi, mmseqs2, parsing, uniprot, blastdb or all")
    
    args = parser.parse_args()
    
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
        args.working_dir, 
        args.TAXON, 
        args.probe_path,
        skip_ncbi=skip_options['ncbi'],
        skip_mmseqs2=skip_options['mmseqs2'],
        skip_uniprot=skip_options['uniprot'],
        skip_parsing=skip_options['parsing'],
        skip_blastdb=skip_options['blastdb']
    )

if __name__ == "__main__":
    main()
