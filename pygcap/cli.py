import argparse
from .main import find_gene_cluster

def main():
    parser = argparse.ArgumentParser(description="Find Gene Cluster")
    parser.add_argument("working_dir", type=str, help="Working directory path")
    parser.add_argument("TAXON", type=str, help="Taxon identifier")
    parser.add_argument("probe_path", type=str, help="Probe file path")
    
    parser.add_argument("--skip", nargs='+', help="Options to skip steps, e.g., ncbi mmseqs2 parsing uniprot or all")
    
    args = parser.parse_args()
    
    skip_options = {
        'all': False,
        'ncbi': False,
        'mmseqs2': False,
        'parsing': False,
        'uniprot': False,
        'blast': False,
    }
    if args.skip:
        for opt in args.skip:
            if opt == 'all':
                skip_options['ncbi'] = True
                skip_options['mmseqs2'] = True
                skip_options['parsing'] = True
                skip_options['uniprot'] = True
                skip_options['blast'] = True
            elif opt in skip_options:
                skip_options[opt] = True
    
    find_gene_cluster(
        args.working_dir, 
        args.TAXON, 
        args.probe_path,
        skip_ncbi=skip_options['ncbi'],
        skip_mmseqs2=skip_options['mmseqs2'],
        skip_uniprot=skip_options['uniprot'],
        skip_parsing=skip_options['parsing'],
        skip_parsing=skip_options['blast']
    )

if __name__ == "__main__":
    main()
