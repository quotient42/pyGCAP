import argparse
from .main import find_gene_cluster

def main():
    parser = argparse.ArgumentParser(description="Find Gene Cluster")
    parser.add_argument("working_dir", type=str, help="Working directory path")
    parser.add_argument("TAXON", type=str, help="Taxon identifier")
    parser.add_argument("probe_path", type=str, help="Probe file path")
    
    args = parser.parse_args()
    
    find_gene_cluster(args.working_dir, args.TAXON, args.probe_path)

if __name__ == "__main__":
    main()
