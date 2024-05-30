import pandas as pd
import requests
import shutil

def fetch_protein_data(accession):
    base_url = "https://www.ebi.ac.uk/proteins/api/proteins"
    url = f"{base_url}/{accession}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        data = response.json()
        
        protein_name = "-"
        gene_name = "-"
        protein_families = "-"
        ec = "-"
        pfam = []
        organism_name = "-"
        organism_id = "-"
        length = 0
        sequence = None
        refseq = "-"
        
        if 'protein' in data:
            protein_name_data = data['protein'].get('recommendedName', {}).get('fullName', "-")
            if protein_name_data == "-":
                protein_name_data = data['protein'].get('submittedName', ["-"])[0]
            protein_name = protein_name_data.get('value', "-")
        else:
            protein_name = "-"

        if 'gene' in data:
            gene_name_data = data['gene'][0].get('name', "-")
            if gene_name_data != "-":
                gene_name = gene_name_data.get('value', "-")
        
        if 'comments' in data:
            for comment in data['comments']:
                if comment.get('type') == 'FUNCTION':
                    protein_families = comment.get('text', "-")[0].get('value', "-")
                elif comment.get('type') == 'CATALYTIC_ACTIVITY':
                    ec = comment.get('reaction', {}).get('ecNumber', "-")
                    
        if 'dbReferences' in data:
            for xref in data['dbReferences']:
                if xref.get('type') == 'Pfam':
                    pfam.append(xref.get('id', "-"))
                elif xref.get('type') == 'RefSeq':
                    refseq_id = xref["id"]
                    if refseq_id.startswith("WP_"):
                        refseq = refseq_id
        
        if 'organism' in data:
            organism_name = data['organism'].get('names', "-")[0].get('value', "-")
            organism_id = str(data['organism'].get('taxonomy', "-"))
        
        if 'sequence' in data:
            length = data['sequence'].get('length', 0)
            sequence = data['sequence'].get('sequence', None)
        
        return {
            'protein_id': refseq,
            'accession': accession,
            'protein_name': protein_name,
            'gene_name': gene_name,
            'ec': ec,
            'organism_name': organism_name,
            'organism_id': organism_id,
            'pfam': pfam,
            'protein_families': protein_families,
            'length': length,
            'sequence': sequence
        }
    except requests.exceptions.RequestException as e:
        print(f"<< [Error] Skip for {accession}: {e}")
        return {}

def save_fasta(df, filename):
    print("<< Processing target data into FASTA file...")

    with open(filename, "w") as file:
        for index, row in df.iterrows():
            if row['protein_id'] == '-' or pd.isna(row['protein_id']):
                continue
            header = f">{row['Probe Name']} {row['protein_id']} [{row['organism_name']}]\n"
            sequence = str(row['sequence'])
            seq_lines = [sequence[i:i+50] for i in range(0, len(sequence), 50)]
            file.write(header)
            for line in seq_lines:
                file.write(line + "\n")

def process_target_data(project_info, target_original_path):
    target_new_path = f"{project_info['data']}/target.tsv"
    meta_file_path = f"{project_info['data']}/metadata_target.tsv"
    fasta_file_path = f"{project_info['data']}/target.fasta"
    
    try:
        shutil.copy(target_original_path, target_new_path)
        print("<< Searching targets from UniprotKB...")
    except FileNotFoundError:
        print("<< [Error] File not found. Please provide a valid file path.")
    
    df = pd.read_csv(target_new_path, sep='\t')

    data_list = []
    for accession in df['accession']:
        if accession == '-':
            data = {}
        else:
            data = fetch_protein_data(accession)
        data_list.append(data)

    df_result = pd.DataFrame(data_list)

    result_df = pd.concat([df, df_result], axis=1)

    result_df.to_csv(meta_file_path, sep='\t', index=False)
    save_fasta(result_df, fasta_file_path)