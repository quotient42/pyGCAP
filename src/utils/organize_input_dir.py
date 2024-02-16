import pandas as pd
import os
import shutil

def organize_input_dir(project_info):
    assembly_df = pd.read_csv(f"{project_info['output']}/tsv/assembly_report_sum.tsv", sep='\t')
    assembly_df = assembly_df[['Accession', 'Genus']]

    grouped_df = assembly_df.groupby('Genus')['Accession'].agg(list).reset_index()
    grouped_df['Total'] = grouped_df['Accession'].apply(len)

    input_dir = project_info['input']
    for index, row in grouped_df.iterrows():
        genus = row['Genus']
        genus_dir = os.path.join(input_dir, genus)

        if not os.path.exists(genus_dir):
            os.makedirs(genus_dir)

    for index, row in grouped_df.iterrows():
        genus = row['Genus']
        accessions = row['Accession']
        
        for accession in accessions:
            source_path = os.path.join(input_dir, accession)
            destination_path = os.path.join(input_dir, genus, accession)
            
            if os.path.exists(source_path):
                shutil.move(source_path, destination_path)
            
    for directory in os.listdir(input_dir):
        directory_path = os.path.join(input_dir, directory)
        if os.path.isdir(directory_path) and not os.listdir(directory_path):
            os.rmdir(directory_path)
        elif os.path.isdir(directory_path) and os.listdir(directory_path):
            tmp_dir = os.path.join(directory_path, 'tmp')
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)

