''' 
    - sum_assembly_report  ˑˑˑˑˑˑ  to create summarized tsv file from assembly report.
    - organize_input_dir  ˑˑˑˑˑˑˑ  to organize input directories based on assembly report summary.

'''

import pandas as pd
import os
import shutil
import random

#=======================================================================================
def sum_assembly_report(project_info):
    raw_df = pd.read_csv(f"{project_info['data']}/assembly_report.tsv", sep='\t')

    completeness_filter = 80
    contamination_filter = 10

    low_completeness, high_contamination = _filter_and_save_completeness(project_info, raw_df, completeness_filter, contamination_filter)
    unique_genus, missclassified_species = _process_and_save_summary(project_info, raw_df)

    print("------------------------------------------")
    print("  *  total number of genomes:\t{}".format(raw_df.shape[0]))
    print("  *  completeness < {}%:\t{}".format(completeness_filter, low_completeness))
    print("  *  contamination > {}%:\t{}".format(contamination_filter, high_contamination))
    print("  *  number of different genus:\t{}".format(unique_genus))
    print("  *  missclassified species:\t{}".format(missclassified_species))
    print("------------------------------------------")

    return raw_df.shape[0]

def _filter_and_save_completeness(project_info, raw_df, completeness_filter, contamination_filter):
    low_completeness = raw_df[raw_df['CheckM completeness'] <= completeness_filter]
    low_completeness = low_completeness[['Assembly Accession', 'Organism Name', 'CheckM completeness']]
    low_completeness.to_csv(f"{project_info['output']}/tsv/completeness_{str(completeness_filter)}.tsv", sep='\t', index=False)

    high_contamination = raw_df[raw_df['CheckM contamination'] >= contamination_filter]
    high_contamination = high_contamination[['Assembly Accession', 'Organism Name', 'CheckM contamination']]
    high_contamination.to_csv(f"{project_info['output']}/tsv/contamination_{str(contamination_filter)}.tsv", sep='\t', index=False)

    return low_completeness.shape[0], high_contamination.shape[0]

def _process_and_save_summary(project_info, raw_df):
    df = raw_df.copy()
    df['Organism Name'] = df['Organism Name'].apply(lambda x: ' '.join(x.split()[:2]))
    df[['Genus', 'Species']] = df['Organism Name'].str.split(n=2, expand=True)

    missclassified_species = df[df['Genus'].str.startswith('[')].shape[0]
    df = df[~df['Genus'].str.startswith('[')]

    df = df.sort_values(by='Organism Name').reset_index(drop=True)
    df = df[['Assembly Accession', 'Genus', 'Species', 'Assembly Stats Total Sequence Length', 
             'CheckM completeness', 'CheckM contamination', 'Annotation Count Gene Total', 'Organism Name']]
    df.columns = ['Accession', 'Genus', 'Species', 'Seq Length', 'Completeness', 'Contamination', 'Annotation Count Gene', 'Name']

    df.to_csv(f"{project_info['output']}/tsv/assembly_report_sum.tsv", sep='\t', index=False)

    genus = df['Genus']
    unique_genus = genus.nunique()
    genus_cnt = genus.value_counts()

    df_genus = pd.DataFrame(genus_cnt.reset_index().values, columns=['genus name', 'total'])
    df_genus['%'] = (df_genus['total'] / df_genus['total'].sum()) * 100
    df_genus['%'] = df_genus['%'].apply(lambda x: round(x, 2))

    df_genus.to_csv(f"{project_info['output']}/tsv/genus_summary.tsv", sep='\t', index=False)

    return unique_genus, missclassified_species

#=======================================================================================
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

#=======================================================================================
def generate_random_color():
    return '#{:06x}'.format(random.randint(0, 0xFFFFFF))

def remove_non_numeric(input_string):
    start_index = 0
    while start_index < len(input_string) and not input_string[start_index].isdigit():
        start_index += 1

    end_index = len(input_string) - 1
    while end_index >= 0 and not input_string[end_index].isdigit():
        end_index -= 1

    return input_string[start_index:end_index + 1]