import pandas as pd

pd.set_option('display.max_rows', None)

def filter_and_save_completeness(project_info, raw_df, completeness_filter, contamination_filter):
    low_completeness = raw_df[raw_df['CheckM completeness'] <= completeness_filter]
    low_completeness = low_completeness[['Assembly Accession', 'Organism Name', 'CheckM completeness']]
    low_completeness.to_csv(f"{project_info['output']}/tsv/completeness_{str(completeness_filter)}.tsv", sep='\t', index=False)

    high_contamination = raw_df[raw_df['CheckM contamination'] >= contamination_filter]
    high_contamination = high_contamination[['Assembly Accession', 'Organism Name', 'CheckM contamination']]
    high_contamination.to_csv(f"{project_info['output']}/tsv/contamination_{str(contamination_filter)}.tsv", sep='\t', index=False)

    return low_completeness.shape[0], high_contamination.shape[0]

def process_and_save_summary(project_info, raw_df):
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

def parse_assembly_report(project_info):
    raw_df = pd.read_csv(f"{project_info['data']}/assembly_report.tsv", sep='\t')

    completeness_filter = 80
    contamination_filter = 10

    low_completeness, high_contamination = filter_and_save_completeness(project_info, raw_df, completeness_filter, contamination_filter)
    unique_genus, missclassified_species = process_and_save_summary(project_info, raw_df)

    print("------------------------------------------")
    print("  *  total number of genomes:\t{}".format(raw_df.shape[0]))
    print("  *  completeness < {}%:\t{}".format(completeness_filter, low_completeness))
    print("  *  contamination > {}%:\t{}".format(contamination_filter, high_contamination))
    print("  *  number of different genus:\t{}".format(unique_genus))
    print("  *  missclassified species:\t{}".format(missclassified_species))
    print("------------------------------------------")

    return raw_df.shape[0]
