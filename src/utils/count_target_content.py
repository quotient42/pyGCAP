import pandas as pd

def count_target_by_species(file_path, target_names):
    with open(file_path, 'r') as f:
        accession = f.readline().split('#')[-1].strip()
        name = f.readline().split('#')[-1].strip()

    df = pd.read_csv(file_path, sep='\t', comment='#')

    target_cnt = df['prediction'].value_counts().reset_index()
    target_cnt.columns = ['protein', 'count']

    target_cnt = pd.DataFrame({'protein': target_names}).merge(target_cnt, on='protein', how='left').fillna(0)
    target_cnt = target_cnt.astype({'count': int})

    result_df = pd.DataFrame(columns=['accession', 'name'] + target_names)
    result_df.loc[0, 'accession'] = accession
    result_df.loc[0, 'name'] = name
    result_df.loc[0, target_cnt['protein']] = target_cnt['count'].values

    return result_df

def count_target_by_genus(output_dir, df):
  df['genus'] = df['name'].str.split().str[0]

  protein_columns = df.columns[2:-1]
  result_df = df.drop('name', axis=1).groupby('genus')[protein_columns].mean().reset_index()

  result_df = result_df.fillna(0)
  result_df = result_df.round(3)
  result_df.to_csv(f"{output_dir}/tsv/target_genus.tsv", sep='\t', index=False)