''' 
    count_blastp_result  ˑˑˑˑˑˑˑ  count & visualize gene contents from blastp result
    count_mmseq_result  ˑˑˑˑˑˑˑˑ  count & visualize gene contents mmseq2 result
    ├── count_mmseq_species  
    └── count_mmseq_genus  
    count_mmseq_cluster  ˑˑˑˑˑˑˑ  count & visualize protein clusters from mmseq2 result
    
'''

import os
import pandas as pd
import sys
import time

from ._visualize import visualize_heatmap
from ._visualize import visualize_freq

#===============================================================================
def count_blastp_result(project_info):
    output_dir = project_info['output']

    result_df = pd.read_csv(f"{project_info['seqlib']}/blast_output1.tsv", sep='\t')
    target_df = pd.read_csv(f"{project_info['data']}/metadata_target.tsv", sep='\t')
    assembly_df = pd.read_csv(f"{project_info['output']}/tsv/assembly_report_sum.tsv", sep='\t')

    target = target_df['Probe Name'].values.tolist()
    target2 = target_df['Prediction'].values.tolist()

    name = assembly_df['Organism Name'].values.tolist() 
    accession = assembly_df['Accession'].values.tolist() 
    
    new_df = pd.DataFrame({
        'accession': accession,
        'name': name,
    })

    for t in target:
        new_df[t] = 0

    for index, row in result_df.iterrows():
        orgid = row['orgid']
        qseqid = row['qseqid']
        if orgid in new_df['accession'].values:
            idx = new_df[new_df['accession'] == orgid].index[0]
            new_df.loc[idx, qseqid] += 1

    new_df.to_csv(f"{output_dir}/tsv/contents_blastp_species1.tsv", sep='\t', index=False)

    target_column_indices = {key: [] for key in target2}
    first_row  = new_df.columns.tolist()

    for key in target_column_indices:
        for index, column_name in enumerate(first_row):
            if key in column_name:
                target_column_indices[key].append(index)

    collapsed_df = new_df[['accession', 'name']].copy()

    for key, indices in target_column_indices.items():
        new_col_values = []
        for index, row in new_df.iterrows():
            value_sum = sum(row[idx] for idx in indices)
            new_col_values.append(value_sum)
        collapsed_df[key] = new_col_values

    collapsed_df = collapsed_df.groupby(['accession', 'name'], as_index=False).sum()
    collapsed_df = collapsed_df.sort_values('name')
    collapsed_df.to_csv(f"{output_dir}/tsv/contents_blastp_species2.tsv", sep='\t', index=False)

    collapsed_df['genus'] = collapsed_df['name'].str.split().str[0]

    protein_columns = collapsed_df.columns[2:-1]
    result_df = collapsed_df.drop('name', axis=1).groupby('genus')[protein_columns].mean().reset_index()

    result_df = result_df.fillna(0)
    result_df = result_df.round(1)
    result_df.to_csv(f"{output_dir}/tsv/contents_blastp_genus.tsv", sep='\t', index=False)

    visualize_heatmap(output_dir, "contents_blastp_species2", 0)
    visualize_heatmap(output_dir, "contents_blastp_genus", 1)

#===============================================================================
def count_mmseq_result(project_info, TAXON):
    print("<< counting target by species and genus level...")
    input_dir = project_info['input']
    output_dir = project_info['output']
    df = pd.read_csv(f"{project_info['data']}/target.tsv", sep='\t')
    target_names = df['Prediction'].drop_duplicates().tolist()

    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    start_time = time.time()

    target_cnt_list_1 = []
    for genus in all_directories:
        genus_df_list = []
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GC" in accession:
                continue
            file_path = f"{input_dir}/{genus}/{accession}/target_seq.tsv"
            tmp = count_mmseq_species(file_path, target_names)
            genus_df_list.append(tmp)
        
        if genus_df_list:
            genus_df = pd.concat(genus_df_list, ignore_index=True)
            target_cnt_list_1.append(genus_df)

    target_cnt_1 = pd.concat(target_cnt_list_1, ignore_index=True)
    target_cnt_1 = target_cnt_1.sort_values("name")
    target_cnt_1.to_csv(f"{output_dir}/tsv/contents_mmseq_species.tsv", sep='\t', index=False)

    count_mmseq_genus(output_dir, target_cnt_1)
    visualize_heatmap(output_dir, "contents_mmseq_species", 0)
    visualize_heatmap(output_dir, "contents_mmseq_genus", 1)
    target_sum = target_cnt_1[target_names]
    target_sum = target_sum.sum(axis=0).reset_index()
    target_sum.to_csv(f"{output_dir}/tsv/target_frequency.tsv", sep='\t', index=False)
    visualize_freq(output_dir, TAXON)

    target_cnt_list_2 = []
    for genus in all_directories:
        genus_df_list = []
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GC" in accession:
                continue
            file_path = f"{input_dir}/{genus}/{accession}/target_merge.tsv"
            tmp = count_mmseq_species(file_path, target_names)
            genus_df_list.append(tmp)
        
        if genus_df_list:
            genus_df = pd.concat(genus_df_list, ignore_index=True)
            target_cnt_list_2.append(genus_df)

    target_cnt_2 = pd.concat(target_cnt_list_2, ignore_index=True)
    target_cnt_2 = target_cnt_2.sort_values("name")
    target_cnt_2.to_csv(f"{output_dir}/tsv/contents_merge_species.tsv", sep='\t', index=False)

    count_merge_genus(output_dir, target_cnt_2)
    visualize_heatmap(output_dir, "contents_merge_species", 0)
    visualize_heatmap(output_dir, "contents_merge_genus", 1)

    end_time = time.time()
    total = end_time - start_time
    print(f"   └── counting done (elapsed time: {round(total / 60, 3)} min)")

def count_mmseq_species(file_path, target_names):
    with open(file_path, 'r') as f:
        accession = f.readline().split('#')[-1].strip()
        name = f.readline().split('#')[-1].strip()

    df = pd.read_csv(file_path, sep='\t', comment='#')

    target_cnt = df['Prediction'].value_counts().reset_index()
    target_cnt.columns = ['protein', 'count']

    target_cnt = pd.DataFrame({'protein': target_names}).merge(target_cnt, on='protein', how='left').fillna(0)
    target_cnt = target_cnt.astype({'count': int})

    result_df = pd.DataFrame(columns=['accession', 'name'] + target_names)
    result_df.loc[0, 'accession'] = accession
    result_df.loc[0, 'name'] = name
    result_df.loc[0, target_cnt['protein']] = target_cnt['count'].values

    return result_df

def count_mmseq_genus(output_dir, df):
  df['genus'] = df['name'].str.split().str[0]

  protein_columns = df.columns[2:-1]
  result_df = df.drop('name', axis=1).groupby('genus')[protein_columns].mean().reset_index()
  result_df = result_df.fillna(0)
  result_df = result_df.round(3)
  result_df.to_csv(f"{output_dir}/tsv/contents_mmseq_genus.tsv", sep='\t', index=False)

def count_merge_genus(output_dir, df):
  df['genus'] = df['name'].str.split().str[0]

  protein_columns = df.columns[2:-1]
  result_df = df.drop('name', axis=1).groupby('genus')[protein_columns].mean().reset_index()
  result_df = result_df.fillna(0)
  result_df = result_df.round(3)
  result_df.to_csv(f"{output_dir}/tsv/contents_merge_genus.tsv", sep='\t', index=False)

#===============================================================================
def count_mmseq_cluster(project_info):
    print("<< counting target by species and genus level...")
    input_dir = project_info['input']
    output_dir = project_info['output']

    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    seq_lib = pd.read_csv(f"{project_info['seqlib']}/seqlib.tsv", sep='\t')
    seq_lib['Prediction'] = seq_lib['Prediction'].astype(str)
    seq_lib['RepSeq'] = seq_lib['RepSeq'].astype(str)
    seq_lib['Total'] = seq_lib['Total'].astype(int)
    seq_lib['SeqList'] = seq_lib['SeqList'].astype(str)

    seq_lib['no'] = seq_lib.groupby('Prediction').cumcount() + 1
    seq_lib['Prediction'] = seq_lib['Prediction'] + '-' + seq_lib['no'].astype(str)
    seq_lib.drop(['no', 'Total'], axis=1, inplace=True)

    genus_rep = seq_lib.copy()

    start_time = time.time()

    for genus in all_directories:
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GC" in accession:
                continue
            file_path = f"{input_dir}/{genus}/{accession}/target_seq.tsv"
        
            data_df = pd.read_csv(file_path, sep='\t', comment='#')
            seq_lib = pd.concat([seq_lib, pd.DataFrame({accession: [0] * len(seq_lib)})], axis=1)
            seq_lib.loc[seq_lib['RepSeq'].isin(data_df['repseq']), accession] = 1

            genus_rep = pd.concat([genus_rep, pd.DataFrame({genus: [0] * len(genus_rep)})], axis=1)
            genus_rep.loc[genus_rep['RepSeq'].isin(data_df['repseq']), genus] = 1
    				
        seq_lib = pd.concat([seq_lib], axis=1)
        genus_rep = pd.concat([genus_rep], axis=1)

    seq_lib.drop(['RepSeq', 'SeqList'], axis=1, inplace=True)

    seq_lib = seq_lib.transpose()
    seq_lib.rename(columns=seq_lib.iloc[0], inplace=True)
    seq_lib = seq_lib.reset_index()

    seq_lib = seq_lib[1:]
    seq_lib.to_csv(f"{output_dir}/tsv/repseq_species.tsv", sep='\t', index=False)

    genus_rep.drop(['RepSeq', 'SeqList'], axis=1, inplace=True)
    genus_rep = genus_rep.transpose()
    genus_rep.rename(columns=genus_rep.iloc[0], inplace=True)
    genus_rep = genus_rep.reset_index()
    genus_rep = genus_rep[1:]
    genus_rep = genus_rep.groupby('index').mean().reset_index()
    genus_rep.to_csv(f"{output_dir}/tsv/repseq_genus.tsv", sep='\t', index=False)    

    end_time = time.time()
    total = end_time - start_time
    print(f"   └── counting done (elapsed time: {round(total / 60, 3)} min)")
