'''
    search_target  
    └── search_by_species  ˑˑˑˑˑˑ  search target in each species using seqlib

'''

import os
import pandas as pd
import time
import ast

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.max_colwidth', None)

#==============================================================================
def search_target(project_info):
    print("<< searching target based on seqlib...")
    input_dir = project_info['input']
    input_len = project_info['input_len']
    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    seq_lib = pd.read_csv(f"{project_info['seqlib']}/seqlib.tsv", sep='\t')
    seq_lib['RepName'] = seq_lib['RepName'].astype(str)
    seq_lib['Prediction'] = seq_lib['Prediction'].astype(str)
    seq_lib['RepSeq'] = seq_lib['RepSeq'].astype(str)
    seq_lib['Total'] = seq_lib['Total'].astype(int)
    seq_lib['SeqList'] = seq_lib['SeqList'].astype(str)

    progress = 0
    start_time = time.time()

    for genus in all_directories:
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GC" in accession:
                continue
            
            cur_dir = f"{input_dir}/{genus}/{accession}"

            with open(f"{cur_dir}/genome_summary.tsv", 'r', encoding='utf-8') as file:
                lines = file.readlines()
                species = lines[1][1:].strip()

            genome_summary = pd.read_csv(f"{cur_dir}/genome_summary.tsv", sep='\t', comment='#')
            input_info = [accession, cur_dir, species]
            target_seq = search_by_species(input_info, genome_summary, seq_lib)  # target_seq.tsv

            target_blast = pd.read_csv(f"{cur_dir}/target_blast2.tsv", sep='\t', comment='#')

            common_columns = [col for col in target_seq.columns if col not in ['RepName', 'repseq']]
            merged_df = pd.merge(target_blast, target_seq, on=common_columns, how='outer')
            merged_df = merged_df.fillna('------')
            merged_df = merged_df.sort_values('Prediction')

            first_cols = ['Prediction', 'TarName', 'RepName', 'repseq']
            remaining_cols = [col for col in merged_df.columns if col not in first_cols]

            new_order = first_cols + remaining_cols
            merged_df = merged_df[new_order]

            with open(f"{cur_dir}/target_merge.tsv", 'w') as summary_file:
                summary_file.write(f"# {input_info[0]}\n# {input_info[2]}\n")
                merged_df.to_csv(summary_file, sep='\t', index=False)

    end_time = time.time()
    total = end_time - start_time
    print(f"   └── searching done (elapsed time: {round(total / 60, 3)} min)")

def search_by_species(input_info, genome_summary, seq_lib):
    df_gff = genome_summary.copy()

    seq_lib['SeqList'] = seq_lib['SeqList'].apply(_convert_to_list)

    seq_set = set(seq_lib['SeqList'].explode())

    df_gff['protein_id'] = df_gff['protein_id'].str.strip()

    mask = df_gff['protein_id'].isin(seq_set)

    seq_lib_exploded = seq_lib.explode('SeqList')
    seq_lib_exploded = seq_lib_exploded.reset_index(drop=True)

    merged_df = df_gff[mask].merge(seq_lib_exploded, left_on='protein_id', right_on='SeqList', how='inner')
    merged_df = merged_df.drop_duplicates(subset=['RepName', 'Prediction', 'contig', 'protein_id', 'strand', 'start', 'end'])

    merged_df['repseq'] = merged_df['RepSeq'].str.strip()

    hit_df = merged_df[['RepName', 'Prediction', 'contig', 'protein_id', 'strand', 'start', 'end', 'repseq']]
    hit_df = hit_df.sort_values(by='RepName')

    target_seq = f"{input_info[1]}/target_seq.tsv"

    with open(target_seq, 'w') as summary_file:
        summary_file.write(f"# {input_info[0]}\n# {input_info[2]}\n")
        hit_df.to_csv(summary_file, sep='\t', index=False)

    return hit_df

def _convert_to_list(value):
    try:
        return ast.literal_eval(value) if isinstance(value, str) else value
    except (SyntaxError, ValueError):
        return value
