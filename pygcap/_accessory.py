import os
import time
import pandas as pd
import ast
from ._utils import remove_non_numeric

def collect_accessory(project_info):
    print("<< collecting accessory genes...")
    input_dir = project_info['input']
    output_dir = project_info['output']
    seqlib_dir = project_info['seqlib']
    data_dir = project_info['data']
    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    seq_df = pd.read_csv(f"{data_dir}/result_cluster.tsv", sep="\t").copy()
    seq_df['subseq'] = seq_df['subseq'].apply(lambda x: float(remove_non_numeric(x)))
    seq_df = seq_df.sort_values('subseq')

    seq_lib = pd.read_csv(f"{project_info['seqlib']}/seqlib.tsv", sep='\t')
    seq_lib['Prediction'] = seq_lib['Prediction'].astype(str)
    seq_lib['RepSeq'] = seq_lib['RepSeq'].astype(str)
    seq_lib['Total'] = seq_lib['Total'].astype(int)
    seq_lib['SeqList'] = seq_lib['SeqList'].astype(str)

    all_data = pd.DataFrame()

    for genus in all_directories:
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GC" in accession:
                continue
            
            cur_dir = f"{input_dir}/{genus}/{accession}"

            with open(f"{cur_dir}/genome_summary.tsv", 'r', encoding='utf-8') as file:
                lines = file.readlines()
                species = lines[1][1:].strip()

            genome_summary = pd.read_csv(f"{cur_dir}/genome_summary.tsv", sep='\t', comment='#')
            cluster = pd.read_csv(f"{cur_dir}/probe_cluster.tsv", sep='\t', comment='#')
            input_info = [accession, cur_dir, species]

            df = cluster.copy()
            df = df[df['cluster']]

            genome_summary_subset = genome_summary[['protein_id', 'gene', 'product', 'translation']]

            merged_df = pd.merge(df, genome_summary_subset, how='left', on='protein_id')
            merged_df['source'] = species

            merged_df.to_csv(f"{cur_dir}/accessory1.tsv", sep='\t', index=False)

            merged_df = merged_df[merged_df['Prediction'] == "----"]
            merged_df = merged_df[['source', 'protein_id', 'gene', 'product', 'translation']]
            merged_df.columns = ['source', 'Accession', 'gene', 'product', 'translation']
            merged_df.to_csv(f"{cur_dir}/accessory2.tsv", sep='\t', index=False)

            all_data = pd.concat([all_data, merged_df])

    all_data['Accession_tmp'] = all_data['Accession'].apply(lambda x: float(remove_non_numeric(x)))

    all_data['repseq'] = "-"

    for index, row in all_data.iterrows():
            accession = row['Accession_tmp']
            low = 0
            high = len(seq_df) - 1
            while low <= high:
                    mid = (low + high) // 2
                    if seq_df.iloc[mid]['subseq'] == accession:
                            all_data.at[index, 'repseq'] = seq_df.iloc[mid]['repseq']
                            break
                    elif seq_df.iloc[mid]['subseq'] < accession:
                            low = mid + 1
                    else:
                            high = mid - 1

    all_data = all_data[all_data['repseq'] != '-']
    all_data = all_data[['source', 'Accession', 'repseq', 'gene', 'product', 'translation']]
    all_data.to_csv(f"{output_dir}/tsv/accessory.tsv", sep='\t', index=False)

    repseq_counts = all_data['repseq'].value_counts()
    repseq_counts_df = repseq_counts.reset_index()
    
    unique_repseq = all_data.drop_duplicates('repseq')
    repseq_counts_df['product'] = repseq_counts_df['repseq'].map(unique_repseq.set_index('repseq')['product'])

    repseq_counts_df.columns = ['repseq', 'count', 'product']

    merged_df = pd.merge(repseq_counts_df, seq_lib, how='left', left_on='repseq', right_on='RepSeq')

    merged_df['RepName'] = [f"Acc.{i}" for i in range(1, len(merged_df) + 1)]
    merged_df = merged_df[['RepName', 'repseq', 'count', 'product']]

    merged_df.to_csv(f"{output_dir}/tsv/accessory_rep.tsv", sep='\t', index=False)

    # construct seqlib_acc
    mmseq_df = pd.read_csv(f"{seqlib_dir}/mmseq_rep.tsv", sep='\t')
    merged_df = pd.merge(merged_df, mmseq_df, left_on='repseq', right_on='RepSeq', how='left')

    # 필요한 열만 선택
    merged_df = merged_df[['RepName', 'repseq', 'count', 'product', 'SeqList']]
    merged_df = merged_df.rename(columns={'repseq': 'RepSeq', 'count': 'Total'})

    # 저장
    merged_df.to_csv(f"{output_dir}/tsv/acclib.tsv", sep='\t', index=False)

    print("<< searching accessory genes based on acclib...")
    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    acc_lib = pd.read_csv(f"{output_dir}/tsv/acclib.tsv", sep='\t')
    acc_lib['RepName'] = acc_lib['RepName'].astype(str)
    acc_lib['RepSeq'] = acc_lib['RepSeq'].astype(str)
    acc_lib['Total'] = acc_lib['Total'].astype(int)
    acc_lib['SeqList'] = acc_lib['SeqList'].astype(str)

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
            probe_acc = search_by_species(input_info, genome_summary, acc_lib)
            probe_seq = pd.read_csv(f"{cur_dir}/probe_seq.tsv", sep='\t', comment='#')

            common_columns = [col for col in probe_acc.columns if col != 'Prediction']
            probe_final = pd.merge(probe_seq[common_columns], probe_acc[common_columns], how='outer')

            with open(f"{cur_dir}/probe_final.tsv", 'w') as summary_file:
                summary_file.write(f"# {input_info[0]}\n# {input_info[2]}\n")
                probe_final.to_csv(summary_file, sep='\t', index=False)

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
    merged_df = merged_df.drop_duplicates(subset=['RepName', 'contig', 'protein_id', 'strand', 'start', 'end'])

    merged_df['repseq'] = merged_df['RepSeq'].str.strip()

    hit_df = merged_df[['RepName', 'contig', 'protein_id', 'strand', 'start', 'end', 'repseq']]
    hit_df = hit_df.sort_values(by='RepName')

    probe_acc = f"{input_info[1]}/probe_acc.tsv"

    with open(probe_acc, 'w') as summary_file:
        summary_file.write(f"# {input_info[0]}\n# {input_info[2]}\n")
        hit_df.to_csv(summary_file, sep='\t', index=False)

    return hit_df

def _convert_to_list(value):
    try:
        return ast.literal_eval(value) if isinstance(value, str) else value
    except (SyntaxError, ValueError):
        return value
