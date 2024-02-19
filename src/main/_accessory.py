import os
import pandas as pd

def collect_accessory(project_info):
    print("<< searching target based on seq data...")
    input_dir = project_info['input']
    data_dir = project_info['data']
    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    seq_df = pd.read_csv(f"{data_dir}/result_cluster.tsv", sep="\t").copy()
    seq_df['subseq'] = seq_df['subseq'].apply(lambda x: float(x.replace('WP_', '')))
    seq_df = seq_df.sort_values('subseq')

    seq_lib = pd.read_csv(f"{project_info['seqlib']}/seqlib.tsv", sep='\t')
    seq_lib['Prediction'] = seq_lib['Prediction'].astype(str)
    seq_lib['RepSeq'] = seq_lib['RepSeq'].astype(str)
    seq_lib['Total'] = seq_lib['Total'].astype(int)
    seq_lib['Seq_List'] = seq_lib['Seq_List'].astype(str)

    all_data = pd.DataFrame()

    for genus in all_directories:
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GCF" in accession:
                continue
            
            cur_dir = f"{input_dir}/{genus}/{accession}"

            with open(f"{cur_dir}/genome_summary.tsv", 'r', encoding='utf-8') as file:
                lines = file.readlines()
                species = lines[1][1:].strip()

            genome_summary = pd.read_csv(f"{cur_dir}/genome_summary.tsv", sep='\t', comment='#')
            cluster = pd.read_csv(f"{cur_dir}/target_cluster.tsv", sep='\t', comment='#')
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

    all_data['Accession_tmp'] = all_data['Accession'].apply(lambda x: float(x.replace('WP_', '')))

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
    all_data.to_csv(f"{project_info['output']}/tsv/accessory.tsv", sep='\t', index=False)

    repseq_counts = all_data['repseq'].value_counts()
    repseq_counts_df = repseq_counts.reset_index()
    repseq_counts_df.columns = ['repseq', 'count']

    # 새로운 DataFrame을 만들어서 Prediction을 추가합니다.
    merged_df = pd.merge(repseq_counts_df, seq_lib, how='left', left_on='repseq', right_on='RepSeq')

    # 필요한 column만 선택합니다.
    merged_df = merged_df[['repseq', 'count', 'Prediction']]

    merged_df.to_csv(f"{project_info['output']}/tsv/accessory_rep.tsv", sep='\t', index=False)

