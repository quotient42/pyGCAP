import pandas as pd
import ast

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.max_colwidth', None)

def convert_to_list(value):
    try:
        return ast.literal_eval(value) if isinstance(value, str) else value
    except (SyntaxError, ValueError):
        return value

def search_target_by_seq_main(input_info, genome_summary, seq_lib):
    df_gff = genome_summary.copy()

    seq_lib['Seq_List'] = seq_lib['Seq_List'].apply(convert_to_list)

    seq_set = set(seq_lib['Seq_List'].explode())

    df_gff['protein_id'] = df_gff['protein_id'].str.strip()

    mask = df_gff['protein_id'].isin(seq_set)

    seq_lib_exploded = seq_lib.explode('Seq_List')
    seq_lib_exploded = seq_lib_exploded.reset_index(drop=True)

    merged_df = df_gff[mask].merge(seq_lib_exploded, left_on='protein_id', right_on='Seq_List', how='inner')
    merged_df = merged_df.drop_duplicates(subset=['prediction', 'contig', 'protein_id', 'strand', 'start', 'end'])

    merged_df['repseq'] = merged_df['RepSeq'].str.strip()

    hit_df = merged_df[['prediction', 'contig', 'protein_id', 'strand', 'start', 'end', 'repseq']]
    hit_df = hit_df.sort_values(by='prediction')

    target_seq = f"{input_info[1]}/target_seq.tsv"

    with open(target_seq, 'w') as summary_file:
        summary_file.write(f"# {input_info[0]}\n# {input_info[2]}\n")
        hit_df.to_csv(summary_file, sep='\t', index=False)

    return hit_df
