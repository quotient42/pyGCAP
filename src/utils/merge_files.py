import pandas as pd

def merge_gbff_gff(input_info, gbff_summary, gff_summary):
    merged_df = pd.merge(gbff_summary, gff_summary, on='protein_id')
    merged_df = merged_df[['protein_id', 'contig', 'strand', 'start', 'end', 'gene', 'EC_number', 'product', 'translation']]

    genome_summary = f"{input_info[1]}/genome_summary.tsv"

    with open(genome_summary, 'w') as summary_file:
        summary_file.write(f"# {input_info[0]}\n# {input_info[2]}\n")
        merged_df.to_csv(summary_file, sep='\t', index=False)

    return merged_df
