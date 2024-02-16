def check_contig(contig_groups, df_cluster):
    df_cluster['cluster'] = False

    for contig_group in contig_groups:
        contig_group_rows = df_cluster[df_cluster['contig'] == contig_group]

        cur_index = 0
        while cur_index < len(contig_group_rows):
            cur_row = contig_group_rows.iloc[cur_index]

            if cur_row['Prediction'] != '----':
                neighbor_index = cur_index + 1
                while neighbor_index < len(contig_group_rows) and contig_group_rows.iloc[neighbor_index]['Prediction'] == '----':
                    neighbor_index += 1

                if neighbor_index < len(contig_group_rows):
                    neighbor_row = contig_group_rows.iloc[neighbor_index]

                    if neighbor_row.name - cur_row.name <= 5:
                        df_cluster.loc[cur_row.name:neighbor_row.name, 'cluster'] = True

                    cur_index = neighbor_index
                else:
                    cur_index += 1
            else:
                cur_index += 1

    return df_cluster

def check_block(contig_groups, df_cluster):
    df_cluster['block'] = 0
    block_cnt = 0

    for contig_group in contig_groups:
        contig_group_rows = df_cluster[df_cluster['contig'] == contig_group].copy()

        cur_index = 0
        while cur_index < len(contig_group_rows):
            cur_row = contig_group_rows.iloc[cur_index]

            if cur_index > 0:
                prev_row = contig_group_rows.iloc[cur_index - 1]
                if prev_row['cluster'] == False and cur_row['cluster'] == True:
                    block_cnt += 1
                    while cur_index < len(contig_group_rows) and cur_row['cluster']:
                        df_cluster.at[cur_row.name, 'block'] = block_cnt
                        cur_index += 1
                        if cur_index < len(contig_group_rows):
                            cur_row = contig_group_rows.iloc[cur_index]
            else:
                if cur_row['cluster'] == True:
                    block_cnt += 1
                    while cur_index < len(contig_group_rows) and cur_row['cluster']:
                        df_cluster.at[cur_row.name, 'block'] = block_cnt
                        cur_index += 1
                        if cur_index < len(contig_group_rows):
                            cur_row = contig_group_rows.iloc[cur_index]

            cur_index += 1

    df_cluster = df_cluster[(df_cluster['Prediction'] != '----') | (df_cluster['cluster'])]

    return df_cluster

def search_cluster_main(input_info, genome_summary, target_summary):
    df_gff = genome_summary
    df_dcw = target_summary

    df_gff['Prediction'] = '----'
    df_gff['Prediction'] = df_gff['protein_id'].map(df_dcw.set_index('protein_id')['Prediction']).fillna('----')

    contig_groups = df_gff[df_gff['Prediction'] != '----'].groupby('contig').filter(lambda x: len(x) >= 2)['contig'].unique()
    df_cluster = df_gff[df_gff['contig'].isin(contig_groups)].reset_index(drop=True)

    df_cluster = check_contig(contig_groups, df_cluster)
    df_cluster = check_block(contig_groups, df_cluster)
    
    key_order = ["contig", "protein_id", "Prediction", "cluster", "block", "strand", "start", "end"]
    df_cluster = df_cluster[key_order]

    cluster_summary = f"{input_info[1]}/target_cluster.tsv"
    
    with open(cluster_summary, 'w') as summary_file:
        summary_file.write(f"# {input_info[0]}\n# {input_info[2]}\n")
        df_cluster.to_csv(summary_file, sep='\t', index=False)

    return cluster_summary
