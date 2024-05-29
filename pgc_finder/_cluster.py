'''
    cluster_target  ˑˑˑˑˑˑˑˑˑˑˑˑˑ  to cluster & visualize target gene cluster
    ├── search_cluster  ˑˑˑˑˑˑˑˑˑ  to search gene cluster in each species data
    └── adjust_cluster  ˑˑˑˑˑˑˑˑˑ  to adjust cluster data for visualization
    classify_cluster_type  ˑˑˑˑˑˑ  to classify clusters into three type (C/F/D)

'''

import os
import pandas as pd
import time
import random

from ._visualize import visualize_cluster
from ._visualize import visualize_cluster_type

#===============================================================================
def cluster_target(project_info):
    print("<< clustering target genes...")
    df = pd.read_csv(f"{project_info['data']}/target.tsv", sep='\t')
    target_list = df['Prediction'].drop_duplicates().tolist()
    init_color_dict(target_list)

    input_dir = project_info['input']
    input_len = project_info['input_len']
    output_dir = project_info['output']

    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    start_time = time.time()

    for genus in all_directories:
        if ("output" in genus) or ("seqlib" in genus):
            continue
        laps_start = time.time()
        path = os.path.join(input_dir, genus)
        if not os.path.exists(f"{path}/tmp"):
            os.makedirs(f"{path}/tmp")

        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GC" in accession:
                continue
            cur_dir = f"{input_dir}/{genus}/{accession}"

            with open(f"{cur_dir}/genome_summary.tsv", 'r') as file:
                lines = file.readlines()
                species = lines[1][1:].strip()

            genome_summary = pd.read_csv(f"{cur_dir}/genome_summary.tsv", sep='\t', comment='#')
            target_list = pd.read_csv(f"{cur_dir}/target_merge.tsv", sep='\t', comment='#')
            input_info = [accession, cur_dir, species]

            cluster_summary = search_cluster(input_info, genome_summary, target_list)
            adjust_cluster(input_info, cluster_summary, f"{path}/tmp")

        visualize_cluster(f"{path}/tmp", f"{output_dir}/genus", genus, color_dict)

    end_time = time.time()
    total = end_time - start_time
    print(f"   └── clustering done (elapsed time: {round(total / 60, 3)} min)")
    
#===============================================================================
def search_cluster(input_info, genome_summary, target_summary):
    df_gff = genome_summary
    df_dcw = target_summary

    df_gff['Prediction'] = '----'
    df_gff['Prediction'] = df_gff.apply(_update_prediction, df_dcw=df_dcw, axis=1).fillna('----')

    contig_groups = df_gff[df_gff['Prediction'] != '----'].groupby('contig').filter(lambda x: len(x) >= 2)['contig'].unique()
    df_cluster = df_gff[df_gff['contig'].isin(contig_groups)].reset_index(drop=True)

    df_cluster = _check_contig(contig_groups, df_cluster)
    df_cluster = _check_block(contig_groups, df_cluster)
    
    key_order = ["contig", "protein_id", "Prediction", "cluster", "block", "strand", "start", "end"]
    df_cluster = df_cluster[key_order]

    cluster_summary = f"{input_info[1]}/target_cluster.tsv"
    
    with open(cluster_summary, 'w') as summary_file:
        summary_file.write(f"# {input_info[0]}\n# {input_info[2]}\n")
        df_cluster.to_csv(summary_file, sep='\t', index=False)

    return cluster_summary

def _check_contig(contig_groups, df_cluster):
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

def _check_block(contig_groups, df_cluster):
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

def _update_prediction(row, df_dcw):
    protein_id = row['protein_id']
    matching_rows = df_dcw[df_dcw['protein_id'] == protein_id]['Prediction']
    prediction = matching_rows.iloc[0] if not matching_rows.empty else '----'
    return prediction

#===============================================================================
def adjust_cluster(input_info, cluster_summary, fig_input_dir):
    df = pd.read_csv(cluster_summary, sep='\t', comment='#')
    modified_df = _process_cluster_data(df)

    if modified_df.empty:
        return pd.DataFrame()

    modified_groups2 = []
    grouped_df = modified_df.groupby('block')
    for group_name, group_data in grouped_df:
        modified_groups2.append(_fill_blank_rows(group_data))

    modified_df = pd.concat(modified_groups2, ignore_index=True)
    modified_df = modified_df[['block', 'Prediction', 'start', 'end']]
    
    filename = f"{fig_input_dir}/{input_info[0]}.tsv"
    modified_df.to_csv(filename, sep='\t', index=False)

def _process_cluster_data(df):
    df = df[['Prediction', 'block', 'start', 'end']]
    df = df[df['block'] > 0]

    modified_groups = []
    group_end = 0

    grouped_df = df.groupby('block')
    if len(grouped_df) > 0:
        modified_groups = []
        group_end = 0

        for group_name, group_data in grouped_df:
            first_row = group_data.iloc[0]
            block_start = first_row['start']

            group_data['start'] = (group_data['start'] - block_start) // 100 + group_end
            group_data['end'] = (group_data['end'] - block_start) // 100 + group_end

            group_end = group_data.iloc[-1]['end'] + 10

            modified_groups.append(group_data.copy())
        modified_df = pd.concat(modified_groups, ignore_index=True)
    else:
        df.loc[0] = ["----", 0, 0, 0]
        modified_df = df

    return modified_df

def _fill_blank_rows(group_data):
    new_rows = []

    for index in range(len(group_data) - 1):
        current_row = group_data.iloc[index]
        next_row = group_data.iloc[index + 1]

        new_rows.append(current_row.to_dict())
        if current_row['end'] != next_row['start']:
            new_row = {
                'block': current_row['block'],
                'Prediction': 'blank',
                'start': current_row['end'],
                'end': next_row['start']
            }
            new_rows.append(new_row)

    new_rows.append(group_data.iloc[-1].to_dict())

    return pd.DataFrame(new_rows)

#===============================================================================
def classify_cluster_type(project_info):
	print("<< classifying target gene cluster...")
	input_dir = project_info['input']
	output_dir = project_info['output']
	all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]

	start_time = time.time()

	cluster_df = pd.read_csv(f"{output_dir}/tsv/contents_mmseq_species.tsv", sep='\t', comment='#').copy()
	cluster_df.insert(1, 'cluster', "-")

	for genus in all_directories:
			for accession in os.listdir(os.path.join(input_dir, genus)):
					if not "GC" in accession:
							continue
					
					filename = f"{input_dir}/{genus}/{accession}/target_cluster.tsv"
					df = pd.read_csv(filename, sep='\t', comment='#')

					protein_cnt = len(df[(df['Prediction'] != '----') & (df['cluster'] == True)])
					df = df[df['block'] > 0]
					block_sizes = []
					for block_num in df['block'].unique():
							Predictions_in_block = df[df['block'] == block_num]['Prediction'].tolist()
							non_dash_count = sum(1 for pred in Predictions_in_block if pred != '----')
							block_sizes.extend([non_dash_count] * len(Predictions_in_block))

					df['block_size'] = block_sizes
					block_cnt = df['block_size'].max()

					if protein_cnt < 9:
						result = "D"
					elif block_cnt >= 9:
						result = "C"
					else:
						result = "F"
					
					cluster_df.loc[cluster_df['accession'] == accession, 'cluster'] = result

	split_names = cluster_df['name'].str.split(' ', n=1, expand=True)
	cluster_df['genus'] = split_names[0]
	cluster_df['species'] = split_names[1]

	cluster_df = cluster_df[['genus', 'species', 'cluster']]
	cluster_df.to_csv(f"{output_dir}/tsv/cluster_type_species.tsv", sep='\t', index=False)

	conserved_count = cluster_df[cluster_df['cluster'] == 'C'].groupby('genus').size()
	fragmented_count = cluster_df[cluster_df['cluster'] == 'F'].groupby('genus').size()
	disrupted_count = cluster_df[cluster_df['cluster'] == 'D'].groupby('genus').size()

	cluster_df['conserved'] = cluster_df['genus'].map(conserved_count).fillna(0).astype(int)
	cluster_df['fragmented'] = cluster_df['genus'].map(fragmented_count).fillna(0).astype(int)
	cluster_df['disrupted'] = cluster_df['genus'].map(disrupted_count).fillna(0).astype(int)
	cluster_df = cluster_df[['genus', 'conserved', 'fragmented', 'disrupted']]
	unique_genus_df = cluster_df.groupby('genus').first().reset_index()
	total_counts = unique_genus_df[['conserved', 'fragmented', 'disrupted']].sum(axis=1)

	unique_genus_df['c_ratio'] = round(unique_genus_df['conserved'] / total_counts, 2)
	unique_genus_df['f_ratio'] = round(unique_genus_df['fragmented'] / total_counts, 2)
	unique_genus_df['d_ratio'] = round(unique_genus_df['disrupted'] / total_counts, 2)

	unique_genus_df.fillna(0, inplace=True)

	unique_genus_df.to_csv(f"{output_dir}/tsv/cluster_type_genus.tsv", sep='\t', index=False)

	visualize_cluster_type(output_dir, unique_genus_df)

	end_time = time.time()
	total = end_time - start_time
	print(f"   └── classification done (elapsed time: {round(total / 60, 3)} min)")

#===============================================================================
color_dict = {
    "color_1": "#8ED3C7",
    "color_2": "#F99FB5",
    "color_3": "#2A9D8F",
    "color_4": "#BEBADB",
    "color_5": "#FB7F72",
    "color_6": "#81B1D3",
    "color_7": "#FDB364",
    "color_8": "#B2DE68",
    "color_9": "#E4B7CD",
    "color_10": "#B6C8FE",
    "color_11": "#BD7FBB",
    "color_12": "#CDEBC3",
    "color_13": "#DFC27C",
    "color_14": "#C7E8FA",
    "color_15": "#A95C68",
    "color_16": "#FCBBA2",
    "color_17": "#F7DA84",
    "blank": "#F0F0F0",
    "etc": "#C0C0C0"
}

def is_gray_scale(hex_color):
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return r == g == b

def generate_unique_hex_color():
    while True:
        hex_color = ''.join([random.choice('0123456789ABCDEF') for _ in range(6)])
        if not is_gray_scale(hex_color) and hex_color not in color_dict.values():
            hex_color = "#" + str(hex_color)
            return hex_color

def init_color_dict(key_list):
    global color_dict
    existing_keys = [key for key in color_dict.keys() if key not in ['blank', 'etc']]
    existing_colors = [color_dict[key] for key in existing_keys]

    for i, key in enumerate(key_list):
        if i < 17:
            color_dict[key] = existing_colors[i]
        else:
            unique_color = generate_unique_hex_color()
            color_dict[key] = unique_color

    keys_to_delete = [key for key in color_dict if key.startswith('color_') and key not in key_list]
    for key in keys_to_delete:
        del color_dict[key]
    