import pandas as pd

def process_data(df):
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

def fill_blank_rows(group_data):
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

def process_and_concat_data(file_path):
    df = pd.read_csv(file_path, sep='\t', comment='#')
    modified_df = process_data(df)

    if modified_df.empty:
        return pd.DataFrame()

    modified_groups2 = []
    grouped_df = modified_df.groupby('block')
    for group_name, group_data in grouped_df:
        modified_groups2.append(fill_blank_rows(group_data))

    modified_df = pd.concat(modified_groups2, ignore_index=True)
    modified_df = modified_df[['block', 'Prediction', 'start', 'end']]

    return modified_df

def adjust_cluster_main(input_info, cluster_summary, fig_input_dir):
    result_df = process_and_concat_data(cluster_summary)
    
    filename = f"{fig_input_dir}/{input_info[0]}.tsv"
    result_df.to_csv(filename, sep='\t', index=False)
