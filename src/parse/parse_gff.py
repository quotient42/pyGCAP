import pandas as pd

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.max_colwidth', None)

def process_line(line, key_count, summarized_list):

    if line.startswith("#"):
        return key_count, summarized_list

    fields = line.strip().split("\t")
    key_count += 1

    attributes = fields[-1]
    attr_dict = dict(item.split('=', 1) for item in attributes.strip(' []').split(';'))
    protein_id_value = attr_dict.get('protein_id', '-')

    summarized_dict = {
        'contig': fields[0],
        # 'type': fields[1],
        'start': fields[3],
        'end': fields[4],
        # 'score': fields[5],
        'strand': fields[6],
        # 'phase': fields[7],
        'protein_id': protein_id_value,
    }
    summarized_list[key_count] = summarized_dict

    return key_count, summarized_list

def read_and_process_file(gff_input):
    key_count = 0
    summarized_list = {}

    with open(gff_input, 'r') as file:
        for line in file:
            key_count, summarized_list = process_line(line, key_count, summarized_list)

    return summarized_list

def parse_gff_main(gff_input):
    data_dict = read_and_process_file(gff_input)

    summarized_list = list(data_dict.values())
    df = pd.DataFrame(summarized_list)
    filtered_df = df[df['protein_id'] != '-']
    
    return filtered_df 
