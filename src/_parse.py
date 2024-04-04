''' 
    parse_genome
    ├── parse_gbff  ˑˑˑˑˑˑˑˑˑˑ  to parse gbff file.
    ├── parse_gff  ˑˑˑˑˑˑˑˑˑˑˑ  to parse gff file.
    └── merge_gbff_gff  ˑˑˑˑˑˑ  to merge parsed genome data (gbff + gff).

'''
 
import os
import time

import sys
sys.path.append("./info/")
from gbff_info import gbff_key_setting

import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.max_colwidth', None)

#===== parse_genome =================================================================
def parse_genome(base_info):
    print("<< parsing genome data...")
    base_dir_path = base_info['input']
    base_dir_len = base_info['input_len']
    all_directories = [d for d in os.listdir(base_dir_path) if os.path.isdir(os.path.join(base_dir_path, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    progress = 0
    start_time = time.time()

    for genus in all_directories:
        laps_start = time.time()
        for accession in os.listdir(os.path.join(base_dir_path, genus)):
            if accession == "tmp":
                continue
            progress += 1
            print(f"   ├── processing... ({progress}/{base_dir_len}): {accession}") 
            
            cur_dir = f"{base_dir_path}/{genus}/{accession}"
            gbff_input = f"{cur_dir}/genomic.gbff"
            gff_input = f"{cur_dir}/genomic.gff"
            
            input_info = [accession, cur_dir]

            species, gbff_summary = parse_gbff(gbff_input)
            input_info.append(species)
            gff_summary = parse_gff(gff_input)
            genome_summary = merge_gbff_gff(input_info, gbff_summary, gff_summary)

    end_time = time.time()
    total = end_time - start_time
    print(f"   └── parsing done (elapsed time: {round(total / 60, 3)} min)")
  
#===== parse_gbff ======================================================================
def parse_gbff(gbff_input):
    species, feature_list = _get_gbff_info_and_split(gbff_input)
    parsed_regions = [_parse_region(region) for region in feature_list]
    original_dict_list = [_make_dict(region) for region in parsed_regions]
    main_dict = _filter_dict(original_dict_list)

    df = pd.DataFrame.from_dict(main_dict, orient='index')
    key_order = ["protein_id", "gene", "EC_number", "product", "translation"]
    df = df[key_order]

    return species, df   

def _get_gbff_info_and_split(file_path):
    species = ""
    split_regions = []

    with open(file_path, 'r') as file:
        species_flag = False
        lines = file.readlines()
        feature_list = []
        split_regions = []

        for i, line in enumerate(lines):
            if "FEATURES" in line:
                contig_line_index = None
                for j in range(i + 1, len(lines)):
                    if 'CONTIG' in lines[j]:
                        contig_line_index = j
                        break
                    elif 'ORIGIN' in lines[j]:
                        contig_line_index = j
                        break

                if contig_line_index is not None:
                    feature_list.extend(lines[i + 1:contig_line_index])
                    split_regions.append(''.join(feature_list).strip())
            elif species_flag is False and "ORGANISM" in line:
                tmp = line.split()
                species = tmp[1] + " " + tmp[2]
                species_flag = True
            else:
                continue
        
        split_regions = _split_into_regions(feature_list)
        
    return species, split_regions

def _split_into_regions(lines):
    split_regions = []
    cur_region = ""

    for line in lines:
        line = line.strip(' ')
        if line.startswith("gene  "):
            if cur_region:
                split_regions.append(cur_region)
            cur_region = line
        else:
            cur_region += line

    if cur_region:
        split_regions.append(cur_region)

    return split_regions[1:]

def _parse_region(region):
    tmp_region = ""
    tmp_CDS = ""

    for line in region.split('\n'):
        line = line.lstrip()
        tmp_region += line + '\n'

        if line.startswith("CDS  "):
            tmp_CDS = line

    tmp_region_list = tmp_region.split('\n')
    i = 1
    while i < len(tmp_region_list):
        if '=' not in tmp_region_list[i]:
            tmp_region_list[i - 1] += ' ' + tmp_region_list[i]
            tmp_region_list.pop(i)
        else:
            i += 1
    
    tmp_region_list.pop(0)

    position = tmp_CDS.split(" ")[-1]
    tmp_region_list.append(f'/position="{position}"')
    
    return tmp_region_list

def _make_dict(parsed_regions):
    result_dict = {}

    for region in parsed_regions:
        key_start = region.find('/') + 1
        key_end = region.find('=')
        key = region[key_start:key_end].strip()

        value_start = key_end + 1
        value_end = len(region)
        value = region[value_start:value_end].strip().strip('"')

        if key == "translation":
          result_dict[key] = value.replace(" ", "")
        else:
          result_dict[key] = value

    return result_dict
    
def _filter_dict(original_dict_list):
    main_dict = {}

    for original_dict in original_dict_list:
        protein_id = original_dict.get('protein_id')
        if protein_id:
            main_dict[protein_id] = {}

            for key, include in gbff_key_setting.items():
                if include:
                    value = original_dict.get(key)
                    main_dict[protein_id][key] = value if value is not None else "-"
                if key == "translation" and '"' in value:
                    main_dict[protein_id][key] = value.split('"')[0]

    return main_dict

#===== parse_gff ======================================================================
def parse_gff(gff_input):
    data_dict = _read_and_process_file(gff_input)

    summarized_list = list(data_dict.values())
    df = pd.DataFrame(summarized_list)
    filtered_df = df[df['protein_id'] != '-']
    
    return filtered_df 

def _process_line(line, key_count, summarized_list):

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

def _read_and_process_file(gff_input):
    key_count = 0
    summarized_list = {}

    with open(gff_input, 'r') as file:
        for line in file:
            key_count, summarized_list = _process_line(line, key_count, summarized_list)

    return summarized_list

#=======================================================================================
def merge_gbff_gff(input_info, gbff_summary, gff_summary):
    merged_df = pd.merge(gbff_summary, gff_summary, on='protein_id')
    merged_df = merged_df[['protein_id', 'contig', 'strand', 'start', 'end', 'gene', 'EC_number', 'product', 'translation']]

    genome_summary = f"{input_info[1]}/genome_summary.tsv"

    with open(genome_summary, 'w') as summary_file:
        summary_file.write(f"# {input_info[0]}\n# {input_info[2]}\n")
        merged_df.to_csv(summary_file, sep='\t', index=False)

    return merged_df