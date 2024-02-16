'''
This file contains functionality to parse GenBank files and extract relevant information.

Functions:
- get_gbff_info_and_split: Extracts species information and splits the input GenBank file into regions based on the FEATURES section.
- split_into_regions: Splits the input lines into regions based on the 'gene' feature.
- parse_region: Parses a region into a list of key-value pairs.
- make_dict: Creates a dictionary from parsed regions.
- filter_dict: Filters the dictionary based on predefined key settings.
- parse_gbff_main: Parses the input GenBank file and returns species information and a DataFrame of parsed regions.

'''

import pandas as pd
import sys

from gbff_info import gbff_key_setting

def get_gbff_info_and_split(file_path):
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

                if contig_line_index is not None:
                    feature_list.extend(lines[i + 1:contig_line_index])
                    split_regions.append(''.join(feature_list).strip())
            elif species_flag is False and "ORGANISM" in line:
                tmp = line.split()
                species = tmp[1] + " " + tmp[2]
                species_flag = True
            else:
                continue
        
        split_regions = split_into_regions(feature_list)
        
    return species, split_regions

def split_into_regions(lines):
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

def parse_region(region):
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

def make_dict(parsed_regions):
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
    
def filter_dict(original_dict_list):
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

def parse_gbff_main(gbff_input):
    species, feature_list = get_gbff_info_and_split(gbff_input)
    parsed_regions = [parse_region(region) for region in feature_list]
    original_dict_list = [make_dict(region) for region in parsed_regions]
    main_dict = filter_dict(original_dict_list)

    df = pd.DataFrame.from_dict(main_dict, orient='index')
    key_order = ["protein_id", "gene", "EC_number", "product", "translation"]
    df = df[key_order]

    return species, df 
