import os
import pandas as pd
import sys
import time

sys.path.append('../info/')
from target_info import target_list

def count_target(project_info):
    print("<< counting target by species and genus level...")
    input_dir = project_info['input']
    output_dir = project_info['output']

    all_directories = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    all_directories = [d for d in all_directories if d != 'output' and not d.startswith('output/')]

    seq_lib = pd.read_csv(f"{project_info['seqlib']}/seqlib.tsv", sep='\t')
    seq_lib['Prediction'] = seq_lib['Prediction'].astype(str)
    seq_lib['RepSeq'] = seq_lib['RepSeq'].astype(str)
    seq_lib['Total'] = seq_lib['Total'].astype(int)
    seq_lib['Seq_List'] = seq_lib['Seq_List'].astype(str)

    seq_lib['no'] = seq_lib.groupby('Prediction').cumcount() + 1
    seq_lib['Prediction'] = seq_lib['Prediction'] + '-' + seq_lib['no'].astype(str)
    seq_lib.drop(['no', 'Total'], axis=1, inplace=True)

    genus_rep = seq_lib.copy()

    start_time = time.time()

    target_names = []
    for target in target_list:
        target_names.append(target['protein'])

    for genus in all_directories:
        for accession in os.listdir(os.path.join(input_dir, genus)):
            if not "GCF" in accession:
                continue
            file_path = f"{input_dir}/{genus}/{accession}/target_seq.tsv"
        
            data_df = pd.read_csv(file_path, sep='\t', comment='#')
            seq_lib = pd.concat([seq_lib, pd.DataFrame({accession: [0] * len(seq_lib)})], axis=1)
            seq_lib.loc[seq_lib['RepSeq'].isin(data_df['repseq']), accession] = 1

            genus_rep = pd.concat([genus_rep, pd.DataFrame({genus: [0] * len(genus_rep)})], axis=1)
            genus_rep.loc[genus_rep['RepSeq'].isin(data_df['repseq']), genus] = 1
    				
        seq_lib = pd.concat([seq_lib], axis=1)
        genus_rep = pd.concat([genus_rep], axis=1)

    seq_lib.drop(['RepSeq', 'Seq_List'], axis=1, inplace=True)

    seq_lib = seq_lib.transpose()
    seq_lib.rename(columns=seq_lib.iloc[0], inplace=True)
    seq_lib = seq_lib.reset_index()

    seq_lib = seq_lib[1:]
    seq_lib.to_csv(f"{output_dir}/tsv/repseq_species.tsv", sep='\t', index=False)

    genus_rep.drop(['RepSeq', 'Seq_List'], axis=1, inplace=True)
    genus_rep = genus_rep.transpose()
    genus_rep.rename(columns=genus_rep.iloc[0], inplace=True)
    genus_rep = genus_rep.reset_index()
    genus_rep = genus_rep[1:]
    genus_rep = genus_rep.groupby('index').mean().reset_index()
    genus_rep.to_csv(f"{output_dir}/tsv/repseq_genus.tsv", sep='\t', index=False)    

    end_time = time.time()
    total = end_time - start_time
    print(f"   └── counting done (elapsed time: {round(total / 60, 3)} min)")
    print("----------------------------------------------------------")

