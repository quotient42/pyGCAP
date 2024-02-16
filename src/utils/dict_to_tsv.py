import pandas as pd
from target_info import target_list

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.max_colwidth', None)

df = pd.DataFrame(target_list)
df.to_csv(f"../../data/target_dictionary.tsv", sep='\t', index=False)
