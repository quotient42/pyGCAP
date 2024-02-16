import os
import pandas as pd
from Bio.Blast import NCBIXML
import time
import sys

sys.path.append("../seqlib/")
from collect_data import collect_target_data
from search_repseq import search_repseq
from construct_lib import construct_lib

def create_seqlib(project_info):
	print("<< creating seqlib...")
	data_dir = project_info['data']
	seqlib_dir = project_info['seqlib']

	# collect_target_data(project_info)

	search_repseq(project_info)
	construct_lib(project_info)