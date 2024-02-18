import time
import sys

sys.path.append("../seqlib/")
from make_target_fasta import make_target_fasta
from collect_data import collect_target_data
from blast_data import blast_data
from search_repseq import search_repseq
from construct_lib import construct_lib

def create_seqlib(project_info):
	start_time = time.time()

	make_target_fasta(project_info)
	collect_target_data(project_info)
	blast_data(project_info)
	search_repseq(project_info)
	construct_lib(project_info)

	end_time = time.time()
	total = end_time - start_time
	print(f"   seqlib created (elapsed time: {round(total / 60, 3)} min)")
	print("----------------------------------------------------------")
