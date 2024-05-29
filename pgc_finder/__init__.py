__version__ = '0.0.27'

from ._accessory import *
from ._blast import *
from ._cluster import *
from ._count import *
from ._parse import *
from ._profile import *
from ._search import *
from ._seqlib import *
from ._target import *
from ._utils import *
from ._visualize import *

import pkg_resources

DATA_PATH = 'data/target_sample.tsv'

def get_data_path():
    return pkg_resources.resource_filename(__name__, DATA_PATH)