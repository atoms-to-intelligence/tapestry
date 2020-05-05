# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8


# Each new algorithm added to the algos/ directory must be added to this file
# by importing like this
from . import zeros

# Dictionary containing a mapping of algorithm name to corresponding function. 
#
# ANY NEW ALGORITHM MUST BE ADDED TO THIS DICTIONARY
#
# Below, "ZEROS" is the name of the algorithm. The algorithm name must not conflict with
# any of the names used in cs.py (function decode_lasso). If you have an
# algorithm with some name "FOO", then corresponding name "combined_COMP_FOO"
# is also valid. It pre-processes the matrix using COMP and deletes rows and
# columns from A and y which are definitely 0 according to COMP, and sends
# this reduced matrix. In your own algorithm, do not perform such
# pre-processing, rather, use algo "combined_COMP_FOO" to run "FOO" with this
# pre-processed data

# See zeros.py for the description of function signature and return value
algo_dict = {
    "ZEROS" : zeros.zeros,
    "ZEROS_cv" : zeros.zeros_cv,
    # e.g.
    # "SBL" : sbl.sbl,
    # "SBL_cv" : sbl.sbl_cv, etc
    }

