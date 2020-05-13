# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8


# Each new algorithm added to the algos/ directory must be added to this file
# by importing like this
from . import zeros
from . import clustered_sbl

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
    "SBL_clustered" : clustered_sbl.clustered_sbl,
    "SBL_majority_10_0.3" : (lambda params : clustered_sbl.sbl_multi(params, 'majority', 10, 0.3)),
    "SBL_majority_10_0.5" : (lambda params : clustered_sbl.sbl_multi(params, 'majority', 10, 0.5)),
    "SBL_majority_10_0.7" : (lambda params : clustered_sbl.sbl_multi(params, 'majority', 10, 0.7)),
    "SBL_majority_100_0.3" : (lambda params : clustered_sbl.sbl_multi(params, 'majority', 100, 0.3)),
    "SBL_majority_100_0.5" : (lambda params : clustered_sbl.sbl_multi(params, 'majority', 100, 0.5)),
    "SBL_majority_100_0.7" : (lambda params : clustered_sbl.sbl_multi(params, 'majority', 100, 0.7)),
    "SBL_union_10" : (lambda params : clustered_sbl.sbl_multi(params, 'union', 10, 0.7)),
    "SBL_intersection_10" : (lambda params : clustered_sbl.sbl_multi(params, 'intersection', 10, 0.7)),
    "SBL_union_100" : (lambda params : clustered_sbl.sbl_multi(params, 'union', 100, 0.7)),
    "SBL_intersection_100" : (lambda params : clustered_sbl.sbl_multi(params, 'intersection', 100, 0.7)),
    # e.g.
    # "SBL" : sbl.sbl,
    # "SBL_cv" : sbl.sbl_cv, etc
    }

