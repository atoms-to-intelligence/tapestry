# Global configuration parameters which are expected to change less frequently

# Needed for paths
import os

# Trade off precision for recall in CS algorithms
prefer_recall = False


# This is the directory which contains the code. This is used to locate matrix
# data and experimental test data in matrices.py and
# experimental_data_manager.py. 
#
# The experiment code is added as a git submodule to the app backend
# repository under folder "compute". The app backend code changes
# config.root_dir to "compute" before importing any other file from the
# experiments code. This needs to be "." in order to be able to execute our
# code from the current working dir from command-line.
#
# Could have gone for os.environ in retrospect but this works as well.
#root_dir = '.'


# This is where printable pdfs for each matrix are kept
#mat_pdf_dir = os.path.join(root_dir, 'mat_pdfs')

# set_root_dir is now called for root dir changes. Either call this from your
# code or do the changes in this function yourself
def set_root_dir(new_root_dir):
  global root_dir, mat_pdf_dir, mat_dir, kirkman_dir, sts_dir, extra_mat_dir,\
          data_dir
  root_dir = new_root_dir
  mat_pdf_dir = os.path.join(root_dir, 'mat_pdfs')
  mat_dir = os.path.join(root_dir, 'mats')
  kirkman_dir = os.path.join(mat_dir, "kirkman")
  sts_dir = os.path.join(mat_dir, "sts")
  extra_mat_dir = os.path.join(mat_dir, 'extra')
  data_dir = os.path.join(root_dir, 'data')

set_root_dir('.')


########     Following configs are for the app backend only     ########


# Decoding algorithm to be used with the app backend.
app_algo = 'COMP'


# Cycle time cutoff. Samples with greater than or equal to this value are
# considered negative
cycle_time_cutoff = 40


# This is the probability of replication of RNA. Number of molecules after t
# cycles is y_0*(1+p)**t
p = 0.95


########     Following config is for synthetic experiments only     ########


# Flip +ve bits of y with this prob
bit_flip_prob = 0.


# Exponential Gaussian or Variable Gaussian
noise_model = 'exponential_gaussian'
#noise_model = 'variable_gaussian'

# Pickle files containing stats. Stats are first written to tmp and then
# finally copied
stats_pickle_name = "expt_stats.p.gz"
stats_pickle_tmp_name = "expt_stats_temp.p.gz"
stats_pickle = os.path.join(root_dir, stats_pickle_name) 
stats_pickle_tmp = os.path.join(root_dir, stats_pickle_tmp_name)

# Stats directory
#stats_dir_name = "expt_stats"
stats_dir = os.path.join(root_dir, "expt_stats")

