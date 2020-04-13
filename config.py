# Global configuration parameters which are expected to change less frequently

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
root_dir = '.'

########     Following configs are for the app backend only     ########


# Decoding algorithm to be used with the app backend.
app_algo = 'COMP'


# Cycle time cutoff. Samples with greater than or equal to this value are
# considered negative
cycle_time_cutoff = 35


# This is the probability of replication of RNA. Number of molecules after t
# cycles is y_0*(1+p)**t
p = 0.95


########     Following config is for synthetic experiments only     ########


# Flip +ve bits of y with this prob
bit_flip_prob = 0.


# Exponential Gaussian or Variable Gaussian
#noise_model = 'exponential_gaussian'
noise_model = 'variable_gaussian'

