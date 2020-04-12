# Global configuration parameters which are expected to change less frequently

# Trade off precision for recall in CS algorithms
prefer_recall = False

########     Following configs are for the app backend only     ########

root_dir = '.'

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

