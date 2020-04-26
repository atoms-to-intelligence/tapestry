# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
# Utility methods to parse stats dictionary


# Use bootstrapping to compute confidence intervals
#
# Get a list of 1000 expts
#
# compute precision, recall, specificity etc for
def parse_stats_and_get_confidence_intervals(explist, k=3, n_batches=120):
  # 120 batches of 1000 experiments, gotten using resampling
  batches = make_many_batches(explist, n_batches)

  # 120 aggregate stats such as precision, recall etc
  # each individual stat is a dict containing all stats
  stats = [compute_stats_for_batch(batch) for batch in batches]

  # Mean and confidence intevals for each stat computed of 120 batches of 1000
  # expts.
  #
  # means and intervals are dicts keyed by stat name such as precision, recall
  # interval gives a tuple (low, high), both inclusive.
  means, intervals = get_means_and_intervals(stats)
  
  return means, intervals


# Make many batches of size 1000 from 1000 expts by resampling
def make_many_batches(explist, n_batches):
  pass


# Computes stats for one batch of experiments
def compute_stats_for_batch(explist):
  pass


