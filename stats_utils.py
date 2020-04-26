# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
# Utility methods to parse stats dictionary

from cs_expts import CSExpts

import numpy as np

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
  stats_list = [compute_stats_for_batch(batch) for batch in batches]


  # Mean and confidence intevals for each stat computed of 120 batches of 1000
  # expts.
  #
  # means and intervals are dicts keyed by stat name such as precision, recall
  # interval gives a tuple (low, high), both inclusive.
  intervals = get_intervals(stats_list)


  # This is stats computed on the actual list of experiments
  original = compute_stats_for_batch(explist)
  return original, intervals


# Make many batches of size 1000 from 1000 expts by resampling
def make_many_batches(explist, n_batches):
  return np.random.choice(explist, size=(n_batches, len(explist)))


# Computes stats for one batch of experiments
def compute_stats_for_batch(batch):
  # We already have the CSExpts class - we'll just use that to compute stats
  agg_expt = CSExpts
  for expt in batch:
    agg_expt.add_stats(expt.tp, expt.fp, expt.fn, expt.uncon_negs,
        expt.determined, expt.overdetermined,
        expt.surep, expt.unsurep, expt.wrongly_undetected, expt.score)

  # This method computes more stats. It returns something but we don't care
  agg_expt.return_stats(len(batch))

  agg_stats = dict(agg_expt.__dict__)

  # remove unneeded attributes of CSExpts
  del agg_stats['single_expts']
  
  return agg_stats


# Compute confidence intervals for each stat using the list of stats
def get_intervals(stats_list):
  pass

if __name__ == '__main__':
  print(make_many_batches([1, 2, 3], 10))


