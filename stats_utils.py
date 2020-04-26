# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
# Utility methods to parse stats dictionary

import config
from cs_expts import CSExpts

import numpy as np
import os
import gzip
import pickle
import shutil
import json

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
  intervals = get_intervals(stats_list, k, keys=["precision"])


  # This is stats computed on the actual list of experiments
  original = compute_stats_for_batch(explist)
  return original, intervals


# Make many batches of size 1000 from 1000 expts by resampling
def make_many_batches(explist, n_batches):
  return np.random.choice(explist, size=(n_batches, len(explist)))


# Computes stats for one batch of experiments
def compute_stats_for_batch(batch):
  expt = batch[0]
  # We already have the CSExpts class - we'll just use that to compute stats
  agg_expt = CSExpts('dummy_name', expt["n"], expt["d"], expt["t"], expt["mr"],
    expt["algo"], expt["mlabel"])
  for expt in batch:
    agg_expt.add_stats(expt["tp"], expt["fp"], expt["fn"], expt["uncon_negs"],
        expt["determined"], expt["overdetermined"],
        expt["surep"], expt["unsurep"], expt["wrongly_undetected"],
        expt["score"])

  # This method computes more stats. It returns something but we don't care
  agg_expt.return_stats(len(batch))

  agg_stats = dict(agg_expt.__dict__)

  # remove unneeded attributes of CSExpts
  del agg_stats['single_expts']
  
  return agg_stats


# Compute confidence intervals for each stat using the list of stats
def get_intervals(agg_stats_list, k, keys=None):
  assert agg_stats_list
  if not keys:
    # Get which stats we are going to compute the confidence intervals for
    keys = agg_stats_list[0].keys()
  
  intervals = {}
  for key in keys:
    stats = sorted([agg_stats[key] for agg_stats in agg_stats_list])
    # Drop the first k and the last k to get interval
    low = stats[k]
    high = stats[-k-1]
    intervals[key] = (low, high)

  return intervals

class PickleManager:
  def __init__(self, main_pickle, tmp_pickle):
    self.main_pickle = main_pickle
    self.tmp_pickle = tmp_pickle
    self.default = {}


  # helper to load pickle
  def load_pickle(self, name):
    with gzip.open(self.main_pickle, "rb") as f:
      stats = pickle.load(f)
    return stats


  # helper to check and load stats dict
  def get_stats_dict(self):
    if os.path.exists(self.tmp_pickle):
      raise ValueError("Temporary pickle file found. Please check if this "
          "contains valid data")

    if os.path.exists(self.main_pickle):
      stats = self.load_pickle(self.main_pickle)
    else:
      stats = self.default
    return stats


  # Saves to tmp file first, then copies the tmp file onto the old file. Then
  # deletes the tmp file
  def carefully_save_stats(self, stats):
    with gzip.open(self.tmp_pickle, "wb") as f:
      pickle.dump(stats, f)
    shutil.copy2(self.tmp_pickle, self.main_pickle)
    os.remove(self.tmp_pickle)


if __name__ == '__main__':
  #print(make_many_batches([1, 2, 3], 10))
  pm = PickleManager(config.stats_pickle, config.stats_pickle_tmp)
  stats = pm.get_stats_dict()
  explist = stats["optimized_M_1"]["COMP"][2]
  original, intervals = parse_stats_and_get_confidence_intervals(explist, k=3, n_batches=120)

  s = json.dumps(original["precision"])
  print(s)
  s = json.dumps(intervals)
  print(s)

