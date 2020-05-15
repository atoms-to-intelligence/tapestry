# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
# Utility methods to parse stats dictionary

import sys
sys.path.append(".")

import numpy as np
import json
import re

from core import config
from core.cs_expts import CSExpts
from utils.pickle_manager import stats_manager, expt_stats_pickle_manager
from core.matrices import kirkman_mlabels, MDict


# Use bootstrapping to compute confidence intervals
#
# Get a list of 1000 expts
#
# compute precision, recall, specificity etc for
def parse_stats_and_get_confidence_intervals(explist, k=3, n_batches=120,
    keys=None):
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
  #intervals = get_intervals(stats_list, k, keys=["precision"])
  intervals = get_intervals(stats_list, k, keys)


  # This is stats computed on the actual list of experiments
  original = compute_stats_for_batch(explist)
  convert_dict_to_simple_type(original)
  
  combined = combined_original_and_intervals(original, intervals)
  return combined


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
        expt["score"], expt['rmse'])

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
    if key in ["n", "t", "d", "mr", "algo", "name", "mlabel"]:
      continue
    #print("Processing key: ", key)
    stats = sorted([agg_stats[key] for agg_stats in agg_stats_list])
    # Drop the first k and the last k to get interval
    low = stats[k]
    high = stats[-k-1]
    low = convert_scalar_to_simple_type(low)
    high = convert_scalar_to_simple_type(high)
    intervals[key] = (low, high)

  return intervals

def combined_original_and_intervals(original, intervals):
  combined = {}
  for key in intervals:
    combined[key] = (original[key], intervals[key][0], intervals[key][1])
  return combined

# Converts each item in dict from numpy to python types. This makes it
# serializable with json
def convert_dict_to_simple_type(d):
  for key in d:
    try:
      d[key] = d[key].item()
      #print(f"Converted {key} to python")
    except:
      #print(f"Could not convert {key} to python")
      pass

def convert_scalar_to_simple_type(val):
  try:
    return val.item()
  except:
    return val

def get_stats_for_deployed_matrices():
  #pm = PickleManager(config.stats_pickle, config.stats_pickle_tmp)
  #stats = pm.get_stats_dict()
  #explist = stats["optimized_M_1"]["COMP"][2]

  from get_test_results import MSizeToLabelDict
  tups = MSizeToLabelDict.values()
  sizes = MSizeToLabelDict.keys()
  labels = [tup[0] for tup in tups]
  ds = [tup[1] for tup in tups]
  algos = ['COMP', 'SBL', 'combined_COMP_NNOMP_random_cv', 'combined_COMP_l1ls_cv']

  get_stats_for_these_matrices(sizes, labels, ds, algos)

def get_stats_for_kirkman_matrices():
  labels = kirkman_mlabels
  sizes = [f"{M.shape[0]}x{M.shape[1]}" for M in [MDict[label] for label in
    labels]]
  ds = [ 3 for label in labels]
  algos = ['COMP', 'SBL']
  get_stats_for_these_matrices(sizes, labels, ds, algos)

def get_stats_for_these_matrices(sizes, labels, ds, algos):
  for size, label, d_range in zip(sizes, labels, ds):
    for algo in algos:
      for d in d_range:
        print(f"Matrix: {size}, Algo: {algo}, d = {d}")
        explist = stats_manager.load(label, algo, d)[0:100]
        #print(f'size: {size}, batch len: ', len(explist))
        combined = parse_stats_and_get_confidence_intervals(explist, k=3,
            n_batches=120, keys=['precision', 'recall', 'specificity', 'avg_surep', 'avg_unsurep',
            'avg_fp', 'avg_fn'])
        s = json.dumps(combined, indent=2)
        pat = re.compile(r"\d+\.\d{4,}")
        def mround(match):
          return "{:.3f}".format(float(match.group()))

        sys.stdout.write(re.sub(pat, mround, s) + '\n')
        #print(s)

# Migrate from using single PickleManager to a directory structure
def migrate_stats_from_dict_to_directory():
  stats = expt_stats_pickle_manager.load_stats()
  for mlabel in stats:
    for algo in stats[mlabel]:
      for d in stats[mlabel][algo]:
        stats_manager.save(mlabel, algo, d, stats[mlabel][algo][d])

def get_stats_2_pct_prevalance_matrices():
  labels = [
      "optimized_M_27_117_kirkman",
      "optimized_M_48_320_kirkman",
      #"optimized_M_60_500_kirkman",
      "optimized_M_75_800_kirkman[:, :500]",
      "optimized_M_192_5120_social_golfer[:, :1024]",
      ]
  sizes = [
      "27x117",
      "48x320",
      '75x500', 
      '192x1024',
      ]
  ds = [ 5, 11, 17, 30 ]
  algos = ['combined_COMP_SBL']
  get_stats_for_these_matrices(sizes, labels, ds, algos)

def get_stats_clinical_trials_TMH():
  labels = [
      'optimized_M_36_180_kirkman[:, :72]',
      'optimized_M_45_105_kirkman[:, :105]',
      ]
  sizes = [
      '36x72',
      '45x105',
      ]
  ds = [ [7, 9, 12], [10, 12, 15] ]
  algos = [ 'precise_SBL_combined_COMP_SBL' ]
  get_stats_for_these_matrices(sizes, labels, ds, algos)


if __name__ == '__main__':
  get_stats_clinical_trials_TMH()
  #get_stats_2_pct_prevalance_matrices()
  #migrate_stats_from_dict_to_directory()
  #get_stats_for_deployed_matrices()
  #get_stats_for_kirkman_matrices()
  #print(make_many_batches([1, 2, 3], 10))
  #pm = PickleManager(config.stats_pickle, config.stats_pickle_tmp)
  #stats = pm.get_stats_dict()
  #explist = stats["optimized_M_1"]["COMP"][2]
  #combined = parse_stats_and_get_confidence_intervals(explist, k=3,
  #    n_batches=120, keys=['precision', 'recall', 'specificity'])

  #s = json.dumps(combined, indent=2)
  #print(s)
else:
  raise ValueError("This is not a library file. Can't import this. Use standalone.")

