# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
from core import config

import os
import gzip
import pickle
import shutil

# Manages a single stats pickle. Agnostic to the data structure being stored
# in the pickle. Can be used to store entire stats in a dictionary.
class PickleManager:
  def __init__(self, main_pickle, tmp_pickle, default):
    self.main_pickle = main_pickle
    self.tmp_pickle = tmp_pickle
    self.default = default


  # helper to load pickle
  def load_pickle(self, name):
    with gzip.open(name, "rb") as f:
      stats = pickle.load(f)
    return stats


  # helper to check and load stats dict
  def load_stats(self):
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


# Manage experimental stats using a directory based structure instead of
# having a single pickle for all stats
#
# Top-level stats dir: expt_stats/
# One dir for each matrix: e.g. optimized_M_45_330_kirkman/
# Inside it, one dir for each algo: e.g. COMP/ or SBL/
# Inside it, one dir for each 'd': e.g. 1/ 2/ 3/ etc
# Inside it, one pickle containing list of SingleExpt.__dict__, 1000 expts.
class DirectoryStatsManager:
  def __init__(self, stats_dir):
    self.stats_dir = stats_dir

  # Save stats to pickle in directory
  # stats_dir/mlabel/algo/d/expt_stats.p.gz
  def save(self, mlabel, algo, d, explist):
    dirs = os.path.join(self.stats_dir, mlabel, algo, str(d))
    os.makedirs(dirs, exist_ok=True)
    main_pickle = os.path.join(dirs, config.stats_pickle_name)
    tmp_pickle = os.path.join(dirs, config.stats_pickle_tmp_name)
    pm = PickleManager(main_pickle, tmp_pickle, [])
    pm.carefully_save_stats(explist)


  # load stats list from pickle in directory
  # stats_dir/mlabel/algo/d/expt_stats.p.gz
  def load(self, mlabel, algo, d):
    self.do_directory_checks(mlabel, algo, d)
    dirs = os.path.join(self.stats_dir, mlabel, algo, str(d))
    main_pickle = os.path.join(dirs, config.stats_pickle_name)
    tmp_pickle = os.path.join(dirs, config.stats_pickle_tmp_name)
    pm = PickleManager(main_pickle, tmp_pickle, [])
    return pm.load_stats()

  def do_directory_checks(self, mlabel, algo, d):
    if not os.path.exists(self.stats_dir):
      raise ValueError(f"Directory {self.stats_dir} does not exist")

    mlabel_dir = os.path.join(self.stats_dir, mlabel)
    if not os.path.exists(mlabel_dir):
      raise ValueError(f"Directory {mlabel_dir} does not exist")

    algo_dir = os.path.join(mlabel_dir, algo)
    if not os.path.exists(algo_dir):
      raise ValueError(f"Directory {algo_dir} does not exist")

    d_dir = os.path.join(algo_dir, str(d))
    if not os.path.exists(d_dir):
      raise ValueError(f"Directory {d_dir} does not exist")

expt_stats_pickle_manager = PickleManager(config.stats_pickle,
        config.stats_pickle_tmp, default={})

stats_manager = DirectoryStatsManager(config.stats_dir)

