# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
import os
import gzip
import pickle
import shutil

class PickleManager:
  def __init__(self, main_pickle, tmp_pickle, default):
    self.main_pickle = main_pickle
    self.tmp_pickle = tmp_pickle
    self.default = default


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

# Manage experimental stats using a directory based structure instead of
# having a single pickle for all stats
#
# Top-level stats dir: expt_stats/
# One dir for each matrix: e.g. optimized_M_45_330_kirkman/
# Inside it, one dir for each algo: e.g. COMP/ or SBL/
# Inside it, one dir for each 'd': e.g. 1/ 2/ 3/ etc
# Inside it, one pickle containing list of SingleExpt.__dict__, 1000 expts.
class DirectoryStatsManager:
  pass


expt_stats_pickle_manager = PickleManager(config.stats_pickle,
        config.stats_pickle_tmp, default={})
