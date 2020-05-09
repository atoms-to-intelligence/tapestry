# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8

if __name__=='__main__':
  import sys
  # This is a hack needed so that you can directly run this file as 
  # python inbuilt_algos/nnompcv.py.
  # 
  # This is needed because we want to import "matrices" from one level up, but python
  # does not know about one level up. Please only use this hack in test files
  # and not in files which are meant to be imported from other code.
  sys.path.append(".")

from utils.matrix_utils import *

d = {}
extra_mlabels = []
load_extra_mats(d, extra_mlabels)
print(extra_mlabels)
print([ (key, d[key].shape) for key in d])

