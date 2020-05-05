# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8

import sys
sys.path.append(".")

from matrix_utils import *

if __name__ == '__main__':
  d = {}
  extra_mlabels = []
  load_extra_mats(d, extra_mlabels)
  print(extra_mlabels)
  print([ (key, d[key].shape) for key in d])

