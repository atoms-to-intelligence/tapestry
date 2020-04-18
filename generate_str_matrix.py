# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
#from matrices import *
from get_test_results import mat_codenames, MSizeToLabelDict, MLabelToMatrixDict
import itertools as it
import config

import numpy as np
import pandas as pd
import os

C = [f'{t[0]}{t[1]}' for t in it.product('ABCDEFGH',range(1,13))]
#print(C, len(C))

# For each matrix, print out the well number
def idx_to_well(idx, tests):
  if tests < 48:
    return C[2*idx]
  else:
    return C[idx]


def print_well_nums(tests):
  for idx in range(tests):
    print(idx, idx_to_well(idx, tests))

#print_well_nums(48)
#print_well_nums(63)

def adj_list_as_csv(M, name):
  t = M.shape[0]
  n = M.shape[1]
  adj_list = {}
  for j in range(n):
    sample = j + 1
    tests, = np.nonzero(M[:,j])
    wells = [idx_to_well(item, t) for item in tests]
    assert len(wells) < 10
    wells = wells + ['']*(10 - len(wells))
    adj_list[str(sample)] = wells
  df = pd.DataFrame(data=adj_list)
  df.to_csv(name, sep=',', index=False)

csv_dir = os.path.join(config.root_dir, "csv")
#mlabel = "optimized_M_16_40_ncbs"
for msize in MSizeToLabelDict:
  mlabel, _1, _2 = MSizeToLabelDict[msize]
  mcodename = mat_codenames[mlabel]
  mfilename = mlabel + ".csv" 
  mfilename = os.path.join(csv_dir, mfilename)
  adj_list_as_csv(MLabelToMatrixDict[mlabel], mfilename)
  with open(mfilename, "a") as f:
    f.write(f'{msize},{mcodename}\n')


