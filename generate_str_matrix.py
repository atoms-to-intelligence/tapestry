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

def adj_list_as_csv_block_format(M, name, block_size=15, max_wells=3):
  t = M.shape[0]
  n = M.shape[1]
  adj_list = {}
  for j in range(n):
    sample = j + 1
    tests, = np.nonzero(M[:,j])
    wells = [idx_to_well(item, t) for item in tests]
    assert len(wells) <= max_wells
    if len(wells) < max_wells:
      wells = wells + ['']*(max_wells - len(wells))

    # Add gap of two rows
    wells = wells + ['']*2
    if j < block_size:
      adj_list[str(sample)] = wells
    else:
      lst = [sample] + wells
      col = j % block_size + 1
      adj_list[str(col)].extend(lst)

  # Now extend the last few blocks
  max_len = len(adj_list["1"])
  for col in [str(j) for j in range(1, block_size + 1)]:
    l = len(adj_list[col])
    assert l <= max_len
    if l < max_len:
      adj_list[col].extend(['']*(max_len - l))

  df = pd.DataFrame(data=adj_list)
  df.to_csv(name, sep=',', index=False, quoting=1)

def create_all_csvs():
  #mlabel = "optimized_M_16_40_ncbs"
  for msize in MSizeToLabelDict:
    mlabel, _1, _2 = MSizeToLabelDict[msize]
    create_csv(mlabel)


def create_csv(mlabel, msize, block=False):
  csv_dir = os.path.join(config.root_dir, "csv")
  mcodename = mat_codenames[mlabel]
  if not block:
    mfilename = mlabel + ".csv" 
  else:
    mfilename = mlabel + "_block.csv" 
  mfilename = os.path.join(csv_dir, mfilename)
  if not block:
    adj_list_as_csv(MLabelToMatrixDict[mlabel], mfilename)
  else:
    adj_list_as_csv_block_format(MLabelToMatrixDict[mlabel], mfilename)
  with open(mfilename, "a") as f:
    f.write(f'{msize},{mcodename}\n')

create_csv('optimized_M_3', "24x60", block=True)

#csv_dir = os.path.join(config.root_dir, "csv")
#mlabel = "optimized_M_3"
#mcodename = mat_codenames[mlabel]
#mfilename = mlabel + "_block.csv" 
#mfilename = os.path.join(csv_dir, mfilename)
#adj_list_as_csv_block_format(MLabelToMatrixDict[mlabel], mfilename)

