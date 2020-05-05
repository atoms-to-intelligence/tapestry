# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
import pandas as pd
import numpy as np
import math

import os

import config

from matrices import MDict

# Standard format
# 1 entry per line
def read_standard_cts(data_file):
  #data_dir = config.data_dir
  #data_file = os.path.join(data_dir, fname)
  with open(data_file, "r") as f:
    cts = [ float(item.strip()) for item in f.readlines() ]
  cts = np.array(cts)
  cts[cts == 0] = config.cycle_time_cutoff + 1
  return cts

def read_harvard_data_cts():
  t = 24
  harvard_data_file = os.path.join(config.data_dir, 'harvard', "harvard_test1.csv")
  df = pd.read_csv(harvard_data_file)
  #cts = df.values[:, 0]
  fl = df.values[:, 1:]

  assert t == fl.shape[1]
  # config.cycle_time_cutoff may keep changing. Let's just make this 100.
  #cts = np.zeros(t) + config.cycle_time_cutoff
  cts = np.zeros(t) + 100
  for j in range(t):
    for i in range(fl.shape[0]):
      if fl[i, j] >= 1000:
        fl1 = fl[i-1, j]
        fl2 = fl[i, j]

        ct = i + math.log(1000/fl1) / math.log(fl2/fl1)
        #print(j, ct)
        cts[j] = ct
        break

  # Ground truth positives. Indices start from 1.
  pos_idx = [10, 28]
  return pos_idx, cts

def get_israel_data_cts():
  data_file = os.path.join(config.root_dir, "israel_2_positives_cts.txt")
  with open(data_file, "r") as f:
    cts = [float(item.strip()) for item in f.readlines()]
    print('positive ys:', np.argwhere(cts) + 1)
  cts = np.array(cts)
  cts[cts == 0] = config.cycle_time_cutoff + 1
  return cts

def parse_israel_matrix():
  A = np.zeros((48,384), dtype=np.int32)
  data_file = os.path.join(config.mat_dir, "israel_48_384_matrix.txt")
  with open(data_file) as f:
    for line in f.readlines()[1:]:
      words = line.strip().split()
      row = int(words[0]) - 1
      col = int(words[-1]) - 1
      #print(row, col)
      A[row, col] = 1
  outfile = os.path.join(config.mat_dir, "israel_48_384.txt")
  np.savetxt(outfile, A, fmt="%d")
  return A

# Returns randomly generated bool_x and cycle times.
# Needs matrix size and  matrix label
#
# Process of generation is as follows:
#
# 1. First choose sparsity 'd' in range 1...5.
# 2. Generate bool_x array of size 'n'
#
def get_random_fake_test_data(mat_size, mat_label):
  d = np.random.randint(0, 10)
  t = int(mat_size.split('x')[0])
  n = int(mat_size.split('x')[1])
  pos_idx = np.random.choice(list(range(n)), size=d)
  bool_x = np.zeros(n)
  bool_x[pos_idx] = 1
  x = np.random.rand(n) * bool_x
  M = MDict[mat_label]
  assert t == M.shape[0]
  assert n == M.shape[1]
  y = np.matmul(M, x)
  bool_y = (y > 0).astype(np.int32)
  ct1 = np.random.randint(1, config.cycle_time_cutoff, t)
  ct2 = np.random.randint(config.cycle_time_cutoff, 50, t)
  cts = ct1 * bool_y + ct2 * (1 - bool_y)
  return [i + 1 for i in pos_idx], cts
  
if __name__ == '__main__':
  cts = get_israel_data_cts()
  print(cts, len(cts))
  A = parse_israel_matrix()
  #print(A)
  print(np.sum(A))
  print(np.sum(A, axis=1))
  print(np.sum(A, axis=0))
  from app_utils import get_test_results
  from get_test_results import get_result_string_from_lists

  n = 384
  config.app_algo = 'SBL'
  sure_list, unsure_list, neg_list, x = get_test_results(A, cts)
  #print('sure_list:', sure_list)
  #print('unsure_list:', unsure_list)
  #print('neg_list:', neg_list)
  #print('x:', x)
  result_string = get_result_string_from_lists(sure_list, unsure_list,
      neg_list, x, n)
  print(result_string)
  #pos_idx, cts = read_harvard_data_cts()
  #print('Harvard data (24x60)')
  #print('Positive indicies:', pos_idx)
  #print('Cycle times:', "[", ", ".join(["%.1f" % item for item in cts]), "]")

  #print('\nFake Data (16x40)\n')
  #for i in range(10):
  #  pos_idx, cts = get_random_fake_test_data('16x40', 'optimized_M_16_40_ncbs')
  #  print('Positive indicies:', pos_idx)
  #  print('Cycle times:', "[", ", ".join(["%.1f" % item for item in cts]), "]")

