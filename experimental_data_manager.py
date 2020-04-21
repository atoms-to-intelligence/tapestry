import pandas as pd
import numpy as np
import math

import os

import config

from matrices import MDict

def read_harvard_data_cts():
  t = 24
  harvard_data_file = os.path.join(config.root_dir, "harvard_test1.csv")
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
  pos_idx, cts = read_harvard_data_cts()
  print('Harvard data (24x60)')
  print('Positive indicies:', pos_idx)
  print('Cycle times:', "[", ", ".join(["%.1f" % item for item in cts]), "]")

  print('\nFake Data (16x40)\n')
  for i in range(10):
    pos_idx, cts = get_random_fake_test_data('16x40', 'optimized_M_16_40_ncbs')
    print('Positive indicies:', pos_idx)
    print('Cycle times:', "[", ", ".join(["%.1f" % item for item in cts]), "]")

