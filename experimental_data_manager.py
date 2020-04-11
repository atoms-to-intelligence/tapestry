import pandas as pd
import numpy as np
import math

import config

def read_harvard_data_cts():
  t = 24
  df = pd.read_csv("harvard_test1.csv")
  #cts = df.values[:, 0]
  fl = df.values[:, 1:]

  assert t == fl.shape[1]
  cts = np.zeros(t) + config.cycle_time_cutoff
  for j in range(t):
    for i in range(fl.shape[0]):
      if fl[i, j] >= 1000:
        fl1 = fl[i-1, j]
        fl2 = fl[i, j]

        ct = i + math.log(1000/fl1) / math.log(fl2/fl1)
        #print(j, ct)
        cts[j] = ct
        break

  return cts


