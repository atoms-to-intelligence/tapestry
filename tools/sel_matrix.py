import sys
sys.path.append(".")

from core.cs_expts import *
import pickle

import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.3f')

# Out of a list of matrices, choose the best one
#
# Generate 1000 x's to do the expt. The same x will be used to evaluate all
# the matrices
# list of stats will be returned for each matrix
def find_best_matrix(n, d, t, Ms, num_expts=1000):
  xs = [create_infection_array_with_num_cases(n, d) for i in range(num_expts)]
  stats = []
  for M in Ms:
    _, item = do_many_expts(n, d, t, num_expts, xs, M,
        cross_validation=False,
        add_noise=True,
        algo='NNOMP_random_cv',
        mr=11)
    stats.append(item)

  return stats

n = 40
d = 2
t = 16
#Ms = [np.random.binomial(1, 0.5, size=(t, n)) for i in range(10)]
#Ms = [optimized_M for i in range(10)]
Ms = []
Ms.append(optimized_M_1)
Ms.append(optimized_M_2)
ss = []
for d in range(1, 11):
  stats = find_best_matrix(n, d, t, Ms, num_expts=1)
  ss.extend(stats)
  #for stat, M in zip(stats, Ms):
  #  if stat['precision'] >= 0.80 and stat['recall'] >= 0.90:
  #    print(stat)

for stat in ss:
  s = json.dumps(stat)
  print(s)

#print(stat)
