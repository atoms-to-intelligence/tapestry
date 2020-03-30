from cs import *
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
    item = do_many_expts(n, d, t, num_expts, xs, M,
        cross_validation=False,
        add_noise=True)
    #print(item)
    stats.append(item)

  return stats

n = 60
d = 2
t = 24
#Ms = [np.random.binomial(1, 0.5, size=(t, n)) for i in range(10)]
#Ms = [optimized_M for i in range(10)]
Ms = []
Ms.append(optimized_M_3)
#Ms.append(optimized_M_4)
ss = []
for d in range(1,3):
  stats = find_best_matrix(n, d, t, Ms, num_expts=1000)
  ss.extend(stats)
  #for stat, M in zip(stats, Ms):
  #  if stat['precision'] >= 0.80 and stat['recall'] >= 0.90:
  #    print(stat)

for stat in ss:
  s = json.dumps(stat)
  print(s)

#print(stat)
