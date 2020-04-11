import numpy as np

import config
from cs import CS

# Get results given a matrix M and cycle times cts
def get_test_results(M, cts):
  y = _get_y_from_cts(cts)
  sure_list, unsure_list, neg_list, x = _get_infected_lists(M, y)
  return sure_list, unsure_list, neg_list, x


# Get y (viral loads in each test) from cycle times
#
#
def _get_y_from_cts(cts):
  # First get the min cycle time. This corresponds to the test with the
  # maximum viral load. This sample will have the least amount of variance
  # in cycle time, so good to choose this as the baseline
  ctmin = np.min(cts)      

  # Now filter out the positive tests
  bool_y = (cts < config.cycle_time_cutoff).astype(np.int32)
  cts = cts * bool_y

  y = (1 + config.p) ** (ctmin - cts)
  y = y * bool_y

  print('cts:', cts)
  print('y:', y)
  print('bool_y:', bool_y)

  return y


# This function calls the decoding algorithm in class CS
# The algorithm chosen can be configured using config.app_algo
def _get_infected_lists(M, y):
  assert M is not None

  n = M.shape[1]
  t = M.shape[0]

  assert t == len(y)

  # unused params
  arr = np.zeros(n)
  mr = None
  d = 1
  s = 0.5
  l = 0.1


  bool_y = (y > 0).astype(np.int32)
  cs = CS(n, t, s, d, l, arr, M, mr)

  if config.app_algo == 'COMP':
    infected, infected_dd, score, tp, fp, fn, surep, unsurep,\
        num_infected_in_test = \
        cs.decode_comp_new(bool_y, compute_stats=False)
    x = np.zeros(n)
  else:
    x, infected, infected_dd, prob1, prob0, score, tp, fp, fn, uncon_negs, determined,\
        overdetermined, surep, unsurep, wrongly_undetected,\
        num_infected_in_test = cs.decode_lasso(y, config.app_algo, prefer_recall=False,
            compute_stats=False)

  _detect_discrepancies_in_test(t, bool_y, num_infected_in_test)
  sure_list, unsure_list, neg_list = _get_lists_from_infected(infected,
      infected_dd, n)
  
  return sure_list, unsure_list, neg_list, x

# Detects discrepancies which can happen due to experimental error or using
# wrong matrix, or poor algorithm performance.
def _detect_discrepancies_in_test(t, bool_y, num_infected_in_test):
  for test in range(t):
    if bool_y[test] > 0 and num_infected_in_test[test] == 0:
      print('y[%d] is infected but no infected people found' % test)
    if bool_y[test] == 0 and num_infected_in_test[test] > 0:
      print('y[%d] is not infected but infected people found' % test)

# Get list of sure, unsure and negative people
def _get_lists_from_infected(infected, infected_dd, n):
  sure_list = []
  unsure_list = []
  neg_list = []
  for i in range(n):
    if infected_dd[i] == 1:
      sure_list.append(i+1)
    elif infected[i] == 1:
      unsure_list.append(i+1)
    else:
      neg_list.append(i+1)

  return sure_list, unsure_list, neg_list

if __name__ == '__main__':
  import pandas as pd
  import math
  from matrices import optimized_M_3

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

  # Test using Harvard and ncbs data
  def run_tests():
    pass

  def harvard_test():
    cts = read_harvard_data_cts()
    M = optimized_M_3
    print("Test on Harvard Dataset")
    sure_list, unsure_list, neg_list, x = get_test_results(M, cts)
    print("Sure list:", sure_list)
    print("Unsure list:", unsure_list)
    print("Neg list:", neg_list)
    print("x:", x)

  config.app_algo = "COMP"
  print("\nUsing algorithm %s \n" % config.app_algo)
  harvard_test()

  config.app_algo = "SBL"
  print("\nUsing algorithm %s \n" % config.app_algo)
  harvard_test()

  config.app_algo = "combined_COMP_NNOMP_random_cv"
  print("\nUsing algorithm %s \n" % config.app_algo)
  harvard_test()
