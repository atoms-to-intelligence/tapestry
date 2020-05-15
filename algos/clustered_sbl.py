# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8

# Imports the precise variant of SBL using clustering method instead of
# thresholding

from inbuilt_algos import sbl

import numpy as np

# Use thresholding method 'cluster' for SBL to get very precise results
def clustered_sbl(params):
  A = params["A"]
  y = params["y"]
  t = A.shape[0]
  n = A.shape[1]
  #print(t, n)
  #print(len(y))
  assert t == len(y)
  x_est = sbl.sbl(A, y, thresholding_method='cluster') 
  res = {'x_est' : x_est}
  #print("config.noise_model = ", config.noise_model)
  return res

# clustered SBL has some variability in results. 
#
# Repeat above precise SBL many times to reduce variability
def sbl_multi(params, selection_method, n_tries, frac):
  A = params["A"]
  y = params["y"]
  t = A.shape[0]
  n = A.shape[1]
  #print(t, n)
  #print(len(y))
  assert t == len(y)
  sum_x_est = np.zeros(n, dtype=np.int32)
  and_x_est = np.ones(n, dtype=np.int32)
  for i in range(n_tries):
    tmp_x_est = sbl.sbl(A, y, thresholding_method='cluster') 
    tmp_x_est = (tmp_x_est > 0).astype(np.int32)
    sum_x_est += tmp_x_est
    and_x_est = and_x_est * tmp_x_est
  
  union_x_est = (sum_x_est > 0).astype(np.int32)
  # Only choose an x if chosen by 80% of tests
  maj_x_est = np.zeros(n)
  maj_x_est[sum_x_est >= n_tries * frac] = 1
  #res = {'x_est' : sum_x_est}
  if selection_method == 'majority':
    res = {'x_est' : maj_x_est}
  elif selection_method == 'union':
    res = {'x_est' : union_x_est}
  elif selection_method == 'intersection':
    res = {'x_est' : and_x_est}
  else:
    raise ValueError('wrong selection method: %s' % selection_method)

  #print("config.noise_model = ", config.noise_model)
  return res
