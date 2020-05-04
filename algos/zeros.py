# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8

import config
import numpy as np

# Params is a dictionary containing the following
def zeros(params):
  A = params["A"]
  y = params["y"]
  t = A.shape[0]
  n = A.shape[1]
  #print(t, n)
  #print(len(y))
  assert t == len(y)

  res = {'x_est' : np.zeros(n)}
  #print("config.noise_model = ", config.noise_model)
  return res

