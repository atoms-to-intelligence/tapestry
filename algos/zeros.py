# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8

# Sample algorithm file

# config.py is one directory up. You can import it and other needed files from
# the parent directory like this. Typically you won't need anything
from core import config

# Any standard import works as usual
import numpy as np

# params is a dictionary containing the following
#
# "A" - the matrix A. t x n numpy.ndarray
# "y" - the (noisy) vector of results. length t nump.ndarray
# 
# The number of tests 't' and number of samples 'n' are not explicitly
# provided. They are to be inferred from the matrix shape.
#
# Return: dictionary res containing the following:
# "x_est" - numpy vector of length n
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

# If your algorithm has some hyper-parameters which should be cross-validated,
# please put the cross-validation code here. Signature of this is identifcal
# to that of above. If your algo name is "FOO" then add this as "FOO_cv" in
# __init__.py in the same folder
def zeros_cv(params):
  res = {'x_est' : np.zeros(n)}
  #print("config.noise_model = ", config.noise_model)
  return res

