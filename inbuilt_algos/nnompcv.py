# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
import numpy as np
import math
from sklearn.preprocessing import normalize
from scipy.optimize import nnls

def nnomp(A, Acv, y, ycv, d, cv=False):
  r = y
  T1 = np.empty(0,dtype='int')
  #print(A.shape)
  n = A.shape[1]
  x1 = np.zeros(n)
  Au = normalize(np.transpose(A))
  epscv = np.linalg.norm(ycv)
  minepscv = math.inf
  minx = x1
  min_i = 0
  errors = []
  for i in range(d):
    a = np.matmul(Au,r)
    m = np.amax(a)
    I = np.where(a == m)
    if (m>0):
      T1 = np.union1d(T1, [I])
      np.put(x1,T1,nnls(A[: , T1], y)[0])
      r = y-np.matmul(A,x1)
      if(cv):
        rcv = ycv - np.matmul(Acv, x1)
        epscv = np.linalg.norm(rcv)
        errors.append(epscv)
        if(epscv < minepscv):
          minepscv = epscv
          minx = x1
          min_i = i
  if(cv):
    # Propagate error to higher values of d
    l = len(errors)
    if l == 0:
      errors = [0] * d
    elif l < d:
      errors.extend([errors[-1]] * (d - l))
      assert len(errors) == d

    return minx, minepscv, min_i, errors
  else:
    return x1

