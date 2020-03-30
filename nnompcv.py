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
        if(epscv < minepscv):
          minepscv = epscv
          minx = x1
  if(cv):
    return minx
  else:
    return x1

if __name__=='__main__':
  from matrices import *
  A=optimized_M_3
  print(np.shape(A))
  print('**')
  x = np.random.rand(A.shape[1])
  x[np.random.permutation(A.shape[1])[0:A.shape[1]-5]]=0
  y=np.matmul(A,x)
  xe=nnomp(A,0,y,0,8)
  print(np.nonzero(xe))
  print(np.nonzero(x))





