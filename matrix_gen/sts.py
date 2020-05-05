# Code to generate STS based matrices
import numpy as np
def sts(S, n=None):
  assert S % 3 == 0
  assert (S - 3) % 6 == 0
  T = int(S*(S-1)//6)
  if n is None:
    n = T
  assert n <= T
  K = S//3
  k = (K-1)//2
  s = int(S)
  arr = np.arange(K)
  mat = np.zeros((K,K))
  for i in range(K):
    mat[:,i] = np.copy(arr)
    arr = np.roll(arr,-1)
  mat = mat.astype('int')
  lis = [(mat[iter,iter],iter) for iter in range(K)]
  d = dict(lis)
  l = []
  for i in range(K):
    l.append([i*3,i*3+1,i*3+2])
  for i in range(K):
    for j in range(i+1,K):
      l.append([i*3,j*3,d[int(mat[i,j])]*3+1])
      l.append([i*3+1,j*3+1,d[int(mat[i,j])]*3+2])
      l.append([i*3+2,j*3+2,d[int(mat[i,j])]*3])
  A = np.zeros((s,T))
  c = 0
  for el in l:
    for e in el:
      A[e,c] = 1
    c = c+1
  A = A*1.0
  return A[:, :n]
  
if __name__ == '__main__':
  A = sts(93)
  print(A.shape)
  print(A)

  A = sts(93, 960)
  print(A.shape)
  print(A)
