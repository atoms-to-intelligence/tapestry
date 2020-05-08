import nnompcv
if __name__=='__main__':
  import sys
  # This is a hack needed so that you can directly run this file as 
  # python inbuilt_algos/nnompcv.py.
  # 
  # This is needed because we want to import "matrices" from one level up, but python
  # does not know about one level up. Please only use this hack in test files
  # and not in files which are meant to be imported from other code.
  sys.path.append(".")
  from core.matrices import *
  A=optimized_M_3
  print(np.shape(A))
  print('**')
  x = np.random.rand(A.shape[1])
  x[np.random.permutation(A.shape[1])[0:A.shape[1]-5]]=0
  y=np.matmul(A,x)
  xe=nnompcv.nnomp(A,0,y,0,8)
  print(np.nonzero(xe))
  print(np.nonzero(x))

