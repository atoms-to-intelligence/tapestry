import numpy as np
import math
from sklearn.preprocessing import normalize
from scipy.optimize import nnls

def nnomp(A, Acv, y, ycv, d):
	r = y
	T1 = np.empty(0,dtype='int')
	#print(A.shape)
	n = A.shape[1]
	x1 = np.zeros(n)
	Au = normalize(np.transpose(A))
	for i in range(d):
		a = np.matmul(Au,r)
		m = np.amax(a)
		I = np.where(a == m)
		if (m>0):
			T1 = np.union1d(T1, [I])
			np.put(x1,T1,nnls(A[: , T1], y)[0])
			r = y-np.matmul(A,x1)
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





