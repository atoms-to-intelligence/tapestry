# Code to generate STS based matrices
import numpy as np
def sts(S):
	T = int(S*(S-1)//6)
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
	return A
	

