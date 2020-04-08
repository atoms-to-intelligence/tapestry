import l1_ls_nonneg
import numpy as np
import math

def l1ls_cv(A, y, sigval, tau):
  
  [m,n] = A.shape
  mr = math.ceil(0.9*m);
  mcv = m-mr;
  
  yr = y[0:mr-1]
  ycv = y[mr:m-1]
  Ar = A[0:mr-1,:]
  Acv = A[mr:m,:]

  lambda_min = max(sigval * math.sqrt(math.log(n)) - 5, 0)
  lambda_max = sigval * math.sqrt(math.log(n)) + 5
  lambdas = np.arange(lambda_min, lambda_max, 0.05)
  l_lambdas = len(lambdas)
  cv_error = np.zeros((l_lambdas, 1))

  for i in range(l_lambdas):
    [x_est, status, hist] = l1_ls_nonneg.l1ls_nonneg(Ar, yr, lambdas[i], tar_gap=0.001, quiet=1)
    cv_error[i] = np.linalg.norm(ycv - np.matmul(Acv, x_est), 2) / mcv

  minind = np.argmin(cv_error)
  [x_est, status, hist] = l1_ls_nonneg.l1ls_nonneg(Ar, yr, lambdas[minind], tar_gap=0.001, quiet=1)
  x_est[x_est < tau] = 0

  return x_est

def l1ls(A, y, l, tau):
  [x_est, status, hist] = l1_ls_nonneg.l1ls_nonneg(A, y, l, tar_gap=0.001, quiet=1)
  x_est[x_est < tau] = 0

  return x_est
