import numpy as np
import math

def sbl(A, y, sigval, tau, eps=1e-4):
  '''
  Refer to SLide 31 of http://math.iisc.ernet.in/~nmi/Chandra%20Murthy.pdf for algorithm
  Inputs:
  A = measurement matrix
  y = measurements
  sigval = variance of noise in measurement
  tau = threshold on signal
  '''
  [m,n] = A.shape
  
  # Pre-Processing
  # As all A_ij and x_j are positive, for any y_i=0 implies that for all j s.t A_ij=1, x_j=0. This reduces problem dimension.
  nz_index = (np.where(y != 0))[0]
  z_index = (np.where(y == 0))[0]
  red_y = y[nz_index]

  [r,c] = np.where(A[z_index,:] != 0)
  ind_zero_x = np.unique(c)
  ind = np.arange(0, n)
  ind_nonzero_x = np.setxor1d(ind,ind_zero_x)
  #red_x = x[ind_nonzero_x]
  red_A = (A[nz_index,:])[:,ind_nonzero_x]
  red_n = (ind_nonzero_x.shape)[0]

  # Sparse Bayesian Learning
  if red_n == 0:
    x_est = np.zeros(n)
    #print('inside if')
  else:
    #E-step
      #   mu is estimated mean of posterior distribution x|y, and so is the estimated red_x computed iteratively
      #   Sigma is variance of posterior distribution of x|y
    # M-step
      #   Gamma is the prior variance of x, inv_Gamma is saved as E-step only requires the inverse
    inv_Gamma = np.identity(red_n)
    mu_old = np.ones(red_n)
    mu = np.zeros(red_n)

    #print('inside else')
    while np.sum(mu) == 0 or np.linalg.norm(mu_old - mu, 2) / np.linalg.norm(mu, 2) > eps:
      #print('inside loop')
      mu_old = mu;
      inv_Sigma = np.matmul(np.transpose(red_A), red_A) / (sigval*sigval) + inv_Gamma
      Sigma = np.linalg.inv(inv_Sigma)
      mu = np.matmul(Sigma, np.matmul(np.transpose(red_A), red_y)) / (sigval*sigval)
      gamma = np.square(mu) + np.diag(Sigma)
      inv_Gamma = np.diag(1 / gamma)

    x_est = np.zeros(n)
    x_est[ind_nonzero_x] = mu
    x_est[x_est < tau] = 0

  #print(x_est)
  return x_est
