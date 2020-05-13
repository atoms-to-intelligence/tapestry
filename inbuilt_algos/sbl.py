# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
if __name__=='__main__':
  raise ValueError('This is a library file. Import this in a test file and run.')

import numpy as np
import math
import kmeans1d

# thresholding_method is either 'tau' or 'cluster'
def sbl(A, y, eps=1e-3, thresholding_method='tau'):
  '''
  Refer to SLide 31 of http://math.iisc.ernet.in/~nmi/Chandra%20Murthy.pdf for algorithm
  Inputs:
  A = measurement matrix
  y = measurements
  sigval = variance of noise in measurement
  tau = threshold on signal
  '''
  # Doing this preprocessing inside sbl() function itself
  tau = 0. #0.01 * np.min(y/np.sum(A, axis=-1))

  y_max = np.max(y)
  assert y_max >= 0
  if y_max > 0:
    A = A / y_max
    y = y / y_max

    pos_y = y[y>0.]
    pos_A = A[y>0.]
    sigval = np.std(pos_y/np.sum(pos_A, axis=-1))
  else:
    # sigma should be 0 but this will mess with the algo. Set it to some small
    # value
    sigval = 0.1
  y = np.array(y, dtype=np.float64)
  A = np.array(A, dtype=np.float64)

  [m,n] = A.shape
  
  # Pre-Processing
  # As all A_ij and x_j are positive, for any y_i=0 implies that for all j s.t A_ij=1, x_j=0. This reduces problem dimension.
  #nz_index = (np.where(y != 0))[0]
  #z_index = (np.where(y == 0))[0]
  #red_y = y[nz_index]

  #[r,c] = np.where(A[z_index,:] != 0)
  #ind_zero_x = np.unique(c)
  #ind = np.arange(0, n)
  #ind_nonzero_x = np.setxor1d(ind,ind_zero_x)
  #red_x = x[ind_nonzero_x]
  #red_A = (A[nz_index,:])[:,ind_nonzero_x]
  #red_n = (ind_nonzero_x.shape)[0]
  
  red_A = A
  red_n = n
  red_y = y
  ind_nonzero_x = np.arange(n)

  # Sparse Bayesian Learning
  # Corner cases are when 0 samples or all ys are 0
  if red_n == 0 or np.all(red_y == 0):
    x_est = np.zeros(n)
  else:
    #E-step
      #   mu is estimated mean of posterior distribution x|y, and so is the estimated red_x computed iteratively
      #   Sigma is variance of posterior distribution of x|y
    # M-step
      #   Gamma is the prior variance of x, inv_Gamma is saved as E-step only requires the inverse
    inv_Gamma = np.identity(red_n)
    gamma = np.ones(n)
    mu_old = np.ones(red_n)
    mu = np.zeros(red_n)
    variance = sigval * sigval
    #print('inside else')
    while np.sum(mu) == 0 or np.linalg.norm(mu_old - mu, 2) / np.linalg.norm(mu, 2) > eps:
      #print('inside loop')
      mu_old = mu;
      inv_Sigma = np.matmul(np.transpose(red_A), red_A) / (variance) + inv_Gamma
      Sigma = np.linalg.inv(inv_Sigma)
      mu = np.matmul(Sigma, np.matmul(np.transpose(red_A), red_y)) / (variance)
      #mu[mu<0] = 0
      err = red_y-np.dot(red_A,mu)
      variance = (np.sum(err*err) + variance * np.sum(1.0 - (1.0/gamma) * np.diag(Sigma))) / m
      gamma = np.square(mu) + np.diag(Sigma)
      inv_Gamma = np.diag(1 / gamma)

    #pos_mu = mu[mu>0]
    #min_pos_mu = np.min(pos_mu)
    #mu[mu<=0] = min_pos_mu
    #log_mu = np.log(mu)
    #clusters, centroids = kmeans1d.cluster(log_mu, 2)
    #lower_centroid = centroids[0]
    #clusters = np.array(clusters)
    #lower_cluster = mu[clusters == 0]
    #lower_cluster_std = np.std(lower_cluster)
    #tau = lower_centroid - lower_cluster_std
    #log_tau = lower_centroid - lower_cluster_std
    #mu[log_mu < log_tau] = 0
    #print('\nmu = ', mu)
    #print('\nlog_mu = ', np.log(mu))
    #print('clusters = ', clusters)
    #print('centroids = ', centroids)

    assert thresholding_method in [ 'tau', 'cluster' ]

    if thresholding_method == 'cluster':
      clusters, centroids = kmeans1d.cluster(mu, 2)
      mu = mu * np.array(clusters)

    x_est = np.zeros(n)
    x_est[ind_nonzero_x] = mu
    if thresholding_method == 'tau':
      x_est[x_est < tau] = 0

  #print(x_est)
  return x_est
