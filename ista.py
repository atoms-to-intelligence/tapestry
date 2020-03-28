import numpy as np
import math
from sklearn.linear_model import Lasso, LassoLars, LassoCV, LassoLarsCV

from comp import create_infection_array_with_num_cases, COMP

import pylops
import spgl1 as spg
import time
from math import sqrt

from scipy import linalg

def soft_thresh(x, l):
    return np.sign(x) * np.maximum(np.abs(x) - l, 0.)


def ista(A, b, l, maxit):
    x = np.zeros(A.shape[1])
    pobj = []
    L = linalg.norm(A) ** 2  # Lipschitz constant
    time0 = time.time()
    for _ in range(maxit):
        x = soft_thresh(x + np.dot(A.T, b - A.dot(x)) / L, l / L)
        this_pobj = 0.5 * linalg.norm(A.dot(x) - b) ** 2 + l * linalg.norm(x, 1)
        pobj.append((time.time() - time0, this_pobj))

    times, pobj = map(np.array, zip(*pobj))
    return x, pobj, times

def fista(A, b, l, maxit):
    x = np.zeros(A.shape[1])
    pobj = []
    t = 1
    z = x.copy()
    L = linalg.norm(A) ** 2
    time0 = time.time()
    for _ in range(maxit):
        xold = x.copy()
        z = z + A.T.dot(b - A.dot(z)) / L
        x = soft_thresh(z, l / L)
        t0 = t
        t = (1. + sqrt(1. + 4. * t ** 2)) / 2.
        z = x + ((t0 - 1.) / t) * (x - xold)
        this_pobj = 0.5 * linalg.norm(A.dot(x) - b) ** 2 + l * linalg.norm(x, 1)
        pobj.append((time.time() - time0, this_pobj))

    times, pobj = map(np.array, zip(*pobj))
    return x, pobj, times


# Use compressed sensing to solve 0.5*||Mx - y||^2 + l * ||x||_1
class CS(COMP):
  def __init__(self, n, t, s, d, l, arr):
    super().__init__(n, t, s, d, arr)
    self.create_conc_matrix_from_infection_array(arr)
    self.l = l

  # Multiply actual conc matrix to M. This is the part done by mixing samples
  # and qpcr
  def get_quantitative_results(self, conc):
    conc = np.expand_dims(conc, axis=-1)
    #print(self.M.shape, conc.shape)
    return np.matmul(self.M.T, conc).flatten()

  # Initial concentration of RNA in each sample
  def create_conc_matrix_from_infection_array(self, arr):
    #conc = 1 + np.random.poisson(lam=5, size=self.n)
    conc = np.random.randint(low=1, high=11, size=self.n)
    #conc = np.ones(self.n)
    self.conc = conc * arr # Only keep those entries which are 

  # Solve the CS problem using Lasso
  #
  # y is results
  def decode_lasso(self, results):
    #lasso = LassoLars(alpha=self.l)
    #lasso = Lasso(alpha=self.l, max_iter=1000)
    #lasso = LassoCV(n_alphas=100)
    #lasso.fit(self.M.T, results)
    temp_mat=(self.M.T).astype(float)
    #print(type(results))
    
    
    
    #print('best lambda = ', lasso.alpha_)
    #answer = lasso.coef_
    #Using ISTA
    temp_mat=pylops.MatrixMult(temp_mat)
    answer = pylops.optimization.sparsity.FISTA(temp_mat, results, 1000, eps=self.l, tol=0, returninfo=True)[0]
    #Using OMP
    #answer = pylops.optimization.sparsity.OMP(temp_mat, results, 20000, sigma=2e-5)[0]
    #Using spgl1 lasso
    #[answer,r,g,info] = spg.spg_lasso(temp_mat,results,1e-4);
    #Using ISTA function
    #answer = fista(temp_mat, results, self.l, 10000)[0]

    score = math.sqrt(np.linalg.norm(answer - self.conc) / self.d)
    infected = (answer != 0.).astype(np.int32)

    # Compute stats
    tpos = (infected * self.arr)
    fneg = (1 - infected) * self.arr
    fpos = infected * (1 - self.arr)
    
    tp = sum(tpos)
    fp = sum(fpos)
    fn = sum(fneg)
    
    return score, tp, fp, fn

  def decode_qp(self, results):
    pass

class CSExpts:
  def __init__(self, name):
    self.name = name
    self.no_error = 0
    self.total_tp = 0
    self.total_fp = 0
    self.total_fn = 0

  # Find results using qPCR
  def do_single_expt(self, i, cs, x):
    results = cs.get_quantitative_results(x)
    score, tp, fp, fn = cs.decode_lasso(results)
    print('%s iter = %d score: %.2f' % (self.name, i, score), 'tp = ', tp, 'fp =', fp, 'fn = ', fn)
    if fp == 0 and fn == 0:
      self.no_error += 1
    self.total_tp += tp
    self.total_fp += fp
    self.total_fn += fn

  def print_stats(self, num_expts):
    precision = self.total_tp / float(self.total_tp + self.total_fp)
    recall = self.total_tp / float(self.total_tp + self.total_fn)
    print('******', self.name, 'Statistics', '******')
    print('No errors in %d / %d cases' % (self.no_error, num_expts))
    print('precision = %.3f, recall = %.3f' % (precision, recall))
    print('total tp =', self.total_tp, 'total fp = ', self.total_fp, 'total fn = ', self.total_fn)

def main():
  # Test width. Max number of parallel tests available.
  t = 384

  # Group size
  n = 1000

  # Number of infections. Sparsity
  d = 20

  # Test assignment probability. Probability that a person gets assigned to a
  # test
  s = 500. / 1000

  # lambda for regularization
  l = 500

  num_expts = 50
  quant_csexpts = CSExpts('Quant    ')
  ber_csexpts = CSExpts('Bernoulli')
  for i in range(num_expts):
    arr = create_infection_array_with_num_cases(n, d)
    cs = CS(n, t, s, d, l, arr)
    quant_csexpts.do_single_expt(i, cs, cs.conc)
    ber_csexpts.do_single_expt(i, cs, arr)
  quant_csexpts.print_stats(num_expts)
  ber_csexpts.print_stats(num_expts)


main()