import numpy as np
import math
from sklearn.linear_model import Lasso, LassoLars, LassoCV, LassoLarsCV

from comp import create_infection_array_with_num_cases, COMP

# Use compressed sensing to solve 0.5*||Mx - y||^2 + l * ||x||_1
class CS(COMP):
  def __init__(self, n, t, s, d, l, arr):
    super().__init__(n, t, s, d, arr)
    self.create_conc_matrix_from_infection_array(arr)
    self.l = l

  # Multiply actual conc matrix to M. This is the part done by mixing samples
  # and qpcr
  def get_quantitative_results(self):
    conc = np.expand_dims(self.conc, axis=-1)
    #print(self.M.shape, conc.shape)
    return np.matmul(self.M.T, conc).flatten()

  # Initial concentration of RNA in each sample
  def create_conc_matrix_from_infection_array(self, arr):
    conc = 1 + np.random.poisson(lam=5, size=self.n)
    self.conc = conc * arr # Only keep those entries which are 

  # Solve the CS problem using Lasso
  #
  # y is results
  def decode_lasso(self, results):
    #lasso = LassoLars(alpha=self.l)
    lasso = Lasso(alpha=self.l, max_iter=1000)
    #lasso = LassoCV(n_alphas=100)
    lasso.fit(self.M.T, results)
    #print('best lambda = ', lasso.alpha_)
    answer = lasso.coef_
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

def main():
  # Test width. Max number of parallel tests available.
  t = 96*4

  # Group size
  n = 1000

  # Number of infections. Sparsity
  d = 40

  # Test assignment probability. Probability that a person gets assigned to a
  # test
  s = 0.5

  # lambda for regularization
  l = 0.1

  no_error = 0
  num_expts = 1000
  total_tp = 0
  total_fp = 0
  total_fn = 0
  for i in range(num_expts):
    arr = create_infection_array_with_num_cases(n, d)
    cs = CS(n, t, s, d, l, arr)
    results = cs.get_quantitative_results()
    score, tp, fp, fn = cs.decode_lasso(results)
    print('iter = %d score: %.2f' % (i, score), 'tp = ', tp, 'fp =', fp, 'fn = ', fn)
    if fp == 0 and fn == 0:
      no_error += 1
    total_tp += tp
    total_fp += fp
    total_fn += fn

  precision = total_tp / float(total_tp + total_fp)
  recall = total_tp / float(total_tp + total_fn)
  print('No errors in %d / %d cases' % (no_error, num_expts))
  print('precision = %.3f, recall = %.3f' % (precision, recall))
  print('total fp = ', total_fp, 'total fn = ', total_fn)

main()
