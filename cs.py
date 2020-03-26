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
    print(self.M.shape, conc.shape)
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
    score = math.sqrt(np.linalg.norm(answer - self.conc) / 20)
    infected = (answer != 0.).astype(np.int32)

    # Compute stats
    tpos = (infected * self.arr)
    fneg = (1 - infected) * self.arr
    fpos = infected * (1 - self.arr)
    
    tp = sum(tpos)
    fp = sum(fpos)
    fn = sum(fneg)
    
    print('score: %.2f' % score, 'tp = ', tp, 'fp =', fp, 'fn = ', fn)
    return tp, fp, fn

def main():
  # Test width. Max number of parallel tests available.
  t = 96

  # Group size
  n = 1000

  # Number of infections. Sparsity
  d = 8

  # Test assignment probability. Probability that a person gets assigned to a
  # test
  s = 0.5

  # lambda for regularization
  l = 0.1

  arr = create_infection_array_with_num_cases(n, d)
  cs = CS(n, t, s, d, l, arr)
  results = cs.get_quantitative_results()
  print(results.shape)
  tp, fp, fn = cs.decode_lasso(results)

main()
