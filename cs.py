import numpy as np
import math
from sklearn.linear_model import Lasso, LassoLars, LassoCV, LassoLarsCV

from comp import create_infection_array_with_num_cases, COMP

import json
import pandas as pd
import os

np.set_printoptions(precision=2)

optimized_M = np.array([[0., 1., 0., 1., 0., 1., 1., 0., 0., 0., 1., 1., 0., 1., 1., 0., 0., 0.,       0., 0., 0., 1., 0., 1., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 1., 1.,         0., 1., 1., 0.],
        [1., 0., 0., 1., 0., 0., 0., 1., 0., 1., 0., 0., 1., 1., 0., 0., 0., 1.,         1., 1., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 1.,         1., 1., 0., 0.],
        [0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 1., 1., 0., 0., 0., 0., 1., 1.,         1., 0., 1., 0., 0., 0., 1., 1., 0., 1., 1., 0., 1., 1., 1., 0., 0., 0.,         0., 1., 0., 0.],
        [1., 1., 1., 1., 0., 0., 1., 0., 0., 0., 1., 1., 0., 1., 1., 0., 0., 0.,         0., 1., 0., 0., 0., 0., 1., 0., 1., 0., 0., 1., 0., 1., 0., 0., 0., 0.,         0., 0., 1., 1.],
        [1., 1., 0., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,         1., 1., 1., 0., 1., 0., 0., 1., 1., 1., 0., 0., 1., 0., 0., 1., 0., 0.,         0., 1., 0., 1.],
        [0., 0., 1., 1., 1., 1., 1., 1., 0., 1., 0., 0., 0., 1., 1., 0., 0., 0.,         0., 0., 1., 1., 0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 1., 1., 1., 0.,         0., 0., 0., 1.],
        [0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 1., 0., 1., 0., 1., 1., 0., 0.,         1., 0., 1., 0., 0., 0., 1., 0., 0., 1., 0., 1., 1., 0., 1., 0., 0., 1.,         1., 0., 1., 0.],
        [1., 1., 1., 0., 1., 1., 0., 0., 1., 1., 0., 0., 0., 0., 0., 1., 0., 0.,         0., 1., 0., 0., 1., 0., 1., 1., 0., 0., 0., 0., 0., 1., 0., 1., 1., 0.,         1., 0., 1., 0.],
        [0., 1., 0., 0., 0., 1., 0., 1., 1., 1., 0., 1., 0., 0., 0., 1., 0., 1.,         0., 0., 0., 0., 0., 1., 0., 0., 1., 0., 1., 0., 0., 0., 1., 1., 1., 0.,         1., 1., 0., 0.],
        [0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0.,         1., 1., 0., 1., 1., 1., 1., 1., 1., 0., 0., 0., 1., 0., 0., 0., 1., 1.,         1., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 0., 1., 0., 1., 1., 0.,         1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 0., 0., 1.,         0., 0., 0., 1.],
        [1., 0., 1., 0., 0., 1., 0., 0., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1.,         0., 0., 0., 1., 0., 1., 0., 0., 1., 1., 0., 1., 0., 1., 1., 1., 0., 0.,         0., 0., 0., 1.],
        [0., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 1.,         0., 1., 1., 0., 0., 1., 1., 1., 0., 0., 1., 1., 0., 1., 0., 0., 0., 0.,         1., 0., 0., 0.],
        [0., 0., 0., 0., 1., 0., 1., 1., 0., 1., 1., 0., 1., 1., 0., 0., 1., 0.,         0., 0., 0., 1., 1., 0., 1., 0., 0., 1., 1., 1., 0., 0., 1., 0., 0., 0.,         0., 0., 1., 0.],
        [0., 1., 0., 0., 1., 0., 0., 1., 1., 0., 0., 0., 1., 0., 0., 1., 0., 1.,         0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 1., 0., 0., 1., 0., 1.,         1., 0., 0., 1.],
        [1., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 1., 1., 0., 1., 1., 1.,
          1.,         1., 0., 1., 0., 1., 0., 0., 0., 1., 1., 0., 0., 0., 0.,
          0., 1., 1., 1.,         0., 0., 0., 0.]])

print(optimized_M.shape)
#sys.exit(1)

# Use compressed sensing to solve 0.5*||Mx - y||^2 + l * ||x||_1
class CS(COMP):
  def __init__(self, n, t, s, d, l, arr, M=None):
    super().__init__(n, t, s, d, arr)
    if M is not None:
      self.M = M.T
    self.create_conc_matrix_from_infection_array(arr)
    self.l = l

  # Multiply actual conc matrix to M. This is the part done by mixing samples
  # and qpcr
  def get_quantitative_results(self, conc, add_noise=False):
    conc = np.expand_dims(conc, axis=-1)
    #print(self.M.shape, conc.shape)
    y = np.matmul(self.M.T, conc).flatten()
    sigval = 0.
    if add_noise:
      sigval = 0.01*np.median(np.absolute(y))
      y = y + np.random.normal(0., sigval, len(y))
    return y, sigval

  # Initial concentration of RNA in each sample
  def create_conc_matrix_from_infection_array(self, arr):
    #conc = 1 + np.random.poisson(lam=5, size=self.n)
    conc = np.random.randint(low=1, high=32769, size=self.n) / 32768.
    #conc = np.ones(self.n)
    self.conc = conc * arr # Only keep those entries which are 

  # Solve the CS problem using Lasso
  #
  # y is results
  def decode_lasso(self, results):
    #lasso = LassoLars(alpha=self.l)
    lasso = Lasso(alpha=self.l, max_iter=10000)
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

  def decode_lasso_for_cv(self, l, train_Ms, train_ys, test_Ms, test_ys):
    scores = []
    for train_M, train_y, test_M, test_y in zip(train_Ms, train_ys, test_Ms,
        test_ys):
      #print('Doing lasso with')
      #print(train_M.shape, train_y.shape, test_M.shape, test_y.shape)
      lasso = Lasso(alpha=l, max_iter=10000)
      lasso.fit(train_M, train_y)
      pred_y = lasso.predict(test_M)
      score = np.linalg.norm(test_y - pred_y) / len(test_y)
      scores.append(score)

    avg_score = np.average(scores)
    max_score = max(scores)
    median_score = np.median(scores)
    min_score = np.min(scores)
    return avg_score
    #return min_score

  def do_cross_validation_get_lambda(self, y, sigval):
    lambda_min = max([sigval*math.sqrt(math.log(self.n))-5,0.01]);
    lambda_max = sigval*math.sqrt(math.log(self.n))+5;
    n_step = math.ceil((lambda_max - lambda_min) / 0.01)
    ll = np.linspace(lambda_min, lambda_max, n_step)
    #a = np.linspace(0.001, 0.01, num=10)
    #ll = np.concatenate([a, 10*a, 100*a, 1000*a, 10000*a, 100000*a])
    #ll = np.concatenate([a, 10*a, 100*a])
    #ll = np.linspace(0.001, 1., 1000)

    train_Ms = []
    test_Ms = []
    train_ys = []
    test_ys = []
    M = self.M.T
    # We'll do leave one out cross-validation
    for r in range(1):
      train_M = np.delete(M, r, axis=0)
      test_M = np.expand_dims(M[r], axis=0)
      train_y = np.delete(y, r, axis=0)
      test_y = np.array([y[r]])

      train_Ms.append(train_M)
      train_ys.append(train_y)
      test_Ms.append(test_M)
      test_ys.append(test_y)

    scores = []
    for l in ll:
      score = self.decode_lasso_for_cv(l, train_Ms, train_ys, test_Ms, test_ys)
      scores.append(score)
    scores = np.array(scores)
    idx = np.argmin(scores)
    #print(idx)
    self.l = ll[idx]
    #print('lambdas = ', ll)
    #print('scores = ', scores)
    print('Choosing lambda = %.4f' % self.l, 'score = %.4f' % score)
    return self.l

  def decode_qp(self, results):
    pass

  def print_matrix(self):
    pass

  def pickle_dump(self, filename):
    pass

class CSExpts:
  def __init__(self, name):
    self.name = name
    self.no_error = 0 # number of expts with no error
    self.no_fp = 0 # number of expts with no false +ve
    self.no_fn = 0 # number of expts with no false -ve
    self.total_tp = 0
    self.total_fp = 0
    self.total_fn = 0

  # Find results using qPCR
  def do_single_expt(self, i, cs, x, cross_validation=True):
    y, sigval = cs.get_quantitative_results(cs.conc, add_noise=True)
    # Now find lambda
    if cross_validation:
      l = cs.do_cross_validation_get_lambda(y, sigval)
    score, tp, fp, fn = cs.decode_lasso(y)
    #print('%s iter = %d score: %.2f' % (self.name, i, score), 'tp = ', tp, 'fp =', fp, 'fn = ', fn)
    if fp == 0 and fn == 0:
      self.no_error += 1
    if fp == 0:
      self.no_fp += 1
    if fn == 0:
      self.no_fn += 1
    self.total_tp += tp
    self.total_fp += fp
    self.total_fn += fn

  def print_stats(self, num_expts):
    precision = self.total_tp / float(self.total_tp + self.total_fp)
    recall = self.total_tp / float(self.total_tp + self.total_fn)
    avg_fp = self.total_fp / (num_expts - self.no_fp)
    avg_fn = self.total_fn / (num_expts - self.no_fn)
    print('******', self.name, 'Statistics', '******')
    print('No errors in %d / %d cases' % (self.no_error, num_expts))
    print('No fp in %d / %d cases' % (self.no_fp, num_expts))
    print('No fn in %d / %d cases' % (self.no_fn, num_expts))
    print('precision = %.6f, recall = %.6f' % (precision, recall))
    print('total tp =', self.total_tp, 'total fp = ', self.total_fp, 'total fn = ', self.total_fn)
    print('avg fp = %.3f' % avg_fp, 'avg fn = %.3f' % avg_fn)

  def return_stats(self, num_expts):
    precision = self.total_tp / float(self.total_tp + self.total_fp)
    recall = self.total_tp / float(self.total_tp + self.total_fn)
    stats = {
        'precision'     : precision,
        'recall'        : recall,
        'total_tp'      : self.total_tp,
        'total_fp'      : self.total_fp,
        'total_fn'      : self.total_fn,
        'no_error'      : self.no_error,
        'expts'         : num_expts,
        }

    return stats

def main(n, d, t, num_expts=1000):
  # Test width. Max number of parallel tests available.
  #t = 12

  # Group size
  #n = 32

  # Number of infections. Sparsity
  #d = 1

  # Test assignment probability. Probability that a person gets assigned to a
  # test
  s = 500. / 1000

  # lambda for regularization
  l = 0.01

  cv_expts = CSExpts('CV')
  nocv_expts = CSExpts('No CV')
  for i in range(num_expts):
    arr = create_infection_array_with_num_cases(n, d)
    #cs = CS(n, t, s, d, l, arr, optimized_M)
    cs = CS(n, t, s, d, l, arr, None)
    nocv_expts.do_single_expt(i, cs, cs.conc, cross_validation=False)
    cv_expts.do_single_expt(i, cs, cs.conc, cross_validation=True)
  nocv_expts.print_stats(num_expts)
  cv_expts.print_stats(num_expts)
  return cv_expts.return_stats(num_expts)

main(40, 2, 16, num_expts=100)

def dump_to_file(filename, stats):
  df = pd.DataFrame.from_dict(stats)
  cols = ['t', 'precision', 'recall', 'total_tp', 'total_fp', 'total_fn',
      'no_error', 'expts']
  df.to_csv(filename, index=False, columns=cols)
  

def do_expts_and_dump_stats():
  tt = {
        32: range(5, 17),
        40: range(12, 30),
        50: range(12, 40),
        60: range(12, 45),
        64: range(12, 45),
        70: range(14, 45),
        80: range(14, 50),
        90: range(14, 50),
        96: range(14, 60),
        100: range(14, 60),
        128: range(14, 64)
      }

  #dd = [1, 2]
  dd = [2]
  #nn = [32, 64, 96, 128]
  #nn= [60, 70, 80, 90]
  nn = [32]

  stats = {}
  algorithm = 'lasso'

  stats_dir = 'stats/%s' % algorithm
  if not os.path.exists(stats_dir):
    os.makedirs(stats_dir)

  for n in nn:
    stats[n] = {}
    for d in dd:
      stats[n][d] = {}
      t_stats = []
      for t in tt[n]:
        print('n = %d, d = %d, t = %d' % (n, d, t))
        item = main(n, d, t)
        item['t'] = t
        stats[n][d][t] = item
        t_stats.append(item)

      # Now print these stats to file
      filename = './stats/%s/stats_n_%d_d_%d.csv' % (algorithm, n, d)
      dump_to_file(filename, t_stats)

  print('Viable scenarios with no errors:')
  for n in nn:
    for d in dd:
      for t in tt[n]:
        item = stats[n][d][t]
        if item['precision'] > 0.999 and item['recall'] > 0.9999:
          print(n, d, item )

