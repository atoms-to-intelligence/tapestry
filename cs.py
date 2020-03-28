import numpy as np
import math
from sklearn.linear_model import Lasso, LassoLars, LassoCV, LassoLarsCV

from comp import create_infection_array_with_num_cases, COMP

import json
import pandas as pd
import os

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

  def do_cross_validation_get_lambda(self, results):
    a = np.linspace(0.0001, 0.0010, num=10)
    print(a)
    ll = np.concatenate(a, 10*a, 100*a, 1000*a, 10000*a)
    for l in ll:
      pass

  def decode_qp(self, results):
    pass

  def print_matrix(self):
    pass

  def pickle_dump(self, filename):
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
    #print('%s iter = %d score: %.2f' % (self.name, i, score), 'tp = ', tp, 'fp =', fp, 'fn = ', fn)
    if fp == 0 and fn == 0:
      self.no_error += 1
    self.total_tp += tp
    self.total_fp += fp
    self.total_fn += fn

  def print_stats(self, num_expts):
    precision = self.total_tp / float(self.total_tp + self.total_fp)
    recall = self.total_tp / float(self.total_tp + self.total_fn)
    #print('******', self.name, 'Statistics', '******')
    print('No errors in %d / %d cases' % (self.no_error, num_expts))
    print('precision = %.6f, recall = %.6f' % (precision, recall))
    print('total tp =', self.total_tp, 'total fp = ', self.total_fp, 'total fn = ', self.total_fn)

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

def main(n, d, t):
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

  num_expts = 1000
  quant_csexpts = CSExpts('Quant    ')
  ber_csexpts = CSExpts('Bernoulli')
  for i in range(num_expts):
    arr = create_infection_array_with_num_cases(n, d)
    cs = CS(n, t, s, d, l, arr)
    quant_csexpts.do_single_expt(i, cs, cs.conc)
    #ber_csexpts.do_single_expt(i, cs, arr)
  #quant_csexpts.print_stats(num_expts)
  return quant_csexpts.return_stats(num_expts)
  #ber_csexpts.print_stats(num_expts)

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
        #print(item)
        stats[n][d][t] = item
        t_stats.append(item)
        #printable = json.dumps(item, indent=4)
        #print(printable)
        #stats['(n=%d, d=%d, t=%d)' % (n, d, t)] = item

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
  #printable = json.dumps(stats, indent=4)
  #print(printable)

