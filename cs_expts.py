from cs import *

import json
import pandas as pd
import os


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
  def do_single_expt(self, i, num_expts, cs, x, cross_validation=True, add_noise=True,
      algo='OMP', noise_magnitude=None):

    y, sigval = cs.get_quantitative_results(cs.conc, add_noise=add_noise,
        noise_magnitude=noise_magnitude)
    bool_y = (y > 0).astype(np.int32)

    # Now find lambda
    if cross_validation and algo == 'lasso':
      l = cs.do_cross_validation_get_lambda(y, sigval)
    elif cross_validation:
      raise ValueError('No cross validation implemented for %s' % algo)
    if algo == 'OMP' or algo == 'lasso' or algo == 'NNOMP' or algo == 'NNOMPCV' \
        or algo == 'NNOMP_loo_cv' or algo == 'NNOMP_random_cv':
      score, tp, fp, fn = cs.decode_lasso(y, algo)
    if algo == 'COMP':
      score, tp, fp, fn = cs.decode_comp_new(bool_y)

    sys.stdout.write('\riter = %d / %d score: %.2f tp = %d fp = %d fn = %d' %
        (i, num_expts, score, tp, fp, fn))
    sys.stdout.flush()
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
    if num_expts == self.no_fp:
      avg_fp = 0
    else:
      avg_fp = self.total_fp / (num_expts - self.no_fp)
    if num_expts == self.no_fn:
      avg_fn = 0
    else:
      avg_fn = self.total_fn / (num_expts - self.no_fn)
    #print('******', self.name, 'Statistics', '******')
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

def do_many_expts(n, d, t, num_expts=1000, xs=None, M=None,
    cross_validation=False,
    add_noise=False,
    algo='OMP',
    noise_magnitude=None,
    mr=None):
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
  l = 1000.

  if xs is not None:
    assert num_expts == len(xs)

  #cv_expts = CSExpts('CV')
  nocv_expts = CSExpts('No CV')
  if mr is None:
    prntstr = ('\nn = %d, d = %d, t = %d\n' % (n, d, t))
  else:
    prntstr = ('\nn = %d, d = %d, t = %d, mr = %d\n' % (n, d, t, mr))
  print(prntstr)
  try:
    for i in range(num_expts):
      if xs is not None:
        x = xs[i]
      else:
        x = create_infection_array_with_num_cases(n, d)
      #cs = CS(n, t, s, d, l, x, optimized_M)
      cs = CS(n, t, s, d, l, x, M, mr)
      nocv_expts.do_single_expt(i, num_expts, cs, cs.conc,
          algo=algo,
          cross_validation=cross_validation,
          add_noise=add_noise,
          noise_magnitude=noise_magnitude)
      #cv_expts.do_single_expt(i, cs, cs.conc, cross_validation=True)
    print('')
    nocv_expts.print_stats(num_expts)
  except KeyboardInterrupt:
    nocv_expts.print_stats(i)
    raise

  nocv_expts.n = n
  nocv_expts.d = d
  nocv_expts.t = t
  nocv_expts.mr = mr
  #cv_expts.print_stats(num_expts)
  stat = nocv_expts.return_stats(num_expts)
  stat['n'] = n
  stat['d'] = d
  stat['t'] = t
  #return stat
  return nocv_expts

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
  algorithm = 'NNOMPCV'

  stats_dir = 'stats/%s' % algorithm
  if not os.path.exists(stats_dir):
    os.makedirs(stats_dir)

  for n in nn:
    stats[n] = {}
    for d in dd:
      stats[n][d] = {}
      t_stats = []
      for t in tt[n]:
        print("")
        print('n = %d, d = %d, t = %d ' % (n, d, t))
        print(" ")
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

def run_many_parallel_expts():
  num_expts = 10
  n = 40
  t = 16
  matrix = optimized_M_2
  mr = 11
  expts = Parallel(n_jobs=1)\
  (\
      delayed(do_many_expts)\
      (
        n, d, t, num_expts=num_expts, M=matrix,\
        add_noise=True,algo='NNOMP_random_cv', mr=mr \
      )\
      for d in range(1, 5)\
  )

  for expt in expts:
    prntstr = ('\nn = %d, d = %d, t = %d\n' % (expt.n, expt.d, expt.t))
    print(prntstr)
    expt.print_stats(num_expts)

def run_many_parallel_expts_mr():
  num_expts = 1000
  expts = Parallel(n_jobs=4)\
  (\
      delayed(do_many_expts)\
      (
        60, 2, 24, num_expts=num_expts, M=optimized_M_3,\
        add_noise=True,algo='NNOMP_random_cv', mr=mr \
      )\
      for mr in range(14, 23)\
  )

  for expt in expts:
    prntstr = ('\nn = %d, d = %d, t = %d, mr = %d\n' % (expt.n, expt.d, expt.t, expt.mr))
    print(prntstr)
    expt.print_stats(num_expts)


if __name__=='__main__':
  #do_many_expts(400, 5, 64, num_expts=100, M=optimized_M_5,
  #    add_noise=True,algo='NNOMP_random_cv', mr=45)
  run_many_parallel_expts()
  #for mr in range(8, 15):
  #  do_many_expts(40, 2, 16, num_expts=1000, M=optimized_M_2,
  #      add_noise=True,algo='NNOMP_random_cv', mr=mr)
  #do_many_expts(40, 4, 16, num_expts=1000, M=optimized_M_2,
  #    add_noise=True,algo='NNOMP_loo_cv')
  #do_many_expts(40, 3, 16, num_expts=1000, M=optimized_M_2,
  #    add_noise=True,algo='NNOMP')

# Following code assumes the same Gaussian noise for each y
#for noise in [0.0001, 0.0002, 0.0004, 0.0008, 0.001, 0.002, 0.004, 0.008, 0.01]:
#  do_many_expts(40, 2, 16, num_expts=1000, M=optimized_M_2, add_noise=True,
#      cross_validation=False, algo='OMP', noise_magnitude=noise)


