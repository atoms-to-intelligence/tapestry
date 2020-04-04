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
    self.uncon_negs = 0

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
    if algo == 'COMP':
      infected, score, tp, fp, fn = cs.decode_comp_new(bool_y)
      uncon_negs = 0
    else:
      infected, prob1, prob0, score, tp, fp, fn, uncon_negs = cs.decode_lasso(y, algo,
          prefer_recall=False)

    sys.stdout.write('\riter = %d / %d score: %.2f tp = %d fp = %d fn = %d' %
        (i, num_expts, score, tp, fp, fn))
    sys.stdout.flush()

    self.add_stats(tp, fp, fn, uncon_negs)

  def add_stats(self, tp, fp, fn, uncon_negs):
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
    self.uncon_negs += uncon_negs

  def print_stats(self, num_expts, header=False):
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

    if header:
      print('******', self.name, 'Statistics', '******')
    print('No errors in %d / %d cases' % (self.no_error, num_expts))
    print('No fp in %d / %d cases' % (self.no_fp, num_expts))
    print('No fn in %d / %d cases' % (self.no_fn, num_expts))
    print('precision = %.6f, recall = %.6f' % (precision, recall))
    print('total tp =', self.total_tp, 'total fp = ', self.total_fp,
        'total fn = ', self.total_fn, 'unconfident negs = ', self.uncon_negs)
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
  s = min(1 / d, 0.5)
  #s = 0.5

  # lambda for regularization
  l = 0.01

  if xs is not None:
    assert num_expts == len(xs)

  if isinstance(algo, list):
    #print('list')
    expts = [CSExpts(item) for item in algo]
    ret_list = True
  else:
    #print('str')
    expts = [CSExpts(algo)]
    algo = [algo]
    ret_list = False

  if mr is None:
    prntstr = ('\nn = %d, d = %d, t = %d\n' % (n, d, t))
  else:
    prntstr = ('\nn = %d, d = %d, t = %d, mr = %d\n' % (n, d, t, mr))
  print(prntstr)
  try:
    for i in range(num_expts):
      if xs is not None:
        arr = xs[i]
      else:
        arr = create_infection_array_with_num_cases(n, d)
      #cs = CS(n, t, s, d, l, x, optimized_M)
      cs = CS(n, t, s, d, l, arr, M, mr)
      for expt, alg in zip(expts, algo):
        expt.do_single_expt(i, num_expts, cs, cs.conc,
            algo=alg,
            cross_validation=cross_validation,
            add_noise=add_noise,
            noise_magnitude=noise_magnitude)
    print('')
    for expt in expts:
      expt.print_stats(num_expts)
  except KeyboardInterrupt:
    for expt in expts:
      expt.print_stats(i)
    raise

  stats = []
  for expt in expts:
    expt.n = n
    expt.d = d
    expt.t = t
    expt.mr = mr
    stat = expt.return_stats(num_expts)
    stat['n'] = n
    stat['d'] = d
    stat['t'] = t
    stats.append(stat)

  if ret_list:
    return expts, stats
  else:
    return expts[0], stats[0]

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
  num_expts = 100
  n = 500
  t = 46
  #matrix = np.random.binomial(1, 0.5, size=(t, n))
  #matrix = optimized_M_16_64_1
  matrix = optimized_M_46_500_1
  #matrix = None

  algos = ['combined_COMP_NNOMP_random_cv']
  #algos = ['combined_COMP_NNOMP_random_cv']
  #algos = ['combined_COMP_NNOMP_random_cv',
  #    'NNOMP_random_cv']
  #algos = [ 'NNOMP', ]
  #algos = ['combined_COMP_NNOMP_random_cv', 'SBL']
  #algos = ['combined_COMP_NNOMP_random_cv',
  #    'NNOMP_random_cv']
  add_noise = True
  d_range = range(1, 11)
  retvals = Parallel(n_jobs=10)\
  (\
      delayed(do_many_expts)\
      (
        n, d, t, num_expts=num_expts, M=matrix,\
        add_noise=add_noise,algo=algos, mr=None \
      )\
      for d in d_range\
  )

  l = len(algos)
  statlist = [[] for i in range(l)]
  idx = range(l)
  for item in retvals:
    expts, stats = item
    for i, expt, stat in zip(idx, expts, stats):
      statlist[i].append(stat)

  print('\nn = %d, t = %d\n' % (expt.n, expt.t))
  for i, algo in enumerate(algos):
    print('\n' + algo + '\n')
    print('\td\tPrecision\tRecall\ttotal_tests\n')
    for stat in statlist[i]:
      total_tests = t + stat['d'] / stat['precision']
      print('\t%d\t%.3f\t\t%.3f\t%.1f' % (stat['d'], stat['precision'],
        stat['recall'], total_tests))

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


def small_test_decode_comp_combined():
  A = np.array(
      [
        [1, 0, 0, 1, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 1, 1],
      ])

  x = np.array([3, 0, 2, 9, 1, 0])
  arr = (x > 0).astype(np.int32)
  y = np.matmul(A, x)

  print('A = ')
  print(A)
  print('x = ', x)
  print('y = ', y)

  n = A.shape[1]
  t = A.shape[0]
  d = 1
  s = 0.5
  l = 0.1
  mr = None
  cs = CS(n, t, s, d, l, arr, A, mr)
  # conc is computed from arr by CS. Set it here manually for testing
  cs.conc = x 
  results, _ = cs.get_quantitative_results(x, add_noise=False)
  print(y)
  print(results)
  assert np.all(y == results)

  infected, score, tp, fp, fn = cs.decode_comp_combined(y, 'NNOMP')
  print('infected:', infected)
  print('tp: %d, fp: %d, fn: %d' % (tp, fp, fn))

  infected, score, tp, fp, fn = cs.decode_lasso(y, 'NNOMP')
  print('without comp infected:', infected)
  print('without comp tp: %d, fp: %d, fn: %d' % (tp, fp, fn))


def large_test_decode_comp_combined(num_expts):
  A = optimized_M_2
  n = A.shape[1]
  t = A.shape[0]
  d = 2
  s = 0.5
  l = 0.01
  mr = None

  algo = 'NNOMP_random_cv'
  with_comp = CSExpts('With COMP and %s' % algo)
  without_comp = CSExpts('Only %s' % algo)
  for i in range(num_expts):
    arr = create_infection_array_with_num_cases(n, d)
    cs = CS(n, t, s, d, l, arr, A, mr)
    y, _ = cs.get_quantitative_results(cs.conc, add_noise=True)
    #print(y)

    #infected, score, tp, fp, fn = cs.decode_comp_combined(y, 'NNOMP_random_cv')
    infected, score, tp, fp, fn = cs.decode_comp_combined(y, algo,
        test=True)
    with_comp.add_stats(tp, fp, fn)
    #print('infected:', infected)
    #print('score %.4f, tp: %d, fp: %d, fn: %d' % (score, tp, fp, fn))

    #infected, score, tp, fp, fn = cs.decode_lasso(y, 'NNOMP_random_cv')
    infected, score, tp, fp, fn = cs.decode_lasso(y, algo)
    without_comp.add_stats(tp, fp, fn)
    #print('without comp infected:', infected)
    #print('without comp score %.4f, tp: %d, fp: %d, fn: %d' % (score, tp, fp, fn))

  with_comp.print_stats(num_expts, header=True)
  without_comp.print_stats(num_expts, header=True)

if __name__=='__main__':
  #large_test_decode_comp_combined(1000)
  #mr = None
  #do_many_expts(200, 6, 46, num_expts=100, M=None,
  #    add_noise=True,algo='combined_COMP_NNOMP_random_cv', mr=mr)
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


