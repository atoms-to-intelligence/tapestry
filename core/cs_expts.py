# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
if __name__=='__main__':
  raise ValueError('Use tools/run_expts.py. This is a library file.')

from matrix_gen import sts
from core.cs import *
from utils.pickle_manager import stats_manager
from utils import stats_helper

import json
import pandas as pd
import os

def specificity(precision, recall, n, d, preds):
  if preds == 0 or (precision == 0 and recall ==0):
    return 1
  else:
    assert precision != 0
    assert recall != 0
  return 1 - (recall * d * ( (1 / precision) - 1) / (n - d)) 

def compute_rmse(x, x_est):
  if np.all(x == 0):
    # if x is 0 vector, then RMSE is either 0 if x_est is also 0 or considered to
    # be 1 if x_est is not 0
    if np.all(x_est == 0):
      return 0.
    else:
      return 1.
  return np.linalg.norm(x - x_est) / np.linalg.norm(x)

class SingleExpt:
  def __init__(self, tp, fp, fn, uncon_negs, determined, overdetermined, surep,
      unsurep, wrongly_undetected, score, rmse, x, bool_x, y, bool_y, x_est,
      infected, infected_dd, n, d, t, mr, algo, mlabel):
    # computed stats
    self.tp = tp
    self.fp = fp
    self.fn = fn
    self.uncon_negs = uncon_negs
    self.determined = determined
    self.overdetermined = overdetermined
    self.surep = surep
    self.unsurep = unsurep
    self.wrongly_undetected = wrongly_undetected
    self.score = score
    self.rmse = rmse

    # Experimental setting and returned value
    self.x_idx = np.array(np.where(x), dtype=np.uint16)
    self.x = np.array(x[self.x_idx], dtype=np.float32)
    self.y_idx = np.array(np.where(y), dtype=np.uint16)
    self.y = np.array(y[self.y_idx], dtype=np.float32)

    self.x_est_idx = np.array(np.where(x_est), dtype=np.uint16)
    self.x_est = np.array(x_est[self.x_est_idx], dtype=np.float32)

    self.infected_idx = np.array(np.where(infected), dtype=np.uint16)
    self.infected_dd_idx = np.array(np.where(infected_dd), dtype=np.uint16)

    # Global settings common for many expts
    self.n = n
    self.d = d
    self.t = t
    self.mr = mr
    self.algo = algo
    self.mlabel = mlabel

    # Not saving the matrix label here because it'd need changes at too many
    # places. We'll have the matrix label as the first key into the dict
    # anyway

class CSExpts:
  def __init__(self, name, n, d, t, mr, algo, mlabel):
    self.n = n
    self.d = d
    self.t = t
    self.mr = mr
    self.name = name
    self.algo = algo
    self.mlabel = mlabel

    self.no_error = 0 # number of expts with no error
    self.no_fp = 0 # number of expts with no false +ve
    self.no_fn = 0 # number of expts with no false -ve
    self.total_tp = 0
    self.total_fp = 0
    self.total_fn = 0
    self.uncon_negs = 0
    self.determined = 0
    self.overdetermined = 0
    self.surep = 0
    self.unsurep = 0
    self.num_expts_2_stage = 0
    self.wrongly_undetected = 0
    self.total_score = 0
    self.total_rmse = 0
    self.min_precision = np.inf
    self.all_precisions = []
    self.all_recalls = []
    self.all_fps = []
    self.all_fns = []
    self.single_expts = []

  # manage stats for a single expt
  def record_single_expt(self, tp, fp, fn, uncon_negs, determined, overdetermined, surep,
      unsurep, wrongly_undetected, score, rmse, x, bool_x, y, bool_y, x_est,
      infected, infected_dd):
    single_expt = SingleExpt(tp, fp, fn, uncon_negs, determined, overdetermined, surep,
      unsurep, wrongly_undetected, score, rmse, x, bool_x, y, bool_y, x_est,
      infected, infected_dd, self.n, self.d, self.t, self.mr, self.algo, self.mlabel)

    self.single_expts.append(single_expt)

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
      infected, infected_dd, score, tp, fp, fn, surep, unsurep,\
          num_infected_in_test = cs.decode_comp_new(bool_y)
      uncon_negs = 0
      determined = 0
      overdetermined = 0
      wrongly_undetected = 0
      x_est = np.zeros(cs.n)
      rmse = 1.
    elif algo == 'SCOMP':
      infected, infected_dd, score, tp, fp, fn, surep, unsurep,\
          num_infected_in_test = cs.decode_scomp(bool_y)
      uncon_negs = 0
      determined = 0
      overdetermined = 0
      wrongly_undetected = 0
      x_est = np.zeros(cs.n)
    else:
      x_est, infected, infected_dd, prob1, prob0, score, tp, fp, fn, uncon_negs, determined,\
          overdetermined, surep, unsurep, wrongly_undetected, \
          num_infected_in_test = cs.decode_lasso(y, algo,
          prefer_recall=False)
      rmse = compute_rmse(cs.conc, x_est)

    # Save this single experiment's settings and stats
    self.record_single_expt(tp, fp, fn, uncon_negs, determined, overdetermined, surep,
      unsurep, wrongly_undetected, score, rmse, cs.conc, (cs.conc > 0).astype(int), y, bool_y, x_est,
      infected, infected_dd)

    sys.stdout.write('\riter = %d / %d rmse: %.2f tp = %d fp = %d fn = %d' %
        (i, num_expts, rmse, tp, fp, fn))
    sys.stdout.flush()
    
    self.add_stats(tp, fp, fn, uncon_negs, determined, overdetermined, surep,
        unsurep, wrongly_undetected, score, rmse)

  def add_stats(self, tp, fp, fn, uncon_negs, determined, overdetermined,
      surep, unsurep, wrongly_undetected, score, rmse):
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
    self.determined += determined
    self.overdetermined += overdetermined
    self.surep += surep
    self.unsurep += unsurep
    self.num_expts_2_stage += (unsurep > 0)
    self.wrongly_undetected += wrongly_undetected
    self.total_score += score
    self.total_rmse += rmse

    # Compute minimum precision
    if tp + fp != 0:
      precision = tp / (tp + fp)
    else:
      precision = 1.
    self.all_precisions.append(precision)
    if precision < self.min_precision:
      self.min_precision = precision

    # Compute per-expt recall
    if tp + fn != 0:
      recall = tp / (tp + fn)
    else:
      recall = 1.
    self.all_recalls.append(recall)

    # Compute per-expt false pos
    self.all_fps.append(fp)
    # Compute per-expt false neg
    self.all_fns.append(fn)

  def print_stats(self, num_expts, header=False):
    preds = self.total_tp + self.total_fp
    actual = self.total_tp + self.total_fn
    if preds != 0:
      precision = self.total_tp / preds
    else:
      if actual > 0:
        precision = 0.
      else:
        precision = 1.

    if actual == 0:
      recall = 1.
    else:
      recall = self.total_tp / actual

    if num_expts == self.no_fp:
      self.avg_fp = 0
    else:
      self.avg_fp = self.total_fp / num_expts
    if num_expts == self.no_fn:
      self.avg_fn = 0
    else:
      self.avg_fn = self.total_fn / num_expts

    if header:
      print('******', self.name, 'Statistics', '******')
    print('No errors in %d / %d cases' % (self.no_error, num_expts))
    print('No fp in %d / %d cases' % (self.no_fp, num_expts))
    print('No fn in %d / %d cases' % (self.no_fn, num_expts))
    print('precision = %.6f, recall = %.6f' % (precision, recall))
    print('total tp =', self.total_tp, 'total fp = ', self.total_fp,
        'total fn = ', self.total_fn, 'unconfident negs = ',
        self.uncon_negs, 'determined = ', self.determined,\
            'overdetermined = ', self.overdetermined)
    print('avg fp = %.3f' % self.avg_fp, 'avg fn = %.3f' % self.avg_fn)

  def return_stats(self, num_expts):
    #precision = self.total_tp / float(self.total_tp + self.total_fp)
    #recall = self.total_tp / float(self.total_tp + self.total_fn)
    preds = self.total_tp + self.total_fp
    actual = self.total_tp + self.total_fn
    if preds != 0:
      precision = self.total_tp / preds
    else:
      if actual > 0:
        precision = 0.
      else:
        precision = 1.
    if actual == 0:
      recall = 1.
    else:
      recall = self.total_tp / actual
    stats = {
        'precision'     : precision,
        'recall'        : recall,
        'total_tp'      : self.total_tp,
        'total_fp'      : self.total_fp,
        'total_fn'      : self.total_fn,
        'no_error'      : self.no_error,
        'expts'         : num_expts,
        }

    if num_expts == self.no_fp:
      self.avg_fp = 0
    else:
      self.avg_fp = self.total_fp / num_expts
    if num_expts == self.no_fn:
      self.avg_fn = 0
    else:
      self.avg_fn = self.total_fn / num_expts

    self.preds = preds
    self.actual = actual
    self.precision = precision
    self.recall = recall
    self.expts = num_expts
    self.avg_unsurep = self.unsurep / num_expts
    self.avg_surep = self.surep / num_expts
    # Specificity is True Negative Rate
    self.specificity = specificity(self.precision, self.recall, self.n, self.d, preds)
    self.avg_score = self.total_score / self.expts
    self.rmse = self.total_rmse / self.expts
    return stats

def do_many_expts(n, d, t, num_expts=1000, xs=None, M=None,
    cross_validation=False,
    add_noise=False,
    algo='OMP',
    noise_magnitude=None,
    mr=None,
    mlabel="dummy_matrix_label"):

  # Test width. Max number of parallel tests available.
  #t = 12

  # Group size
  #n = 32

  # Number of infections. Sparsity
  #d = 1

  # Test assignment probability. Probability that a person gets assigned to a
  # test
  if d > 1:
    s = 1 / d
  else:
    s = 0.5
  #s = 0.5

  # lambda for regularization
  l = 300

  if xs is not None:
    assert num_expts == len(xs)

  if isinstance(algo, list):
    #print('list')
    expts = [CSExpts(item, n, d, t, mr, item, mlabel) for item in algo]
    ret_list = True
  else:
    #print('str')
    expts = [CSExpts(algo, n, d, t, mr, algo, mlabel)]
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
        try:
          expt.do_single_expt(i, num_expts, cs, cs.conc,
              algo=alg,
              cross_validation=cross_validation,
              add_noise=add_noise,
              noise_magnitude=noise_magnitude)
        except:
          print('Algo used:', alg)
          sys.stdout.flush()
          raise

    print('')
    for expt, alg in zip(expts, algo):
      try:
        expt.print_stats(num_expts)
      except:
        print('Algo used:', alg)
        sys.stdout.flush()
        raise

  except KeyboardInterrupt:
    for expt in expts:
      expt.print_stats(i)
    raise

  stats = []
  for expt in expts:
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

# This doesn't work since if p is very small you always find a col
# with no 1's
def get_small_random_matrix_bernoulli(t, n, p):
  zero = True
  while zero:
    zeros = 0
    matrix = np.random.binomial(1, p, size=(t, n))
    colsums = np.sum(matrix, axis = 0)
    for i, val in enumerate(colsums):
      if val == 0:
        print('col %d is 0' %i )
        zeros += 1
    if zeros == 0:
      zero = False
  return matrix

def get_small_random_matrix(t, n, col_sparsity):
  matrix = np.zeros((t, n))
  for col in range(n):
    ones = np.random.choice(range(t), size=col_sparsity, replace=False)
    matrix[ones, col] = 1
  return matrix


# Runs many parallel experiments and save stats
def run_many_parallel_expts_many_matrices(mats, mlabels, d_ranges, algos,
    num_expts, save):
  #if save:
  #  #pm = PickleManager(config.stats_pickle, config.stats_pickle_tmp)
  #  # stats is a 3-deep dictionary
  #  # stats[matrix][algo][d] points to list of 1000 experiments
  #  #stats = pm.get_stats_dict()
  #  #stats = stats_manager.load()
  all_exps_list = []
  for M, label, d_range in zip(mats, mlabels, d_ranges):
    #if save:
    #  if not label in stats:
    #    stats[label] = {}
    n = M.shape[1]
    t = M.shape[0]
    add_noise = True
    matrix = M
    n_jobs = 1
    explist = run_many_parallel_expts_internal(num_expts, n, t, add_noise, matrix, algos,
        d_range, n_jobs, xslist=[None for d in d_range], mlabel=label)
    all_exps_list.append(explist)
    if save:
      for algo, expts in zip(algos, explist):
        #if not algo in stats[label]:
        #  stats[label][algo] = {}
        for d, expt in zip(d_range, expts):
          #stats[label][algo][d] = [item.__dict__ for item in expt.single_expts]
          stats_manager.save(label, algo, d, [item.__dict__ for item in
            expt.single_expts])

      # Now that this matrix is done, we want to save the stats
      # We do this after every matrix so that even if the entire process is
      # cancelled, stats are still saved
      #pm.carefully_save_stats(stats)

  for M, label, explist in zip(mats, mlabels, all_exps_list):
    n = M.shape[1]
    t = M.shape[0]
    print(f"\nt = {t}, n = {n}, matrix = {label}\n")
    for algo, expts in zip(algos, explist):
      print('\n' + algo + f', num_expts = {num_expts}\n')
      print_expts(expts, num_expts, t)
  #return stats


def run_many_parallel_expts():
  #from experimental_data_manager import parse_israel_matrix
  #optimized_M_48_384_israel = parse_israel_matrix()

  num_expts = 100
  t = 45
  n = 105
  add_noise = True
  matrix = optimized_M_45_105_kirkman

  #t = 45
  #n = t * (t - 1) // 6
  #matrix = sts.sts(t)
  #matrix = matrix[:48, :384]
  #t = 48
  #n = 384

  algos = []
  algos.extend(['COMP'])
  #algos.extend(['SCOMP'])
  algos.extend(['SBL'])
  #algos.append('NNOMP')
  #algos.append('combined_COMP_NNOMP')
  #algos.append('NNOMP_random_cv')
  algos.extend(['combined_COMP_NNOMP_random_cv'])
  #algos.append('combined_COMP_SBL')
  #algos.append('l1ls_cv')
  algos.append('combined_COMP_l1ls_cv')
  #algos.append('combined_COMP_l1ls')
  d_range = [5, 8, 12, 15, 17] #list(range(1, 4))
  #d_range = list(range(10, 101, 10))
  #d_range = list(range(10, 101, 10))
  #d_range = [1]
  #d_range.extend([15, 20, 25, 30])
  n_jobs = 1

  run_many_parallel_expts_internal(num_expts, n, t, add_noise, matrix, algos,
      d_range, n_jobs, xslist=[None for d in d_range],
      mlabel="dummy_matrix_label")

# Separate out this function from above so that we can call on many matrices
def run_many_parallel_expts_internal(num_expts, n, t, add_noise, matrix,
    algos, d_range, n_jobs, xslist, mlabel):
  retvals = Parallel(n_jobs=n_jobs, backend='loky')\
  (\
      delayed(do_many_expts)\
      (
        n, d, t, num_expts=num_expts, M=matrix,\
        add_noise=add_noise,algo=algos, mr=None, xs=xs, mlabel=mlabel \
      )\
      for d, xs in zip(d_range, xslist)\
  )

  l = len(algos)
  statlist = [[] for i in range(l)]
  explist = [[] for i in range(l)]
  idx = range(l)
  for item in retvals:
    expts, stats = item
    for i, expt, stat in zip(idx, expts, stats):
      statlist[i].append(stat)
      explist[i].append(expt)

  print('\nn = %d, t = %d\n' % (expt.n, expt.t))
  #for i, algo in enumerate(algos):
  #  print('\n' + algo + '\n')
  #  print('\td\tPrecision\tRecall\ttotal_tests\n')
  #  for stat in statlist[i]:
  #    total_tests = t + stat['d'] / stat['precision']
  #    print('\t%d\t%.3f\t\t%.3f\t%.1f' % (stat['d'], stat['precision'],
  #      stat['recall'], total_tests))

  for i, algo in enumerate(algos):
    print('\n' + algo + '\n')
    #print('\td\tPrecision\tRecall\ttotal_tests\tnum_determined\tnum_overdetermined\n')
    #print('\td\tPrecision\tRecall\tsurep\tunsurep  avg_tests  2_stage\tWrongly_undetected')
    #print('\td\tPrecision\tRecall (Sensitivity) \tSpecificity\tsurep\tunsurep  avg_tests  2_stage')
    print_expts(explist[i], num_expts, t)
  return explist

def print_expts(expts, num_expts, t):
  #print('\td\tPrecision\tRecall (Sensitivity)\tSpecificity\tsurep\tunsurep\tfp\tfn\tRMSE\tmin_precision')
  print('\td\tPrecision\tSensitivity\tSpecificity\tsure\tundetermined\tfalse_pos\tfalse_neg\ttotal_tests')
  for expt in expts:
    if expt.precision == 0:
      total_tests = t + expt.n
    else:
      total_tests = t + expt.d / expt.precision
    #print('\t%d\t%.3f\t\t\t%.3f\t\t%.3f\t\t%4.1f\t%.1f\t%.1f\t%.1f\t%.2f\t%.3f' % (expt.d, expt.precision,
    #  expt.recall, expt.specificity, expt.avg_surep, expt.avg_unsurep,
    #  expt.total_fp / num_expts, expt.total_fn / num_expts, expt.rmse, expt.min_precision))
    print('\t%d\t%.3f\t\t%.3f\t\t%.3f\t\t%4.1f\t%.1f\t\t%.1f\t\t%.1f\t\t%d' % (expt.d, expt.precision,
      expt.recall, expt.specificity, expt.avg_surep, expt.avg_unsurep,
      expt.total_fp / num_expts, expt.total_fn / num_expts, expt.t + expt.avg_unsurep))
    #print('\nAll Precisions:', stats_helper.get_sorted_indices(expt.all_precisions))
    #print('\nAll Recalls:', stats_helper.get_sorted_indices(expt.all_recalls))
    #print('\nAll fps:', stats_helper.get_sorted_indices(expt.all_fps, invert=True))
    #print('\nAll fns:', stats_helper.get_sorted_indices(expt.all_fns, invert=True))
    #print('\nfns: ', np.unique(expt.all_fns, return_counts=True), '\n')
    #lower = expt.all_precisions[5]
    #upper = expt.all_precisions[-1]
    #print('95% Confidence Interval: ', (lower, upper))
    #print('\t%d\t%.3f\t\t%.3f\t%.1f\t\t%3d\t\t%3d' % (expt.d, expt.precision,
    #  expt.recall, total_tests, expt.determined, expt.overdetermined))
    #print('\t%d\t%.3f\t\t\t%.3f\t\t%.3f\t\t%4.1f\t%5.1f\t%7.1f\t%8d\t' % (expt.d, expt.precision,
    #  expt.recall, expt.specificity, expt.avg_surep, expt.avg_unsurep, expt.t +
    #  expt.avg_unsurep, expt.num_expts_2_stage, ))

def compare_sts_vs_kirkman():
  t = 27
  M = optimized_M_27_117_kirkman

  best_ds_kirkman = [run_with_matrix_n(M, t, n) for n in [50, 60, 70, 80, 90, 100, 110,
    117]]
  
  ds_kirkman = [item[0] for item in best_ds_kirkman]
  sps_kirkman = np.array([item[1] for item in best_ds_kirkman])
  pr_kirkman = np.array([item[2] for item in best_ds_kirkman])

  M = sts.sts(27)
  best_ds_sts = [run_with_matrix_n(M, t, n) for n in [50, 60, 70, 80, 90, 100, 110,
    117]]
  ds_sts = [item[0] for item in best_ds_sts]
  sps_sts = np.array([item[1] for item in best_ds_sts])
  pr_sts = np.array([item[2] for item in best_ds_sts])


  print('Kirkman:', ds_kirkman)
  print('STS:\t', ds_sts)
  print('Kirkman:', sps_kirkman)
  print('STS:\t', sps_sts)
  print('Kirkman:', pr_kirkman)
  print('STS:\t', pr_sts)

def compare_different_ns():
  explist = []
  t = 192
  #ns = list(range(400, 1000, 200)) + list(range(1000, 5000, 500)) + [5120]
  ns = list(range(400, 1000, 200))
  M = optimized_M_192_5120_social_golfer
  num_expts = 50
  for n in ns:
    expts = run_with_matrix_n(M, t, n, True, num_expts)
    explist.append(expts)
  for expts, n in zip(explist, ns):
    print(f'n = {n}, t = {t}')
    print_expts(expts, num_expts, t)

def compare_different_mats(mat_list, mat_labels):
  t = mat_list[0].shape[0]
  n = mat_list[0].shape[1]
  num_expts = 1000
  explist = []

  d_range = list(range(1, 16))
  xslist = []
  for d in d_range:
    xslist.append([create_infection_array_with_num_cases(n, d) for i in
      range(num_expts)])

  for M in mat_list:
    expts = run_with_matrix_n(M, t, n, True, num_expts, d_range, xslist)
    explist.append(expts)

  for expts, M, label in zip(explist, mat_list, mat_labels):
    t = M.shape[0]
    n = M.shape[1]
    print(f'n = {n}, t = {t}, matrix = {label}')
    print_expts(expts, num_expts, t)


def run_with_matrix_n(M, t, n, ret_explist=False, num_expts=1, d_range=None,
    xslist=None):
  assert n <= M.shape[1]
  assert t == M.shape[0]

  M = M[:, :n]
  add_noise = True

  algos = ['COMP']
  #d_range = list(range(5, 16)) + list(range(20, 41, 5))
  if not d_range:
    #d_range = list(range(1, 16)) + [20, 25, 30, 35]
    d_range = [5, 10, 15, 20, 25, 30, 35]
    assert not xslist
    xslist = [None for d in d_range]
  n_jobs = 4

  explist = run_many_parallel_expts_internal(num_expts, n, t, add_noise, M,
      algos, d_range, n_jobs, xslist, "dummy_label")
  expts = explist[0]
  sp = [expt.specificity for expt in expts]
  pr = [expt.precision for expt in expts]
  best_d = 0
  best_sp = 0
  best_pr = 0
  for i, d in enumerate(d_range):
    if sp[i] >= 0.945:
      best_d = d
      best_sp = sp[i]
      best_pr = pr[i]
  if ret_explist:
    return expts
  else:
    return best_d, best_sp, best_pr

def run_many_parallel_expts_mr():
  num_expts = 100
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

def generate_expts_deployed_matrices(only_these_labels=None, save=True):
  #labels = ['optimized_M_1', 'optimized_M_2']
  from get_test_results import MSizeToLabelDict
  tups = MSizeToLabelDict.values()
  all_labels = [tup[0] for tup in tups]
  if not only_these_labels:
    labels = all_labels
  else:
    for label in only_these_labels:
      assert label in all_labels
    labels = only_these_labels


  mats = [MDict[label] for label in labels]

  #d_ranges = [list(range(1, 5)) for i in range(len(mats))]
  d_ranges = [ [ tup[1] ] for tup in tups if tup[1] in labels]

  num_expts = 1000
  #algos = ['COMP', 'SBL', 'combined_COMP_NNOMP_random_cv',]
  algos = ['COMP', 'SBL', ]
  #algos = ['combined_COMP_l1ls_cv']
  run_many_parallel_expts_many_matrices(mats, labels, d_ranges, algos,
      num_expts, save)

def run_stats_for_these_matrices(labels, save):
  mats = [MDict[label] for label in labels]
  #d_ranges = [ list(range(1, 16)) + [20, 25, 30, 35, 40] for item  in labels]
  ts = [M.shape[0] for M in mats]
  d_ranges = [ [ 30 ] for t in ts] #list(range(1, 4))
  #d_ranges = [[5, 8, 10, 12, 15, 17, 20] for t in ts] #list(range(1, 4))
  #d_ranges = [[ 15 ] for t in ts] #list(range(1, 4))
  #d_ranges = [ list(range(1, (t // 3) + 1)) for t in ts ] 
  #d_ranges = [list(range(1, 6)) for label in labels]

  num_expts = 20
  #algos = ['COMP', 'SBL', 'combined_COMP_NNOMP_random_cv',
  #    'combined_COMP_l1ls_cv']
  algos = ['COMP', 'combined_COMP_SBL', 'combined_COMP_NNOMP_random_cv']
  run_many_parallel_expts_many_matrices(mats, labels, d_ranges, algos,
      num_expts, save)


