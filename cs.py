# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
import numpy as np
import math
from sklearn.linear_model import Lasso, LassoLars, LassoCV, LassoLarsCV
import pylops
from joblib import Parallel, delayed

from comp import create_infection_array_with_num_cases, COMP
from matrices import *
import nnompcv
import sbl
import l1ls

import config

# Numpy configuration
np.set_printoptions(precision=3)
# Numpy should raise Exception on division by zero so that we can catch programming errors
np.seterr(all='raise')

# Use compressed sensing to solve 0.5*||Mx - y||^2 + l * ||x||_1
class CS(COMP):
  def __init__(self, n, t, s, d, l, arr, M=None, mr=None):
    super().__init__(n, t, s, d, arr)
    if M is not None:
      assert n == M.shape[1]
      assert t == M.shape[0]
      self.M = M.T
      #print(self.M.shape)
    self.create_conc_matrix_from_infection_array(arr)
    self.l = l
    self.mr = mr

  # Multiply actual conc matrix to M. This is the part done by mixing samples
  # and qpcr
  def get_quantitative_results(self, conc, add_noise=False,
      noise_magnitude=None, noise_model='exponential_gaussian'):
    conc = np.expand_dims(conc, axis=-1)
    #print(self.M.shape, conc.shape)
    y = np.matmul(self.M.T, conc).flatten()
    # print('lol')
    sigval = 0.
    if add_noise:
      if noise_magnitude is not None:
        # This error is independent of the magnitude
        sigval = noise_magnitude
        error = np.random.normal(0., sigval)
        y = y + error
        raise ValueError("This noise model is incorrect and hence disabled. "
            "Enable this if you know what you're doing")
      elif config.noise_model == 'variable_gaussian':
        # This error is proportional to magnitude of y
        sigval = 0.01*np.absolute(y)
        error = np.random.normal(0., sigval)
        y = y + error
      elif config.noise_model == 'exponential_gaussian':
        # This noise accounts for cycle time variability
        # Cycle time is assumed to be Gaussian distributed, due to which log
        # of y is Gaussian. Hence 
        p = 0.95
        error = np.random.normal(0., 0.1, size=self.t)
        #print('Original y', y)
        #print('error exponents', error)
        y = y * ((1+p) ** error)
        #print('p:', p, 'y with exponentiated error:', y)
      else:
        raise ValueError('Invalid noise model %s' % noise_model)


      if config.bit_flip_prob > 0 :
        raise ValueError('This is probably a mistake')
        print('before flip y = ', y)
        mask = (y > 0).astype(np.int32)
        flip = np.random.binomial(1, config.bit_flip_prob, self.t)
        flip = flip * mask
        y = y * (1 - flip)
        print('after flip y = ', y)

    return y, sigval

  # Initial concentration of RNA in each sample
  def create_conc_matrix_from_infection_array(self, arr):
    # Fix tau to 0.01 * minimum value we expect in x
    self.tau = 0.01 * 1 / 32768.
    #self.tau = 0.01 * 0.1
    #conc = 1 + np.random.poisson(lam=5, size=self.n)
    conc = np.random.randint(low=1, high=32769, size=self.n) / 32768
    #conc = 0.1 + 0.9 * np.random.rand(self.n)
    #conc = np.random.randint(low=1, high=11, size=self.n) / 10.
    #conc = np.ones(self.n)
    self.conc = conc * arr # Only keep those entries which are 

  # Solve the CS problem using Lasso
  #
  # y is results
  def decode_lasso(self, results, algo='lasso', prefer_recall=False,
      compute_stats=True):
    determined = 0
    overdetermined = 0
    # Add check if system is determined or overdetermined
    if self.t == self.n:
      determined = 1
    elif self.t > self.n:
      overdetermined = 1

    prob1 = None
    prob0 = None
    if algo == 'lasso':
      #lasso = LassoLars(alpha=self.l)
      lasso = Lasso(alpha=self.l, max_iter=10000)
      #lasso = LassoCV(n_alphas=100)
      lasso.fit(self.M.T, results)
      #print('best lambda = ', lasso.alpha_)
      answer = lasso.coef_
    elif algo == 'OMP':
      temp_mat=(self.M.T).astype(float)
      temp_mat=pylops.MatrixMult(temp_mat)
      answer = pylops.optimization.sparsity.OMP(temp_mat, results, 10000,
          sigma=0.001)[0]
    elif algo== 'NNOMP':
      # Max d that can be detected by NNOMP is equal to number of rows
      answer=nnompcv.nnomp(self.M.T.astype('float'),0,results,0, self.t, cv=False)
    elif algo=='NNOMPCV':
      temp_mat = (self.M.T).astype(float)
      mr = math.ceil(0.9*temp_mat.shape[1])
      m = temp_mat.shape[1]
      Ar = temp_mat[0:mr, :]
      Acv = temp_mat[mr+1:m, :]
      yr = results[0:mr]
      ycv = results[mr+1:m]
      #print('yo')
      # Max d that can be detected by NNOMP is equal to number of rows
      answer = nnompcv.nnomp(Ar, Acv, yr, ycv, self.t, cv=True)
    elif algo == 'NNOMP_loo_cv':
      answer, prob1, prob0 = self.decode_nnomp_multi_split_cv(results, 'loo_splits')
    elif algo == 'NNOMP_random_cv':
      # Skip cross-validation for really small cases
      if np.sum(results) == 0:
        answer = np.zeros(self.n)
      elif self.t < 4:

        # Max d that can be detected by NNOMP is equal to number of rows
        answer = nnompcv.nnomp(self.M.T.astype('float'),0,results,0, self.t, cv=False)
      else:
        answer, prob1, prob0 = self.decode_nnomp_multi_split_cv(results, 'random_splits')
    elif algo.startswith('combined_COMP_'):
      #print('Doing ', algo)
      l = len('combined_COMP_')
      secondary_algo = algo[l:]
      answer, infected, prob1, prob0, determined, overdetermined =\
          self.decode_comp_combined(results, secondary_algo,
          compute_stats=compute_stats)
    elif algo == 'SBL':
      A = self.M.T
      y = results
      #assert np.all(y!=0)
      #A1 = A/y[:, None]
      #y1 = np.ones(y.shape)
      # This sigval computation was not correct
      #sigval = 0.01 * np.linalg.norm(y, 2)
      #sigval = 0.01 * np.mean(y1)
      sigval = 0.01 * np.mean(y)
      #answer = sbl.sbl(A1, y1, sigval, self.tau)
      answer = sbl.sbl(A, y, sigval, self.tau)
      #print(answer)
    elif algo == 'l1ls':
      A = self.M.T
      y = results
      if np.all(y == 0):
        answer = np.zeros(self.n)
      else:
        answer = l1ls.l1ls(A, y, self.l, self.tau)
    elif algo == 'l1ls_cv':
      A = self.M.T
      y = results
      sigval = 0.01 * np.mean(y)
      if np.all(y == 0):
        answer = np.zeros(self.n)
      else:
        answer = l1ls.l1ls_cv(A, y, sigval, self.tau)
    else:
      raise ValueError('No such algorithm %s' % algo)

    score = np.linalg.norm(answer - self.conc) / math.sqrt(self.t)
    infected = (answer != 0.).astype(np.int32)

    if prob1 is None:
      assert prob0 is None
      prob1 = np.array(infected)
      prob0 = np.array(1 - infected)
    
    num_unconfident_negatives = 0
    if prefer_recall:
      # Report the unconfident -ves as +ve
      negatives = (infected == 0).astype(np.int32)
      unconfident_negatives = negatives * (prob0 < 0.6).astype(np.int32)
      num_unconfident_negatives = np.sum(unconfident_negatives)
      infected = infected + unconfident_negatives

    # Get definite defects
    y = results
    bool_y = (y > 0).astype(np.int32)
    _infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, _unsurep, _ =\
        self.decode_comp_new(bool_y, compute_stats=compute_stats)

    #print(infected.shape)
    #print(infected_dd.shape)
    # Compare definite defects with ours to detect if our algorithm doesn't
    # detect something that should definitely have been detected
    wrongly_undetected = np.sum(infected_dd - infected_dd * infected)

    infected = (infected + infected_dd > 0).astype(np.int32)

    if compute_stats:
      # Compute stats
      tpos = (infected * self.arr)
      fneg = (1 - infected) * self.arr
      fpos = infected * (1 - self.arr)
      
      tp = sum(tpos)
      fp = sum(fpos)
      fn = sum(fneg)
      
      assert surep <= tp
      unsurep = tp + fp - surep
    else:
      tp = 0
      fp = 0
      fn = 0
      surep = 0
      unsurep = 0


    num_infected_in_test = np.zeros(self.t, dtype=np.int32)
    for test in range(self.t):
      for person in range(self.n):
        if infected[person] > 0 and self.M[person, test] == 1:
          num_infected_in_test[test] += 1

    return answer, infected, infected_dd, prob1, prob0, score, tp, fp, fn,\
        num_unconfident_negatives, determined, overdetermined, surep,\
        unsurep, wrongly_undetected, num_infected_in_test

  def decode_lasso_for_cv(self, train_Ms, train_ys, test_Ms, test_ys,
      algo='lasso', l=None, sigma=None):

    if algo == 'lasso' and l is None:
      raise ValueError('Need l for algo lasso')
    elif algo == 'OMP' and sigma is None:
      raise ValueError('Need sigma for algo OMP')

    scores = []
    for train_M, train_y, test_M, test_y in zip(train_Ms, train_ys, test_Ms,
        test_ys):
      #print('Doing lasso with')
      #print(train_M.shape, train_y.shape, test_M.shape, test_y.shape)
      if algo == 'lasso':
        lasso = Lasso(alpha=l, max_iter=10000)
        lasso.fit(train_M, train_y)
        pred_y = lasso.predict(test_M)
      elif algo == 'OMP':
        pass

      score = np.linalg.norm(test_y - pred_y) / len(test_y)
      scores.append(score)

    avg_score = np.average(scores)
    max_score = max(scores)
    median_score = np.median(scores)
    min_score = np.min(scores)
    return avg_score
    #return min_score

  # Get num random splits with given fraction. Sensing matrix will have at
  # most frac fraction of rows
  def return_random_splits(self, y, num, frac, mr=None):
    if mr is None:
      mr = math.floor(frac * self.t)
    else:
      assert mr < self.t
    r = self.t - mr
    # Following code only works for r > 1
    assert r > 1
    
    train_Ms = []
    test_Ms = []
    train_ys = []
    test_ys = []
    M = self.M.T # Uggh
    for i in range(num):
      perm = np.random.permutation(range(self.t))
      r_idx = perm[:r]
      m_idx = perm[r:]

      train_M = np.delete(M, r_idx, axis=0)
      train_y = np.delete(y, r_idx, axis=0)

      test_M = np.delete(M, m_idx, axis=0)
      test_y = np.delete(y, m_idx, axis=0)

      train_Ms.append(train_M)
      train_ys.append(train_y)
      test_Ms.append(test_M)
      test_ys.append(test_y)

    return train_Ms, train_ys, test_Ms, test_ys

  # Return splits for leave-one-out cross-validation
  def return_loo_cv_splits(self, y):
    train_Ms = []
    test_Ms = []
    train_ys = []
    test_ys = []

    # Unfortunately self.M is n x t so we need to transpose it
    M = self.M.T

    # Each row will be left out once as test_M
    for r in range(self.t):
      train_M = np.delete(M, r, axis=0)
      test_M = np.expand_dims(M[r], axis=0)
      train_y = np.delete(y, r, axis=0)
      test_y = np.array([y[r]])

      train_Ms.append(train_M)
      train_ys.append(train_y)
      test_Ms.append(test_M)
      test_ys.append(test_y)

    return train_Ms, train_ys, test_Ms, test_ys

  # Find best d by cross-validation using these splits
  #
  # Best d is the one found by majority of the splits
  def get_d_nnomp_cv(self, splits, max_d, resolve_method='voting', algo='NNOMP'):
    train_Ms, train_ys, test_Ms, test_ys = splits
    counts = np.zeros(max_d + 1)
    cum_error = np.zeros(max_d)
    # Keeps count of number of times each sample was declared as +ve
    x_ones = np.zeros(self.n)
    num_splits = len(train_Ms)
    for train_M, train_y, test_M, test_y in zip(train_Ms, train_ys, test_Ms,
        test_ys):
      x, error, d, errors = nnompcv.nnomp(train_M, test_M, train_y, test_y,
          max_d, cv=True)
      answer = (x > 0).astype(np.int32)
      x_ones += answer
      counts[d] += 1
      #print('Errors: ', np.array(errors))
      if errors:
        cum_error += errors

    best_d_maj = np.argmax(counts) + 1
    best_d_error = np.argmin(cum_error) + 1
    prob_of_one = x_ones / num_splits
    prob_of_zero = 1 - prob_of_one
    #print('prob of one:', prob_of_one)
    #print('prob of zero:', prob_of_zero)
    if resolve_method == 'voting':
      return best_d_maj, prob_of_one, prob_of_zero
    elif resolve_method == 'error':
      return best_d_error, prob_of_one, prob_of_zero
    else:
      raise ValueError('Invalid resolve method %s' % resolve_method)

  # Do leave one out splits
  # get best d from those splits
  # Run final nnomp algorithm using best d and entire matrix
  def decode_nnomp_multi_split_cv(self, y, method='random_splits'):
    if method == 'random_splits':
      splits = self.return_random_splits(y, 100, frac=0.7, mr=self.mr)
    elif method == 'loo_splits':
      splits = self.return_loo_cv_splits(y)

    best_d, prob1, prob0 = self.get_d_nnomp_cv(splits, max_d=self.t)
    if config.prefer_recall:
      best_d = 2 * best_d
    x = nnompcv.nnomp(self.M.T.astype('float'), 0, y, 0,
        best_d, cv=False)
    return x, prob1, prob0
    

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
      score = self.decode_lasso_for_cv(train_Ms, train_ys, test_Ms, test_ys,
          l=l)
      scores.append(score)
    scores = np.array(scores)
    idx = np.argmin(scores)
    #print(idx)
    self.l = ll[idx]
    #print('lambdas = ', ll)
    #print('scores = ', scores)
    print('Choosing lambda = %.4f' % self.l, 'score = %.4f' % score)
    return self.l

  # Filter out those entries of x which are definitely 0 using COMP.
  # Remove corresponding columns from M.
  def decode_comp_combined(self, y, secondary_algo, test=False,
      compute_stats=True):
    # This assertion is needed because mr depends on number of rows.
    # Since number of rows will change for the internal CS, use frac instead
    assert self.mr == None

    bool_y = (y > 0).astype(np.int32)
    infected_comp, infected_dd, _score, _tp, _fp, _fn, surep, unsurep, _ =\
        self.decode_comp_new(bool_y, compute_stats)

    # Find the indices of 1's above. These will be retained. Rest will be
    # discarded
    #print('Comp output: ', infected_comp)
    non_zero_cols,  = np.nonzero(infected_comp)
    non_zero_rows,  = np.nonzero(y)
    #print('Indices of Non-zero columns:', non_zero_cols)
    #print('Indices of Non-zero rows:', non_zero_rows)

    A = self.M.T
    A = np.take(A, non_zero_cols, axis=1)
    A = np.take(A, non_zero_rows, axis=0)
    #print('y: ', y)
    #print('Non-zero rows:', non_zero_rows)
    #print('Non-zero rows len:', non_zero_rows.shape)
    #print('Shape of remaining A:', A.shape)
    #print('Remaining A: ', A)

    y = y[non_zero_rows]
    #print('Remaining y:', y)
    
    x = self.conc
    x = x[non_zero_cols]
    #print('Remaining x:', x)
    arr = (x>0).astype(np.int32)
    # Now solve using this new A and y

    # parameters d and s do not matter. They are not used in the algorithm
    n = A.shape[1]
    t = A.shape[0]
    d = self.d
    s = self.s
    l = 0.1

    if t == 0:
      assert n == 0
    elif n == 0:
      assert t == 0

    infected = np.zeros(self.n)
    answer = np.zeros(self.n)
    prob1_new = np.zeros(self.n)
    prob0_new = np.ones(self.n)
    determined = 1
    overdetermined = 0
    # Calling internal algo is needed only when there is at least one infection
    if t != 0:
      # Create another CS class to run the secondary algorithm
      # Better to set mr parameter to None since it depends on number of rows
      # and will change for this internal CS object. frac will be used instead
      # for deciding splits
      _cs = CS(n, t, s, d, l, arr, A, mr=None)
      _cs.conc = x

      answer_internal, infected_internal, infected_dd, prob1, prob0, score, tp, fp, fn, _, determined,\
          overdetermined, surep, unsurep, wrongly_undetected, _ =\
          _cs.decode_lasso(y, secondary_algo, compute_stats=compute_stats)
      for ans, val, idx in zip(answer_internal, infected_internal, non_zero_cols):
        infected[idx] = val
        answer[idx] = ans

      for p1, p0, idx in zip(prob1, prob0, non_zero_cols):
        prob1_new[idx] = p1
        prob0_new[idx] = p0

    # tp, fp and fn will be correct for the internal algo
    if test:
      return infected, prob1_new, prob0_new, score, tp, fp, fn
    else:
      return answer, infected, prob1_new, prob0_new, determined, overdetermined

  def decode_qp(self, results):
    pass

  def print_matrix(self):
    pass

  def pickle_dump(self, filename):
    pass

if __name__ == '__main__':
  raise ValueError('Running experiments has been moved to cs_expts.py. '
      'Either use that or sel_matrix.py')

