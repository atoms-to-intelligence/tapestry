# Combinatorial Orthogonal Matching Pursuit. Non-adaptive group-testing
# algorithm
#
# This algorithm only has false positives and no false negatives, provided
# each person is in at least one test.

import numpy as np
import sys

def create_infection_array_with_prob(n, p):
  arr = np.random.binomial(1, p, n)
  #while sum(arr) == 0:
  #  arr = np.random.binomial(1, p, n)
  return arr

def create_infection_array_with_num_cases(n, d):
  arr = np.zeros(n)
  idx = np.random.choice(range(n), d, replace=False)
  arr[idx] = 1
  assert sum(arr) == d
  return arr

# Implements the COMP algorithm
class COMP:
  def __init__(self, n, t, s, d, arr):
    self.n = n
    self.t = t
    self.s = s
    self.d = d
    self.arr = arr
    self.initialize_M()

  # Test matrix
  def initialize_M(self):
    self.M = np.random.binomial(1, self.s, size=(self.n, self.t))

  # Get test results from infection array and test matrix
  # In reality this will happen from group testing
  def get_results(self):
    arr_T = np.expand_dims(self.arr, axis=-1)
    infections = np.sum(arr_T * self.M, axis=0)
    infections = (infections>0).astype(np.int32)
    return infections

  def decode(self, infections):
    infected = np.all(self.M * infections == self.M, axis=1).astype(np.int32)
    total = sum(infected)
    errors = sum(infected != self.arr)
    #print(errors, total)
    assert (total - errors) == sum(self.arr)
    return errors, total

  def decode_comp_new1(self, infections, compute_stats=True):
    n = self.n
    t = self.t
    infected = np.zeros(n)
    for person in range(n):
      is_infected = True
      for test in range(t):
        if self.M[person, test] == 1 and infections[test] == 0:
          is_infected = False

      if is_infected:
        infected[person] = 1

    return infected



  def decode_comp_new(self, infections, compute_stats=True):
    infected = np.all(self.M * infections == self.M, axis=1).astype(np.int32)
    infected_dd = self.get_infected_dd(infected, infections)
    assert np.all(infected_dd - infected <= 0)

    num_infected_in_test = np.zeros(self.t, dtype=np.int32)
    for test in range(self.t):
      for person in range(self.n):
        if infected[person] and self.M[person, test] == 1:
          num_infected_in_test[test] += 1
    
    if compute_stats:
      tpos = (infected * self.arr)
      fneg = (1 - infected) * self.arr
      fpos = infected * (1 - self.arr)

      tp = sum(tpos)
      fp = sum(fpos)
      fn = sum(fneg)
      
      surep = np.sum(infected_dd)
      unsurep = np.sum(infected * (1 - infected_dd))

      assert surep + unsurep == tp + fp
      assert surep <= tp
      assert unsurep >= fp
    else:
      tp = 0
      fp = 0
      fn = 0
      surep = 0
      unsurep = 0

    return infected, infected_dd, 0., tp, fp, fn, surep, unsurep,\
        num_infected_in_test

  # Ax = y
  # x -> infected
  # y -> infections/results
  # A -> matrix
  def get_infected_dd(self, x, y):
    assert len(x) == self.n
    assert len(y) == self.t
    #non_zero_cols,  = np.nonzero(x)
    #non_zero_rows,  = np.nonzero(y)

    A = self.M.T
    #A = np.take(A, non_zero_cols, axis=1)
    #A = np.take(A, non_zero_rows, axis=0)
    
    dd = np.zeros(self.n)
    for row in A:
      row = row * x
      if sum(row) == 1:
        non_zero_col,  = np.nonzero(row)
        dd[non_zero_col[0]] = 1
    return dd

if __name__ == '__main__':
  # Test width. Max number of parallel tests available.
  t = 96

  # Infection probability
  p = 0.001

  # Test failure probability
  q = 0.

  # Group size
  n = 1000

  # Number of infections
  d = 40

  # Test assignment probability. Probability that a person gets assigned to a
  # test
  s = 1. / d

  count = 0
  extra_tests = 0
  false_positives = 0
  num_expts = 10
  for i in range(num_expts):
    #sys.stdout.write('\r')
    sys.stdout.write('\n')
    sys.stdout.write('Iteration %d / %d' % (i+1, num_expts))
    arr = create_infection_array_with_num_cases(n, d)
    #arr = create_infection_array_with_prob(n, p)
    comp = COMP(n, t, s, d, arr)
    results = comp.get_results()
    errors, total = comp.decode(results)
    sys.stdout.write(', errors / total = ' + str(errors) + ' / ' + str(total))
    if errors > 0:
      count += 1
      extra_tests += total
      false_positives += errors
  sys.stdout.write('\n')
  print('Number of cases in which error was made:', count, '/', num_expts)
  print('Number of extra tests done:', extra_tests)
  print('Number of false positives:', false_positives)


