# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
# Combinatorial Orthogonal Matching Pursuit. Non-adaptive group-testing
# algorithm
#
# This algorithm only has false positives and no false negatives, provided
# each person is in at least one test.

import numpy as np
import sys

def scomp(A,y):
  idx = -1
  zeroY = []
  nonZeroY = []
  samps = list(range(A.shape[1]))
  for yTest in y:
    idx = idx+1
    if yTest == 0:
      zeroY.append(idx)
    if yTest > 0:
      nonZeroY.append(idx)

  #yZero = find(~y); subA = A; colsRemove = [];
  removeCols = []
  for i in zeroY:
    row = A[i,:]
    idx = -1
    for cols in row:
      idx = idx+1
      if cols > 0:
        removeCols = list(set().union(removeCols, [idx]))

  nonZeroSamps = list(set(samps) - set(removeCols))
  #print(nonZeroSamps)
  subA = np.copy(A)
  #subA = np.delete(subA, removeCols, 1)
  subA[:,removeCols] = 0
  subA = np.delete(subA, zeroY, 0)
  #print(subA.shape)
  greedyCover = []
  
  surePos = []
  for i in range(subA.shape[0]):
    row = subA[i,:]
    if np.sum(row) == 1:
      idx = -1
      for j in row:
        idx = idx + 1
        if j > 0:
          break
      surePos = list(set().union(surePos,[idx]))
  
  removeRows = []
  for i in range(subA.shape[0]):
    row = subA[i,surePos]
    if np.sum(row) > 0:
      removeRows = list(set().union(removeRows,[i]))
  
  subA = np.delete(subA, removeRows, 0)
  greedyCover = surePos

  while subA.shape[0] > 0:
    idx = np.argmax(np.sum(subA, axis = 0))
    greedyCover = list(set().union(greedyCover, [idx]))
    removeRows = []

    for i in range(subA.shape[0]):
      row = subA[i,idx]
      if row > 0:
        #print(removeRows, i)
        removeRows = list(set().union(removeRows, [i]))
    
    subA = np.delete(subA, removeRows, 0)
  
  return greedyCover


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


  def compute_stats(self, infected, infected_dd, compute=True):
    if compute:
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
    return tp, fp, fn, surep, unsurep

  def compute_num_infected_in_test(self, infected):
    num_infected_in_test = np.zeros(self.t, dtype=np.int32)
    for test in range(self.t):
      for person in range(self.n):
        if infected[person] > 0 and self.M[person, test] == 1:
          num_infected_in_test[test] += 1
    return num_infected_in_test

  def decode_comp_new(self, infections, compute_stats=True):
    assert np.all(np.logical_or(infections == 0, infections == 1))
    infected = np.all(self.M * infections == self.M, axis=1).astype(np.int32)
    infected_dd = self.get_infected_dd(infected, infections)
    assert np.all(infected_dd - infected <= 0)

    num_infected_in_test = self.compute_num_infected_in_test(infected)

    tp, fp, fn, surep, unsurep = self.compute_stats(infected, infected_dd,
        compute=compute_stats)
    return infected, infected_dd, 0., tp, fp, fn, surep, unsurep,\
        num_infected_in_test

  # SCOMP: COMP ka baap
  #
  # 1. Run COMP and remove -ves
  # 2. Get definite defects
  # 3. After getting definite defects, declare those tests as negative and delete
  #    the definite defect columns.
  # 4. Goto step 1
  def decode_scomp(self, infections, compute_stats=True):
    A = self.M.T
    y = infections
    infected_idx = scomp(A, y)
    infected = np.zeros(self.n, dtype=np.int32)
    infected[infected_idx] = 1
    infected_comp = np.all(self.M * infections == self.M, axis=1).astype(np.int32)
    infected_dd = self.get_infected_dd(infected_comp, infections)
    assert np.all(infected_dd - infected <= 0)

    num_infected_in_test = self.compute_num_infected_in_test(infected)

    tp, fp, fn, surep, unsurep = self.compute_stats(infected, infected_dd,
        compute=compute_stats)
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


