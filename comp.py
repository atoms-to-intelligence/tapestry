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
    #for person in range(n):
    #  while sum(self.M[person]) == 0:
    #    self.M[person] = np.random.binomial(1, self.s, size=(self.t))
    #    print('changing tests for person', person)
    # Count of number of people in only 1 test
    #print(sum([item == 1 for item in [sum(row) for row in self.M]]))

  # Get test results from infection array and test matrix
  # In reality this will happen from group testing
  def get_results(self):
    arr_T = np.expand_dims(self.arr, axis=-1)
    infections = np.sum(arr_T * self.M, axis=0)
    infections = (infections>0).astype(np.int32)
    #print(infections)
    return infections
    #positive_tests = np.argwhere(infections > 0)
    #positive_tests =  positive_tests[:,0].tolist()
    #print(positive_tests)
    #return positive_tests
    # Now get results from infections
    #positive_tests = []
    #for test in range(self.t):
    #  positive = 0
    #  for person in range(self.n):
    #    positive += (self.arr[person] * self.M[person, test])
    #  if positive > 0:
    #    positive_tests.append(test)

    #positive_tests = np.array(positive_tests)
    #return positive_tests

  def decode(self, infections):
    infected = np.all(self.M * infections == self.M, axis=1).astype(np.int32)
    total = sum(infected)
    errors = sum(infected != self.arr)
    #print(errors, total)
    assert (total - errors) == sum(self.arr)
    return errors, total

if __name__ == '__main__':
  # Test width. Max number of parallel tests available.
  t = 384

  # Infection probability
  p = 0.001

  # Test failure probability
  q = 0.

  # Group size
  n = 1000

  # Number of infections
  d = 10

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


