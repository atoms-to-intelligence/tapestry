# Combinatorial Orthogonal Matching Pursuit. Non-adaptive group-testing
# algorithm
#
# This algorithm only has false positives and no false negatives, provided
# each person is in at least one test.

import numpy as np

def create_infection_array_with_prob(n, p):
  arr = np.random.binomial(1, p, n)
  while sum(arr) == 0:
    arr = np.random.binomial(1, p, n)
  return arr

def create_infection_array_with_num_cases(n, d):
  arr = np.zeros(n)
  idx = np.random.choice(range(n), d, replace=False)
  arr[idx] = 1
  assert sum(arr) == d
  return arr

# Implements the COMP algorithm
class COMP:
  def __init__(self, n, t, s, arr):
    self.n = n
    self.t = t
    self.s = s
    self.arr = arr
    self.initialize_M()

  # Test matrix
  def initialize_M(self):
    self.M = np.random.binomial(1, s, size=(self.n, self.t))
    for person in range(n):
      while sum(self.M[person]) == 0:
        self.M[person] = np.random.binomial(1, self.s, size=(self.t))
        print('changing tests for person', person)
    # Count of number of people in only 1 test
    print(sum([item == 1 for item in [sum(row) for row in self.M]]))

  # Get test results from infection array and test matrix
  # In reality this will happen from group testing
  def get_results(self):
    # Now get results from infections
    positive_tests = []
    for test in range(self.t):
      positive = 0
      for person in range(self.n):
        positive += (self.arr[person] * self.M[person, test])
      if positive > 0:
        positive_tests.append(test)

    #positive_tests = np.array(positive_tests)
    return positive_tests

  def decode(self, positive_tests):
    errors = 0
    total = 0
    if positive_tests:
      positive_people = np.sum(self.M[:, positive_tests], axis=-1)
      for person in range(n):
        if positive_people[person] == sum(self.M[person]):
          total += 1
          print(person)
          if arr[person] == 0:
            #print('False positive')
            errors += 1
          #else:
            #print('True positive')
    assert (total - errors) == sum(self.arr)
    return errors, total

# Test width. Max number of parallel tests available.
t = 500

# Infection probability
p = 0.001

# Test failure probability
q = 0.

# Group size
n = 10000

# Number of infections
d = 200

# Test assignment probability. Probability that a person gets assigned to a
# test
s = 1. / 20

arr = create_infection_array_with_num_cases(n, d)
comp = COMP(n, t, s, arr)
results = comp.get_results()
errors, total = comp.decode(results)
print('errors / total =', errors, '/', total)

