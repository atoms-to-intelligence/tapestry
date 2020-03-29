from cs import *
import pickle

# Out of a list of matrices, choose the best one
#
# Generate 1000 x's to do the expt. The same x will be used to evaluate all
# the matrices
# list of stats will be returned for each matrix
def find_best_matrix(n, d, t, Ms, num_expts=1000):
  xs = [create_infection_array_with_num_cases(n, d) for i in range(num_expts)]
  stats = []
  for M in Ms:
    item = do_many_expts(n, d, t, num_expts, xs, M,
        cross_validation=False,
        add_noise=True)
    #print(item)
    stats.append(item)

  return stats

n = 40
d = 2
t = 16
Ms = [np.random.binomial(1, 0.5, size=(t, n)) for i in range(100)]
Ms.append(optimized_M)
stats = find_best_matrix(n, d, t, Ms, num_expts=1000)
for stat, M in zip(stats, Ms):
  if stat['precision'] >= 0.80 and stat['recall'] >= 0.90:
    print(stat)

print(stat)