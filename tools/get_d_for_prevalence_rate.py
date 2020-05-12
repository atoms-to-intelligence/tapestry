# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
from scipy.stats import binom

def get_d_for_n_p(n, p):
  b = binom(n, p)
  mean = b.mean()
  std =  b.std()
  max_d = mean + 2 * std
  print("n = %d, d_expect = %.1f, d_std = %.1f, d_max = %.1f" % (n, mean, std, max_d))

for n in [100, 300, 500, 1000]:
  get_d_for_n_p(n, 0.02)
