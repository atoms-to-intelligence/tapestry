from matrices import *
import sys

def convert_printable_format(M, f):
  n = M.shape[1]
  t = M.shape[0]
  #f.write("\n\nSample\t ------->\tTest\n\n")
  f.write("\n\nSample, \t \tTest\n\n")

  for i in range(n):
    sample = i + 1
    #f.write("%2d\t, ------->\t" % sample)
    f.write("%2d, " % sample)
    prev = None
    for j in range(t):
      test = j + 1
      if M[j,i] == 1:
        if prev is not None:
          f.write("%2d,\t" % prev)
        prev = test
    if prev is not None:
      f.write("%2d" % prev)

    f.write("\n")

if __name__ == '__main__':
  M = optimized_M_46_96_1
  convert_printable_format(M, sys.stdout)

