from cs import *

n = 40
t = 16
M = optimized_M_1

def load_cycle_times(filename):
  cts = []
  with open(filename, 'r') as f:
    for line in f:
      word = line.strip()
      val = float(word)
      cts.append(val)
  return np.array(cts)

def print_infected_people(bool_y):

  # unused params
  arr = np.zeros(n)
  mr = None
  d = 1
  s = 0.5
  l = 0.1

  cs = CS(n, t, s, d, l, arr, M, mr)
  #infected, infected_dd, score, tp, fp, fn, surep, unsurep = \
  #    cs.decode_comp_new(bool_y, compute_stats=False)

  infected = cs.decode_comp_new1(bool_y)
  infected_dd = np.zeros(n)
  print(infected)
  #print(infected_dd)

  sure_list = []
  unsure_list = []
  neg_list = []
  for i in range(n):
    if infected_dd[i] == 1:
      sure_list.append(i)
    elif infected[i] == 1:
      unsure_list.append(i)
    else:
      neg_list.append(i)

  print('Surely infected: ', sure_list)
  print('Possibly infected: ', unsure_list)
  print('Not infected: ', neg_list)

convert_matlab_format(M, sys.stdout)
filenames = ['test1.txt', 'test2.txt', 'test3.txt', 'test4.txt', 'test5.txt']
for i, filename in enumerate(filenames):
  print('Test %d:\n' % (i+1))
  cts = load_cycle_times(filename)
  bool_y = (cts < 33).astype(np.int32)
  print('Cycle times: ', cts)
  print('y: ', bool_y)
  print_infected_people(bool_y)
  print('\n')

