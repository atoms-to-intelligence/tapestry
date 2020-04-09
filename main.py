from cs import *

n = 40
t = 16
M = optimized_M_16_40_ncbs
filenames = ['test1.txt', 'test2.txt', 'test3.txt', 'test4.txt', 'test5.txt']

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
  infected, infected_dd, score, tp, fp, fn, surep, unsurep,\
      num_infected_in_test = \
      cs.decode_comp_new(bool_y, compute_stats=False)

  for test in range(t):
    if bool_y[test] > 0 and num_infected_in_test[test] == 0:
      print('Test %d is infected but no infected people found' % test)
    if bool_y[test] == 0 and num_infected_in_test[test] > 0:
      print('Test %d is not infected but infected people found' % test)
    
  #infected = cs.decode_comp_new1(bool_y)
  #infected_dd = np.zeros(n)
  print(infected)
  print(infected_dd)

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

def print_results_COMP():
  #convert_matlab_format(M, sys.stdout)
  filenames = ['test1.txt', 'test2.txt', 'test3.txt', 'test4.txt', 'test5.txt']
  for i, filename in enumerate(filenames):
    print('Test %d:\n' % (i+1))
    cts = load_cycle_times(filename)
    #bool_y = (cts < 33.479).astype(np.int32)
    bool_y = (cts < 33).astype(np.int32)
    print('Cycle times: ', cts)
    print('y: ', bool_y)
    print_infected_people(bool_y)
    print('\n')

def find_max_positive_cycle_time():
  cutoff = 33
  CTS = []
  positive_cts = []
  for i, filename in enumerate(filenames):
    cts = load_cycle_times(filename)
    CTS.append(cts)
    bool_cts = (cts < cutoff).astype(np.int32)
    pos_cts = cts * bool_cts
    positive_cts.append(pos_cts)

  CTS_arr = np.concatenate(CTS)
  #print(CTS_arr.shape)
  #print(CTS_arr)
  bool_y = (CTS_arr < cutoff).astype(np.int32)
  #print(bool_y)
  CTS_arr = CTS_arr * bool_y
  m = np.max(CTS_arr)
  return m, CTS, positive_cts

m, cts_list , positive_cts= find_max_positive_cycle_time()

np.set_printoptions(linewidth=800)
p = 0.95
for cts in positive_cts:
  bool_y = (cts > 0).astype(np.int32)
  y = (1+p) ** (m - cts)
  y = y * bool_y
  print(y)


