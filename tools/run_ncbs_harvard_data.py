# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
import sys
sys.path.append(".")

from core.cs import *
import pandas as pd

def load_cycle_times(filename):
  cts = []
  with open(filename, 'r') as f:
    for line in f:
      word = line.strip()
      val = float(word)
      cts.append(val)
  return np.array(cts)

def print_infected_people(y, bool_y, algo):
  print("---------------------")
  print(algo)
  print("---------------------")

  # unused params
  arr = np.zeros(n)
  mr = None
  d = 1
  s = 0.5
  l = 0.1

  cs = CS(n, t, s, d, l, arr, M, mr)
  if algo == 'COMP':
    infected, infected_dd, score, tp, fp, fn, surep, unsurep,\
        num_infected_in_test = \
        cs.decode_comp_new(bool_y, compute_stats=False)
    x = np.zeros(n)
  else:
    x, infected, infected_dd, prob1, prob0, score, tp, fp, fn, uncon_negs, determined,\
        overdetermined, surep, unsurep, wrongly_undetected,\
        num_infected_in_test = cs.decode_lasso(y, algo, prefer_recall=False,
            compute_stats=False)

  for test in range(t):
    if bool_y[test] > 0 and num_infected_in_test[test] == 0:
      print('y[%d] is infected but no infected people found' % test)
    if bool_y[test] == 0 and num_infected_in_test[test] > 0:
      print('y[%d] is not infected but infected people found' % test)
    
  #infected = cs.decode_comp_new1(bool_y)
  #infected_dd = np.zeros(n)
  #print(infected)
  #print(infected_dd)

  sure_list = []
  unsure_list = []
  neg_list = []
  for i in range(n):
    if infected_dd[i] == 1:
      sure_list.append(i+1)
    elif infected[i] == 1:
      unsure_list.append(i+1)
    else:
      neg_list.append(i+1)

  print('Surely infected: ', sure_list)
  print('Possibly infected: ', unsure_list)
  print('viral loads:', x)

  #print('Not infected: ', neg_list)

def print_results_COMP():
  #convert_matlab_format(M, sys.stdout)
  filenames = ['test1.txt', 'test2.txt', 'test3.txt', 'test4.txt', 'test5.txt']
  for i, filename in enumerate(filenames):
    filename = os.path.join(config.data_dir, "ncbs", filename)
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
  filenames = ['test1.txt', 'test2.txt', 'test3.txt', 'test4.txt', 'test5.txt']
  for i, filename in enumerate(filenames):
    filename = os.path.join(config.data_dir, "ncbs", filename)
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

def run_on_list_of_files():
  global n, t, M
  n = 40
  t = 16
  M = optimized_M_16_40_ncbs
  m, cts_list , positive_cts = find_max_positive_cycle_time()

  #print(m)
  np.set_printoptions(linewidth=800)
  p = 0.93
  ys = []
  bool_ys = []
  for cts in positive_cts:
    bool_y = (cts > 0).astype(np.int32)
    y = (1+p) ** (m - cts)
    y = y * bool_y
    #print(y)
    #print((m - cts) * bool_y)
    ys.append(y)
    bool_ys.append(bool_y)

  #for i, y in enumerate(ys):
  #  print('Test', i + 1)
  #convert_matlab_format(ys, sys.stdout)

  #sys.exit(1)

  #for bool_y in bool_ys:
  #  print_infected_people(bool_y, 'COMP')

  for i, y in enumerate(ys):
    print("\n=====================")
    print('Test %d:' % (i+1))
    bool_y = (y > 0).astype(np.int32)
    print_infected_people(y, bool_y, 'COMP')
    print_infected_people(y, bool_y, 'combined_COMP_NNOMP_random_cv')
    print_infected_people(y, bool_y, 'SBL')
    #print_infected_people(y, 'l1ls')
    #print_infected_people(y, 'NNOMP')

def run_harvard_data():
  global n, t, M
  n = 60
  t = 24
  M = optimized_M_3
  bool_y = [
      1, # D3
      1, # D4
      0, # D5
      0, # D6
      0, # D7
      0, # D8
      0, # D9
      0, # D10
      0, # E3
      0, # E4
      1, # E5
      0, # E6
      0, #E7
      1, #E8
      0, #E9
      0, #E10
      1, #F3
      0, #F4
      0, #F5
      0, #F6
      0, #F7
      0, #F8
      0, #F9
      0, #F10
      ]
  bool_y = np.array(bool_y, dtype=np.int32)
  y = bool_y
  print_infected_people(y, bool_y, 'COMP')

  df = pd.read_csv(os.path.join(config.data_dir, "harvard", "harvard_test1.csv"))
  #cts = df.values[:, 0]
  fl = df.values[:, 1:]

  assert t == fl.shape[1]
  cts = np.zeros(t)
  for j in range(t):
    if bool_y[j] == 0:
      print('Skipping', j)
      continue

    for i in range(fl.shape[0]):
      if fl[i, j] >= 1000:
        fl1 = fl[i-1, j]
        fl2 = fl[i, j]

        ct = i + math.log(1000/fl1) / math.log(fl2/fl1)
        print(j, ct)
        cts[j] = ct
        break

  #np.set_printoptions(linewidth=400)
  p = 0.95
  m = np.max(cts)      
  y = (1+p) ** (m - cts)
  y = y * bool_y
  print('cts:', cts)
  print('y:', y)
  print_infected_people(y, bool_y, 'SBL')
  print_infected_people(y, bool_y, 'combined_COMP_NNOMP_random_cv')

run_on_list_of_files()
run_harvard_data()
