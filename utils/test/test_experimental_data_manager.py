import sys
sys.path.append(".")

from utils.experimental_data_manager import *

if __name__ == '__main__':
  cts = get_israel_data_cts()
  print(cts, len(cts))
  A = parse_israel_matrix()
  #print(A)
  print(np.sum(A))
  print(np.sum(A, axis=1))
  print(np.sum(A, axis=0))
  from utils.app_utils import get_test_results
  from core.get_test_results import get_result_string_from_lists

  n = 384
  config.app_algo = 'SBL'
  sure_list, unsure_list, neg_list, x = get_test_results(A, cts)
  #print('sure_list:', sure_list)
  #print('unsure_list:', unsure_list)
  #print('neg_list:', neg_list)
  #print('x:', x)
  result_string = get_result_string_from_lists(sure_list, unsure_list,
      neg_list, x, n)
  print(result_string)
  #pos_idx, cts = read_harvard_data_cts()
  #print('Harvard data (24x60)')
  #print('Positive indicies:', pos_idx)
  #print('Cycle times:', "[", ", ".join(["%.1f" % item for item in cts]), "]")

  #print('\nFake Data (16x40)\n')
  #for i in range(10):
  #  pos_idx, cts = get_random_fake_test_data('16x40', 'optimized_M_16_40_ncbs')
  #  print('Positive indicies:', pos_idx)
  #  print('Cycle times:', "[", ", ".join(["%.1f" % item for item in cts]), "]")


