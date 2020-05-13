# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
if __name__=='__main__':
  import sys
  # This is a hack needed so that you can directly run this file as 
  # python3 core/test/test_get_test_results.py.
  # 
  # This is needed because we want to import "matrices" from two levels up, but python
  # does not know about one level up. Please only use this hack in test files
  # and not in files which are meant to be imported from other code.
  sys.path.append(".")

from core.get_test_results import *
from utils.experimental_data_manager import get_random_fake_test_data

# Code should be robust to invalid y's
def invalid_y_test():
  mlabel = "optimized_M_16_40_ncbs"
  t = 16
  pos_idx, cts = get_random_fake_test_data('16x40', mlabel)
  print(pos_idx, cts)
  
  # First do test with valid cts
  res = get_test_results(mlabel, cts)
  print(res['result_string'])

  # Now test with invalid cts
  cts_invalid = np.zeros(t) + 100
  cts_invalid[0] = 0.
  cts_invalid[1] = 0.
  cts_invalid[2] = 0.
  cts_invalid[3] = 0.
  #cts[1] = 20.
  #cts[3] = 20.

  # Test with invalid cts
  res = get_test_results(mlabel, cts_invalid)
  print(res['result_string'])

  # More invalid cts. this ensures that some entries remain in the matrix, but A
  # has 0 support for some positive tests
  cts_invalid[4] = 0.
  cts_invalid[5] = 0.

  # Test with invalid cts
  res = get_test_results(mlabel, cts_invalid)
  print(res['result_string'])

invalid_y_test()

#matrix_pdfs_sanity_check()
#sanity_check_for_matrices()
#test_harvard_data()
#api_sanity_checks()
#fake_data_test()
#at_deployment()

#import numpy as np
#test_get_result_string_from_lists()


