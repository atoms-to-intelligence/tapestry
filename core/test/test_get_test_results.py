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

#matrix_pdfs_sanity_check()
#sanity_check_for_matrices()
test_harvard_data()
#api_sanity_checks()
fake_data_test()
#at_deployment()

#import numpy as np
#test_get_result_string_from_lists()


