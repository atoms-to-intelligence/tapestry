# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
import sys
sys.path.append(".")

from utils.app_utils import *
if __name__ == '__main__':
  # These imports and functions are for testing purposes only
  from core.matrices import *
  from utils.experimental_data_manager import read_harvard_data_cts, read_standard_cts

  # Test using Harvard and ncbs data
  def run_tests():
    pos_idx, cts = read_harvard_data_cts()
    print('cts:')
    print(cts)
    y = get_y_from_cts(cts)
    print('y:')
    print(y)

    for algo in config.app_algos:
      print("\nUsing algorithm %s \n" % algo)
      harvard_test(algo)

  def harvard_test(algo):
    print("Test on Harvard Dataset")
    pos_idx, cts = read_harvard_data_cts()
    M = optimized_M_3
    sure_list, unsure_list, neg_list, x = get_test_results(M, cts, algo)
    print("Sure list:", sure_list)
    print("Unsure list:", unsure_list)
    print("Neg list:", neg_list)
    print("x:", x)

  def harvard_30_120_test():
    data_file = os.path.join(config.data_dir, 'harvard', "harvard_30_120.txt") 
    cts = read_standard_cts(data_file)
    M = optimized_M_30_120_kirkman
    for algo in config.app_algos:
      print(algo)
      sure_list, unsure_list, neg_list, x = get_test_results(M, cts, algo)
      print("Sure list:", sure_list)
      print("Unsure list:", unsure_list)
      print("Neg list:", neg_list)
      print("x:", x)


  def harvard_90_1140_test():
    data_file = os.path.join(config.data_dir, 'harvard', "harvard_90_1140.txt") 
    cts = read_standard_cts(data_file)
    M = optimized_M_90_1140_kirkman
    for algo in config.app_algos:
      print(algo)
      sure_list, unsure_list, neg_list, x = get_test_results(M, cts, algo)
      print("Sure list:", sure_list)
      print("Unsure list:", unsure_list)
      print("x:", x)

  run_tests()
  harvard_30_120_test()
  harvard_90_1140_test()


