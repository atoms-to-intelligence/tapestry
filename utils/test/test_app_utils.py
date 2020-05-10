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

    for algo in ['COMP', 'SBL', 'combined_COMP_NNOMP_random_cv']:
      config.app_algo = algo
      print("\nUsing algorithm %s \n" % config.app_algo)
      harvard_test()

  def harvard_test():
    print("Test on Harvard Dataset")
    pos_idx, cts = read_harvard_data_cts()
    M = optimized_M_3
    sure_list, unsure_list, neg_list, x = get_test_results(M, cts)
    print("Sure list:", sure_list)
    print("Unsure list:", unsure_list)
    print("Neg list:", neg_list)
    print("x:", x)

  def harvard_30_120_test():
    data_file = os.path.join(config.data_dir, 'harvard', "harvard_30_120.txt") 
    cts = read_standard_cts(data_file)
    M = optimized_M_30_120_kirkman
    config.app_algo = 'COMP'
    print('COMP')
    sure_list, unsure_list, neg_list, x = get_test_results(M, cts)
    print("Sure list:", sure_list)
    print("Unsure list:", unsure_list)
    print("Neg list:", neg_list)
    print("x:", x)

    config.app_algo = 'combined_COMP_SBL'
    print(config.app_algo)
    sure_list, unsure_list, neg_list, x = get_test_results(M, cts)
    print("Sure list:", sure_list)
    print("Unsure list:", unsure_list)
    print("Neg list:", neg_list)
    print("x:", x)

  def harvard_90_1140_test():
    data_file = os.path.join(config.data_dir, 'harvard', "harvard_90_1140.txt") 
    cts = read_standard_cts(data_file)
    M = optimized_M_90_1140_kirkman
    config.app_algo = 'COMP'
    print('COMP')
    sure_list, unsure_list, neg_list, x = get_test_results(M, cts)
    print("Sure list:", sure_list)
    print("Unsure list:", unsure_list)
    #print("Neg list:", neg_list)

    x_list = []
    np.set_printoptions(threshold=10000)
    config.app_algo = 'combined_COMP_SBL'
    print(config.app_algo)
    sure_list, unsure_list, neg_list, x = get_test_results(M, cts)
    print("Sure list:", sure_list)
    print("Unsure list:", unsure_list)
    print("sure + unsure", sorted(sure_list + unsure_list))
    print("x:", x[x>0])
    x_list.append(np.expand_dims(x, -1))

    x_list = []
    np.set_printoptions(threshold=10000)
    config.app_algo = 'SBL'
    print(config.app_algo)
    sure_list, unsure_list, neg_list, x = get_test_results(M, cts)
    print("Sure list:", sure_list)
    print("Unsure list:", unsure_list)
    print("sure + unsure", sorted(sure_list + unsure_list))
    print("x:", x[x>0])
    x_list.append(np.expand_dims(x, -1))

    #config.app_algo = 'combined_COMP_NNOMP_random_cv'
    #print(config.app_algo)
    #sure_list, unsure_list, neg_list, x = get_test_results(M, cts)
    #print("Sure list:", sure_list)
    #print("Unsure list:", unsure_list)
    #print("sure + unsure", np.where(x)[0] + 1)
    #print("x:", x[x>0])
    #x_list.append(np.expand_dims(x, -1))

    #config.app_algo = 'combined_COMP_l1ls_cv'
    #min_or_max = 'max'
    #print(config.app_algo)
    #sure_list, unsure_list, neg_list, x = get_test_results(M, cts)
    #print("Sure list:", sure_list)
    #print("Unsure list:", unsure_list)
    #print("sure + unsure", np.where(x)[0] + 1)
    #print("x:", x[x>0])
    #x_list.append(np.expand_dims(x, -1))

    #x_combined = np.concatenate(x_list, axis=1)
    #print(x_combined.shape)
    #np.savetxt("x_combined_harvard_90_1140.csv", x_combined, fmt="%.6f",
    #    delimiter=",")


  #run_tests()
  #harvard_30_120_test()
  harvard_90_1140_test()


