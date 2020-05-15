import sys
sys.path.append(".")

from utils.app_utils import *
if __name__ == '__main__':
	from core.matrices import MDict
	from utils.experimental_data_manager import read_standard_cts
	mat_label = sys.argv[1]
	data_file = sys.argv[2]
	cts = read_standard_cts(data_file)
	M = MDict[mat_label]
	for algo in config.app_algos:
		#config.app_algo = algo
		print(algo)
		sure_list, unsure_list, neg_list, x = get_test_results(M, cts, algo)
		print("Sure list:", sure_list)
		print("Unsure list:", unsure_list)
		print("Neg list:", neg_list)
		print("x:", x)
