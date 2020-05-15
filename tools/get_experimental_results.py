import sys
sys.path.append(".")

from utils.app_utils import *

if __name__ == '__main__':
	from utils.experimental_data_manager import read_standard_cts
	from core.get_test_results import get_test_results
	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument('-M', '--matrix')
	parser.add_argument('-ct', '--cycle-times')
	parser.add_argument('-n', default=None)
	args = parser.parse_args()
	cts = read_standard_cts(args.cycle_times)
	
	print(get_test_results(args.matrix, cts, args.n)["result_string"])
