# Get the results for a given test
#
# Input: Cycle time vector. Numpy array of size 't'
# Input: Matrix size on which the test was performed. String such as "46x96",
#       "46x192" etc.

# Import the dictionary of matrix labels to the actual numpy matrices
from matrices import MDict

# Dictionary of matrix size to the label of the matrix which will be used
MSizeDict = {
    "46x96":    "optimized_M_46_96_1",
    "46x192":   "optimized_M_46_192_1",
    }

def get_current_matrix_label_for_size(matrix_size):
  if not matrix_size in MSizeDict:
    raise KeyError('Invalid matrix size: "%s". Correct matrix sizes are: "%s"' %
        (matrix_size, '", "'.join(MSizeDict.keys())))

  return MSizeDict[matrix_size]

def get_matrix_for_label(matrix_label):
  return MDict[matrix_label]

def get_test_results(matrix_size, cycle_times):
  result_string = "This is a stub"
  return result_string

if __name__ == '__main__':
  # Do some sanity check tests
  print(get_current_matrix_label_for_size("46x96"))
  print(get_current_matrix_label_for_size("45x96"))

