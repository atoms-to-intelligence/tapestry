### Get the results for a given test. Implemented in function get_test_results() ###


# Import the dictionary of matrix labels to the actual numpy matrices
#
# This dictionary is **auto-generated** by looking at all the variables named
# optimized_M_* in matrices.py and matrices1.py
from matrices import MDict as MLabelToMatrixDict

# Dictionary of matrix size to the label of the matrix which will be used. 
# This dictionary should be modified by hand.
#
# Matrix Sizes are strings such as "46x96", "46x192"
#
# Matrix labels are the actual variable names used for the matrices in
# matrices.py and matrices1.py
#
# Warning: during deployment, once a matrix size is added, it should never be
# removed. Once a matrix label is added, it should never be removed from the
# code. 
MSizeToLabelDict = {
    "46x96":    "optimized_M_46_96_1",
    "46x192":   "optimized_M_46_192_1",
    }

# Returns a copy of size -> label dictionary
def get_matrix_sizes_and_labels():
  return dict(MSizeToLabelDict)

# Returns the currently used matrix label for a given size
def get_current_matrix_label_for_size(matrix_size):
  if not matrix_size in MSizeToLabelDict:
    raise KeyError('Invalid matrix size: "%s". Correct matrix sizes are: "%s"' %
        (matrix_size, '", "'.join(MSizeToLabelDict.keys())))

  return MSizeToLabelDict[matrix_size]

# Returns the matrix (which is a numpy array of shape (t, n)) corresponding to
# the given label
def get_matrix_for_label(matrix_label):
  if not matrix_label in MLabelToMatrixDict:
    raise KeyError('Invalid matrix label: "%s". Correct matrix labels are: "%s"' %
        (matrix_label, '", "'.join(MLabelToMatrixDict.keys())))
  return MLabelToMatrixDict[matrix_label]

# Get the results for a test
#
# Input: Matrix label on which the test was performed. Such as
# "optimized_M_46_96_1"
#
# Input: Cycle time vector. Numpy array of size 't', which is the number of
# tests or equivalently the number of rows in the matrix
def get_test_results(matrix_label, cycle_times):
  result_string = "This is a stub"
  return result_string

# Go through all the labels in MSizeToLabelDict and determine if the corresponding
# matrices are present MLabelToMatrixDict. Checks if the sizes match up.
def sanity_check_for_matrices():
  for msize in MSizeToLabelDict:
    t, n = [int(item) for item in msize.split("x")]
    mlabel = MSizeToLabelDict[msize]
    assert mlabel in MLabelToMatrixDict
    M = MLabelToMatrixDict[mlabel]
    assert t == M.shape[0]
    assert n == M.shape[1]

  print("All OK")

if __name__ == '__main__':
  sanity_check_for_matrices()
  # Do some sanity check tests
  #print(get_current_matrix_label_for_size("46x96"))
  #print(get_current_matrix_label_for_size("45x96"))

