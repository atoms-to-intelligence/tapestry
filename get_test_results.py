# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
### Get the results for a given test. Implemented in function get_test_results() ###


# The actual get_test_results() is present in app_utils.py
import app_utils
import config


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
#
# Also added the number of infections that the matrix can handle and the
# infection percentage. The value inside the dict is now a 3-tuple instead of
# a str.
MSizeToLabelDict = {
    "16x40":     ("optimized_M_16_40_ncbs", 3, 7.5),
    "21x70":     ("optimized_M_21_70_STS", 4, 6),
    #"24x60":     ("optimized_M_3", 4, 6),
    "45x105":    ("optimized_M_45_105_STS_1", 10, 10),
    "45x195":    ("optimized_M_45_195_STS_1", 10, 5),
    "63x399":    ("optimized_M_63_399_STS_1", 10, 2.5),
    "93x961":    ("optimized_M_93_961_STS_1", 10, 1),
    #"46x96":    ("optimized_M_46_96_1", 10, 10)
    #"46x192":   "optimized_M_46_192_1",
    }


# List of matrix labels and their codenames. Not all matrix labels have
# codenames
mat_codenames = {
    'optimized_M_16_40_ncbs':   'RABBIT',
    'optimized_M_3':            'FOX',  # 24x60
    'optimized_M_21_70_STS':    'BEAR',
    "optimized_M_45_105_STS_1": 'LION',
    "optimized_M_45_195_STS_1": 'TIGER',
    "optimized_M_63_399_STS_1": 'RHINO',
    "optimized_M_93_961_STS_1": 'CROC',
    }


# Returns a copy of size -> label dictionary
def get_matrix_sizes_and_labels():
  return dict(MSizeToLabelDict)


# Returns a copy of label -> matrix dictionary
def get_matrix_labels_and_matrices():
  return dict(MLabelToMatrixDict)

# Returns a copy of label -> matrix codename dictionary
def get_matrix_codenames():
  return dict(mat_codenames)

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
#
# Return result_string, sure_list, unsure_list, neg_list, x
#  * These will be stored in db.
#  * Only result_string will be displayed to the user
#
# sure_list - list of surely positive people
# unsure_list - list of possibly positive people
# neg_list - list of negative samples
# x - viral loads
#
# result_string - nicely formatted string containing above info. Will be seen
#       by the user
def get_test_results(matrix_label, cycle_times):
  M = get_matrix_for_label(matrix_label)

  sure_list, unsure_list, neg_list, x = app_utils.get_test_results(M, cycle_times)

  result_string = get_result_string_from_lists(sure_list, unsure_list, neg_list, x)

  res = {
      "result_string" :  result_string,
      "sure_list" :      sure_list,
      "unsure_list" :    unsure_list,
      "neg_list" :       neg_list, 
      "x" :              x,
      }
  return res


# Composes the result string from the list of surely positives, possibly
# positives, negatives and the x values
def get_result_string_from_lists(sure_list, unsure_list, neg_list, x):
  if not sure_list and not unsure_list:
    s0 = "No Positive Samples\n"
    s3 = "All samples are negative\n"

    s1 = ""
    s2 = ""
  else:
    s0 = ""
    if sure_list:
      s1  = "Surely Positive Samples: %s \n" % \
          ", ".join([str(item) for item in sure_list])
    else:
      s1  = "No Surely Positive Samples\n"

    if unsure_list:
      s2  = "Possibly Positive Samples: %s \n" % \
          ", ".join([str(item) for item in unsure_list])
    else:
      s2  = ""

    if not neg_list:
      s3 = "No surely negative samples detected\n"
    else:
      s3 = "Remaining samples are negative\n"

  x_str = ""
  if config.app_algo != 'COMP':
    x_str = "Detected viral loads: %s\n" % \
            ", ".join(["%.3f" % (item) for item in x])
  result_string = s0 + s1 + s2 + s3 + x_str

  return result_string


# Code to run at deployment time. Checks for errors such as invalid matrix
# labels or sizes. 
#
# *******    If this raises an exception, then deployment should not  ********
# *******            happen and a notification must be sent.          ********
def at_deployment():
  sanity_check_for_matrices()
  api_sanity_checks()
  test_harvard_data()
  fake_data_test()





###########      Internal Code for testing      ##############

def test_get_result_string_from_lists():
  tmp = config.app_algo
  for algo in ['COMP', 'SBL']:
    config.app_algo = algo
    print("Using algo: %s\n" % config.app_algo)

    sure_list = [2, 3, 4]
    unsure_list = [1, 7]
    neg_list = list(range(1, 11))
    for item in (sure_list + unsure_list):
      neg_list.remove(item)
    x = np.random.rand(11)
    mask = np.zeros(11)
    mask[sure_list + unsure_list] = 1
    x = x * mask
    x = x[1:]
    res = get_result_string_from_lists(sure_list, unsure_list, neg_list, x)
    print(res)

    sure_list = []
    unsure_list = [1, 7]
    neg_list = list(range(1, 11))
    for item in (sure_list + unsure_list):
      neg_list.remove(item)
    x = np.random.rand(11)
    mask = np.zeros(11)
    mask[sure_list + unsure_list] = 1
    x = x * mask
    x = x[1:]
    res = get_result_string_from_lists(sure_list, unsure_list, neg_list, x)
    print(res)

    sure_list = [2, 3, 4]
    unsure_list = []
    neg_list = list(range(1, 11))
    for item in (sure_list + unsure_list):
      neg_list.remove(item)
    x = np.random.rand(11)
    mask = np.zeros(11)
    mask[sure_list + unsure_list] = 1
    x = x * mask
    x = x[1:]
    res = get_result_string_from_lists(sure_list, unsure_list, neg_list, x)
    print(res)

    sure_list = []
    unsure_list = []
    neg_list = list(range(1, 11))
    for item in (sure_list + unsure_list):
      neg_list.remove(item)
    x = np.random.rand(11)
    mask = np.zeros(11)
    mask[sure_list + unsure_list] = 1
    x = x * mask
    x = x[1:]
    res = get_result_string_from_lists(sure_list, unsure_list, neg_list, x)
    print(res)

    sure_list = list(range(1, 5))
    unsure_list = list(range(5, 11))
    neg_list = []
    for item in (sure_list + unsure_list):
      if item in neg_list:
        neg_list.remove(item)
    x = np.random.rand(11)
    mask = np.zeros(11)
    mask[sure_list + unsure_list] = 1
    x = x * mask
    x = x[1:]
    res = get_result_string_from_lists(sure_list, unsure_list, neg_list, x)
    print(res)
  config.app_algo = tmp


# Go through all the labels in MSizeToLabelDict and determine if the corresponding
# matrices are present MLabelToMatrixDict. Checks if the sizes match up.
def sanity_check_for_matrices():
  print('Printing out matrix labels for each size...\n')
  for msize in MSizeToLabelDict:
    t, n = [int(item) for item in msize.split("x")]
    mlabel, d, percent = MSizeToLabelDict[msize]
    print(f'{n:3} Samples ({t} tests. Upto {d:2} infections, or '
        f'{percent}% infection rate)')

  print('\nChecking if matrix sizes correspond to an actual label...\n')
  error = False
  for msize in MSizeToLabelDict:
    t, n = [int(item) for item in msize.split("x")]
    mlabel, d, percent = MSizeToLabelDict[msize]
    if mlabel in MLabelToMatrixDict:
      print(msize, mlabel, 'label exists')
      M = MLabelToMatrixDict[mlabel]
      if t == M.shape[0] and  n == M.shape[1]:
        print(msize, mlabel, 'Sizes match')
      else:
        print(msize, mlabel, 'Sizes don\'t match:',
            M.shape, "       <------------")
        error = True
    else:
      print(msize, mlabel, 'label does not exists in MLabelToMatrixDict       <------------')
      error = True

    #print("\nChecking if each matrix for given size has a codename\n")
    if mlabel in mat_codenames:
      print(f'{msize} matrix {mlabel} codename {mat_codenames[mlabel]}')
    else:
      print(f'No codename exist for {msize} matrix {mlabel}       <------------')
      error = True
  if error:
    print("\nGot Error :(\n")
    raise ValueError("Some error in matrix setup")
  else:
    print("\nAll OK\n")
  

def api_sanity_checks():
  error = False
  print('\nSome API checks...\n')
  print('Label for valid size 45x105', get_current_matrix_label_for_size("45x105"))
  try:
    print(get_current_matrix_label_for_size("fslkj"))
    print('Did not get expected error while calling with invalid matrix size')
    error = True
  except KeyError as e:
    print('Got expected error while calling with invalid matrix size:', str(e))

  
  print('Shape of the valid matrix optimized_M_16_40_ncbs', get_matrix_for_label("optimized_M_16_40_ncbs").shape)
  try:
    invalid_label = "ojlkj"
    print('shape of the invalid matrix ', invalid_label,
        get_matrix_for_label(invalid_label).shape)
    print('Did not get expected error while calling with invalid matrix size')
    error = True
  except KeyError as e:
    print('Got expected error while calling with invalid matrix label:', str(e))

  # Dict retrieval API
  print("\nDictionary Retrieval API sanity checks...\n")
  d1 = get_matrix_sizes_and_labels()
  d2 = get_matrix_labels_and_matrices()
  d3 = get_matrix_codenames()
  print(d1)
  print(d2.keys())
  if d1 != MSizeToLabelDict:
    error = True
    print("**********       Not returning MSizeToLabelDict in get_matrix_sizes_and_labels()      *********")

  if d1 is MSizeToLabelDict:
    error = True
    print("**********      Return a copy of MSizeToLabelDict in get_matrix_sizes_and_labels()      ***********")

  if d2 != MLabelToMatrixDict:
    error = True
    print("**********      Not returning MLabelToMatrixDict in get_matrix_labels_and_matrices()      ************")

  if d2 is MLabelToMatrixDict:
    error = True
    print("**********      Return a copy of MLabelToMatrixDict in get_matrix_labels_and_matrices()      ************")

  if d3 != mat_codenames:
    error = True
    print("**********       Not returning mat_codenames in get_matrix_sizes_and_labels()      *********")

  if d3 is mat_codenames:
    error = True
    print("**********      Return a copy of mat_codenames in get_matrix_sizes_and_labels()      ***********")

  l1 = list(d2.keys())[0]
  print(l1, type(d2[l1]), d2[l1].shape)

  if error:
    print("\nGot Error :(\n")
    raise ValueError("Some error in matrix setup")
  else:
    print("\nAll OK\n")


# Tests get_test_results on Harvard experimental data using the currently
# configured algorithm. Useful for sanity check. Should not throw an
# exception.
def test_harvard_data():
  from experimental_data_manager import read_harvard_data_cts
  print('Testing Harvard data with config.app_algo =', config.app_algo)
  pos_idx, cts = read_harvard_data_cts()
  res = get_test_results("optimized_M_3", cts)
  result_string = res["result_string"]
  sure_list = res["sure_list"]
  unsure_list = res["unsure_list"]
  neg_list = res["neg_list"]
  x = res["x"]
  print(result_string)
  print(sure_list)
  print(unsure_list)
  print(neg_list)
  print(x)
  pos_list = res['sure_list'] + res['unsure_list']
  for idx in pos_idx:
    assert idx in pos_list


def fake_data_test():
  from experimental_data_manager import get_random_fake_test_data
  size_to_label_dict = get_matrix_sizes_and_labels()
  label_to_M_dict = get_matrix_labels_and_matrices()
  for msize in size_to_label_dict:
    mlabel, d, percent = size_to_label_dict[msize]
    bool_x, cts = get_random_fake_test_data(msize, mlabel)
    res = get_test_results(mlabel, cts)
    print("Results for data faked for %s matrix %s" % (msize, mlabel))
    print('bool_x:', bool_x)
    print(res)
    pos_list = res['sure_list'] + res['unsure_list']
    for idx in bool_x:
      assert idx in pos_list

if __name__ == '__main__':
  at_deployment()

  #import numpy as np
  #test_get_result_string_from_lists()

  
