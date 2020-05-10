# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
import sys
sys.path.append(".")

from core.matrices import *

if __name__ == '__main__':
  #print(extra_mlabels)
  #for mlabel in extra_mlabels:
  #  print_matrix_stats(extra_mats_dict[mlabel], mlabel)

  #validate_kirkman()
  #sys.exit(1)
  A = np.arange(8).reshape((2,4)) + 1
  #A = sts.sts(45)
  #n_strides = 4
  #B = strided_randomized_matrix(A, n_strides)
  #print(A)
  #print(B)
  #print(A.shape)
  #print(B.shape)
  #sys.exit(1)
  #int_coded_matrix(4, 15)
  #print(int_coded_M_6_63.shape)
  # Uncomment following code if you want to convert these matrices to matlab

  #mat = np.array([[1, 2, 3, 4], [3, 4, 5, 6], [10, 11, 12, 13]])
  #print(mat)
  #convert_matlab_format(optimized_M_46_96_1, sys.stdout,
  #    start_symbol='',
  #    end_symbol='',
  #    row_end_symbol=''
  #    )
  #sys.exit(1)
  #names = ['optimized_M_2.txt', 'optimized_M_3.txt', 'optimized_M_5.txt']
  #for mat, name in zip([optimized_M_2, optimized_M_3, optimized_M_5], names):
  #  with open(name, 'w') as f:
  #    convert_matlab_format(mat, f)

  # Uncomment following to print stats about these matrices

  print(optimized_M_1.shape)
  print(optimized_M_2.shape)
  print(optimized_M_3.shape)
  print(optimized_M_4.shape)
  print(optimized_M_5.shape)
  print(np.sum(optimized_M_1))
  print(np.sum(optimized_M_2))
  print(np.sum(optimized_M_3))
  print(np.sum(optimized_M_4))
  print(np.sum(optimized_M_5))
  #MList = [item for item in dir() if item.startswith("optimized_M_")]
  #for M in [optimized_M_1, optimized_M_2, optimized_M_3, optimized_M_4]:
  variables = globals()
  for m in MDict:
    M = MDict[m]
    row_sums = np.sum(M, axis=1)
    col_sums = np.sum(M, axis=0)
    row_sparsity = np.max(row_sums)
    min_row_sparsity = np.min(row_sums)
    col_sparsity = np.max(col_sums)
    min_col_sparsity = np.min(col_sums)
    total_sparsity = np.sum(row_sums)
    print(m, M.shape, 'Row sparsity: ', row_sparsity, 'Col sparsity: ', col_sparsity,
        'Min Row sparsity: ', min_row_sparsity, 'Min Col sparsity', min_col_sparsity,
        'Total sparsity: ', total_sparsity)

    #print(max(sums))


