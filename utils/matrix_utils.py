# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8

from core import config

import os
import numpy as np
import scipy.io as sio

# Loads the matrices in the provided dict
def load_extra_mats(variables, extra_mlabels):
  for name in sorted(os.listdir(config.extra_mat_dir)):
    full_path = os.path.join(config.extra_mat_dir, name)
    vname = name[:-4]
    extra_mlabels.append(vname)
    #print(name, vname, full_path)
    # is filename .txt or .mat?
    if name[-4:] == ".txt":
      print("txt")
      variables[vname] = np.loadtxt(full_path)
      #variables[vname] = "foo"
    elif name[-4:] == ".mat":
      mat = sio.loadmat(full_path)
      for key in mat:
        if type(mat[key]) == np.ndarray:
          variables[vname] = mat[key]
    else:
      raise ValueError('Invalid extension')

if __name__ == '__main__':
  raise ValueError('This is a library file. Please call test_matrix_utils.py instead')
