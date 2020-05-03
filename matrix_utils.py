# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8

import config
import os

# Loads the matrices in the provided dict
def load_extra_mats(variables, extra_mlabels):
  for name in sorted(os.listdir(config.extra_mat_dir)):
    full_path = os.path.join(config.extra_mat_dir, name)
    vname = name[:-4]
    extra_mlabels.append(vname)
    print(name, vname, full_path)
    # is filename .txt or .mat?
    if name[-4:] == ".txt":
      print("txt")
      #variables[vname] = np.loadtxt(full_path)
      variables[vname] = "foo"
    elif name[-4:] == ".mat":
      print("mat")
      variables[vname] = "bar"
    else:
      raise ValueError('Invalid extension')

if __name__ == '__main__':
  d = {}
  extra_mlabels = []
  load_extra_mats(d, extra_mlabels)
  print(d)
  print(extra_mlabels)
