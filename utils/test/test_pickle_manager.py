# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
if __name__=='__main__':
  import sys
  # This is a hack needed so that you can directly run this file as 
  # python inbuilt_algos/nnompcv.py.
  # 
  # This is needed because we want to import "matrices" from one level up, but python
  # does not know about one level up. Please only use this hack in test files
  # and not in files which are meant to be imported from other code.
  sys.path.append(".")

from utils.pickle_manager import *

#stats_manager.do_directory_checks("optimized_M_3", "COMP", 3)
#x = stats_manager.load("optimized_M_3", "COMP", 3)
#print(x, type(x), len(x))
stats_manager.save("TEST", "COMP", 3, [1, 2, 3, 4, 5])
x = stats_manager.load("TEST", "COMP", 3)
print(x, type(x), len(x))
stats_manager.save("TEST", "COMP", 3, ["apple", "orange"])
x = stats_manager.load("TEST", "COMP", 3)
print(x, type(x), len(x))

stats_manager.save("TEST", "SBL", 3, ["china", "japan"])
x = stats_manager.load("TEST", "COMP", 3)
print(x, type(x), len(x))

stats_manager.save("TEST", "COMP", 5, ["haathi", "ghoda"])
x = stats_manager.load("TEST", "COMP", 5)
print(x, type(x), len(x))

x = stats_manager.load("TEST", "SBL", 3)
print(x, type(x), len(x))

stats_manager.save("HAATHI", "SBL", 3, ["apple", "orange"])
x = stats_manager.load("HAATHI", "SBL", 3)
print(x, type(x), len(x))


