import sys
sys.path.append(".")
from pickle_manager import *

if __name__ == '__main__':
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
  x = stats_manager.load("optimized_M_45_285_kirkman", "COMP", 3)
  print(x, type(x), len(x))


