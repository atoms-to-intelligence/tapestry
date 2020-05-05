# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=2
import sys
sys.path.append(".")

from core import config

import string
import numpy as np
import os

def parse_social_golfer_wolfram():
  row_map = {}
  upper = string.ascii_uppercase
  lower = string.ascii_lowercase
  for i in range(26):
    row_map[upper[i]] = i

  row_map['+'] = 26

  #for j in range(18):
  #  row_map[lower[j]] = j + 18

  print(row_map)

  #M = np.zeros((36, 99), dtype=np.int32)
  #M = np.zeros((18, 153), dtype=np.int32)
  M = np.zeros((27, 117), dtype=np.int32)
  with open(os.path.join(config.unparsed_mat_dir, 'social_golfers', '27x117_social_golfer.txt')) as f:
    lines = f.readlines()
    col = 0
    for line in lines:
      if line.strip() == "":
        continue
      rows = [row_map[ch] for ch in line.strip()]
      print(rows)
      M[rows, col] = 1
      col += 1
      
    print(col)
    assert col == 117

  write_matrix(M, "optimized_M_27_117_social_golfer.txt")

def write_matrix(M, fname):
  print(np.sum(M))
  print(np.sum(M, axis=1))
  print(np.sum(M, axis=0))

  full_path = os.path.join(config.mat_dir, fname)
  np.savetxt(full_path, M, fmt="%d")

def parse_social_golfer_github(t, n, block_size):
  count = 0
  #f = open("48x384_social_golfer.txt")
  #f = open("384_4096.txt")
  f = open(os.path.join(config.unparsed_mat_dir, 'social_golfers', f"{t}_{n}.txt"), "r")
  A = np.zeros((t,n), dtype=np.int32)
  for line in f:
    if line.strip() == '':
      continue
    rows = [int(item) for item in line.strip().split() if item]
    if len(rows) != 12:
      print("Last row")
    print(rows)
    
    block = count // block_size
    
    idx = count % block_size
    
    w1 = block*4
    w2 = w1 + 1
    w3 = w2 + 1
    w4 = w3 + 1
    
    col1 = block_size*w1 + idx
    col2 = block_size*w2 + idx
    col3 = block_size*w3 + idx
    col4 = block_size*w4 + idx
    
    try:
      A[rows[0], col1] = 1
      A[rows[1], col1] = 1
      A[rows[2], col1] = 1

      A[rows[3], col2] = 1
      A[rows[4], col2] = 1
      A[rows[5], col2] = 1

      A[rows[6], col3] = 1
      A[rows[7], col3] = 1
      A[rows[8], col3] = 1

      A[rows[9], col4] = 1
      A[rows[10], col4] = 1
      A[rows[11], col4] = 1
    except:
      print("Last row")
      pass
    count += 1

  print(count*4, n)
  #assert count*4 == n
  write_matrix(A, f"optimized_M_{t}_{n}_social_golfer.txt")

parse_social_golfer_wolfram()
parse_social_golfer_github(96, 1312, 32)
#parse_social_golfer_github(45, 285, 15)


