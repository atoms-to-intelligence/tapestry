# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
import sys
sys.path.append(".")

from get_test_results import at_deployment, mat_codenames
from core.cs_expts import run_many_parallel_expts_internal, do_many_expts
from core.matrices import MDict, validate_kirkman

def test_get_test_results():
  at_deployment()

def test_kirkman_matrices():
  validate_kirkman()

def test_all_matrices_and_algos():
  algos = []
  algos.extend(['COMP'])
  algos.append('NNOMP')
  algos.append('combined_COMP_NNOMP')
  algos.append('NNOMP_random_cv')
  algos.extend(['combined_COMP_NNOMP_random_cv'])
  algos.append('SBL')
  algos.append('combined_COMP_SBL')
  algos.append('l1ls')
  algos.append('combined_COMP_l1ls')
  #algos.append('l1ls_cv')
  #algos.append('combined_COMP_l1ls_cv')
  # List of all deployed matrices past and present is in mat_codenames
  for mlabel in mat_codenames:
    print("Running algos for matrix ", mlabel)
    M = MDict[mlabel]
    n = M.shape[1]
    t = M.shape[0]
    num_expts=1
    add_noise = True
    d_range = list(range(0,11))
    n_jobs = 4 #len(d_range)

    run_many_parallel_expts_internal(num_expts, n, t, add_noise, M, algos,
            d_range, n_jobs, xslist=[None for d in d_range], mlabel=mlabel)
    #for d in d_range:
    #  do_many_expts(n, d, t, num_expts=num_expts, M=M, add_noise=add_noise, algo=algos)

if __name__ == '__main__':
  test_kirkman_matrices()
  test_get_test_results()
  test_all_matrices_and_algos()
