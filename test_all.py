from get_test_results import at_deployment
from cs_expts import run_many_parallel_expts_internal, do_many_expts
from matrices import MDict

def test_get_test_results():
  at_deployment()

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
  for mlabel in MDict:
    print("Running algos for matrix ", mlabel)
    M = MDict[mlabel]
    n = M.shape[1]
    t = M.shape[0]
    num_expts=1
    add_noise = True
    d_range = list(range(1,11))
    n_jobs = len(d_range)

    run_many_parallel_expts_internal(num_expts, n, t, add_noise, M, algos, d_range, n_jobs)
    #for d in d_range:
    #  do_many_expts(n, d, t, num_expts=num_expts, M=M, add_noise=add_noise, algo=algos)

if __name__ == '__main__':
  test_get_test_results()
  test_all_matrices_and_algos()
