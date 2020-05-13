import sys
sys.path.append('.')

from core.cs_expts import *

if __name__=='__main__':
  #large_test_decode_comp_combined(1000)
  #mr = None
  #do_many_expts(200, 6, 46, num_expts=100, M=None,
  #    add_noise=True,algo='combined_COMP_NNOMP_random_cv', mr=mr)
  #compare_different_ns()
  #M = [optimized_M_45_105_STS_1, optimized_M_45_285_social_golfer[:, :105]]
  #mlabels = ['optimized_M_45_105_STS_1', 'optimized_M_45_285_social_golfer[:, :105]']
  #M = [optimized_M_45_195_STS_1, optimized_M_45_285_social_golfer[:, :195]]
  #mlabels = ['optimized_M_45_195_STS_1', 'optimized_M_45_285_social_golfer[:, :195]']

  #compare_different_mats(M, mlabels)
  #run_many_parallel_expts()
  #run_stats_for_these_matrices(
  #    massive_pooling_matrices,
  #    save=False
  #  )
  #print(kirkman_mlabels)
  #run_stats_for_these_matrices(
  #    kirkman_mlabels,
  #    save=True
  #  )
  #run_stats_for_these_matrices(
  #    [
  #        'optimized_M_45_105_kirkman',
  #        'optimized_M_36_180_kirkman',
  #    ],
  #    save=True
  #  )
  run_stats_for_these_matrices1()


  #compare_sts_vs_kirkman()
  #for mr in range(8, 15):
  #  do_many_expts(40, 2, 16, num_expts=1000, M=optimized_M_2,
  #      add_noise=True,algo='NNOMP_random_cv', mr=mr)
  #do_many_expts(40, 4, 16, num_expts=1000, M=optimized_M_2,
  #    add_noise=True,algo='NNOMP_loo_cv')
  #do_many_expts(40, 3, 16, num_expts=1000, M=optimized_M_2,
  #    add_noise=True,algo='NNOMP')

# Following code assumes the same Gaussian noise for each y
#for noise in [0.0001, 0.0002, 0.0004, 0.0008, 0.001, 0.002, 0.004, 0.008, 0.01]:
#  do_many_expts(40, 2, 16, num_expts=1000, M=optimized_M_2, add_noise=True,
#      cross_validation=False, algo='OMP', noise_magnitude=noise)



