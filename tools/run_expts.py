# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8
import sys
sys.path.append('.')

from core.cs_expts import *

def run_stats_for_these_matrices1():
  labels = [
      #"optimized_M_27_117_kirkman",
      #"optimized_M_45_330_STS",
      #"optimized_M_48_320_kirkman",
      #"optimized_M_60_500_kirkman",
      #"optimized_M_60_500_kirkman",
      #"optimized_M_63_546_kirkman",
      #"optimized_M_69_667_kirkman",
      #"optimized_M_75_800_kirkman",
      #"optimized_M_93_1240_kirkman",
      #"optimized_M_192_5120_social_golfer",
      #"optimized_M_36_180_kirkman",
      #"optimized_M_81_918_kirkman",
      #"optimized_M_45_105_kirkman",
      "optimized_M_24_80_kirkman",
      ]
  ns = [
      #117,
      #300,
      #304,
      #300,
      #500,
      #500,
      #504,
      #506,
      #500,
      #1000,
      #1024,
      #72,
      #108,
      #105,
      56,
      ]
  mats = [MDict[label][:, :n] for (label, n) in zip(labels, ns) ]
  for i, n in enumerate(ns):
    labels[i] = f"{labels[i]}[:, :{n}]"
  #d_ranges = [ list(range(1, 16)) + [20, 25, 30, 35, 40] for item  in labels]
  ts = [M.shape[0] for M in mats]
  d_ranges = [
      #[3, 4, 5, 7],
      #[6, 8, 11, 13],
      #[6, 8, 11, 13],
      #[6, 8, 11, 13],
      #[19, 17, 13, 10,],
      #[19, 17, 13, 10,],
      #[19, 17, 13, 10,],
      #[15, 20, 25, 30],
      #[ 30, 25, 20, 15 ],
      #[ 7, 9, 10, 12, 15],
      #[10, 15, 20, 25],
      #[ 10, 12, 15, 17, 20,], 
      #[15, ]
      [ 4, 5, 6, 7]
      ] #list(range(1, 4))
  #d_ranges = [[ 15 ] for t in ts] #list(range(1, 4))
  #d_ranges = [ list(range(1, (t // 3) + 1)) for t in ts ] 
  #d_ranges = [list(range(1, 6)) for label in labels]

  num_expts = 100
  #algos = ['COMP', 'SBL', 'combined_COMP_NNOMP_random_cv',
  #    'combined_COMP_l1ls_cv']
  algos = ['COMP', 'combined_COMP_SBL', 'combined_COMP_SBL_clustered',
          'precise_SBL_combined_COMP_SBL', 'precise_SBL_COMP'] #'combined_COMP_NNOMP_random_cv']
  #algos = ['COMP', 'combined_COMP_SBL_clustered']
  #algos = ['combined_COMP_SBL_majority', 'combined_COMP_SBL_clustered']

#  algos = [
#      "combined_COMP_SBL_clustered",
#      "combined_COMP_SBL_majority_10_0.3",
#      "combined_COMP_SBL_majority_10_0.5",
#      "combined_COMP_SBL_majority_10_0.7",
#      "combined_COMP_SBL_majority_100_0.3",
#      "combined_COMP_SBL_majority_100_0.5",
#      "combined_COMP_SBL_majority_100_0.7",
#      "combined_COMP_SBL_union_10",
#      "combined_COMP_SBL_intersection_10",
#      "combined_COMP_SBL_union_100",
#      "combined_COMP_SBL_intersection_100",
#      ]

  save = True
  run_many_parallel_expts_many_matrices(mats, labels, d_ranges, algos,
      num_expts, save)

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



