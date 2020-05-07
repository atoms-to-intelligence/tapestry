# Compressed Sensing Algorithms for Pooled Sampling and Testing using qPCR 
Code for the algorithms used in the following paper:
https://www.medrxiv.org/content/10.1101/2020.04.23.20077727v2

<!-- vim-markdown-toc GFM -->

* [Installation and Setup](#installation-and-setup)
* [Code Layout](#code-layout)
* [Running Synthetic Experiments](#running-synthetic-experiments)
  * [Data Model](#data-model)
* [Adding Sensing Matrices](#adding-sensing-matrices)
  * [Deployed Matrices](#deployed-matrices)
* [Adding Algorithms](#adding-algorithms)
  * [Available algorithms](#available-algorithms)
  * [Detailed instructions for adding algorithms](#detailed-instructions-for-adding-algorithms)
* [Running Algorithms on Lab Experiments](#running-algorithms-on-lab-experiments)
  * [Experimental data location](#experimental-data-location)
* [Advanced Behaviour / Details](#advanced-behaviour--details)
  * [Statistics](#statistics)

<!-- vim-markdown-toc -->

# Installation and Setup

Clone this repository with

```
git clone https://github.com/ssgosh/group-test
cd group-test
```

Python 3.6 or above must be installed on your system. The python libraries
neeeded for running the scripts in this repo are listed in `requirements.txt`.
Install dependencies using:

```
pip3 install -r requirements.txt
```
# Code Layout

Code is layed out in subdirectories like this:

* `core/`: The core functionality, such as loading matrices, data model, performing
experiments etc are exposed as libraries.
* `tools/`: Executable scripts which use this functionality. These are the
ones you should run.
* `algos/`: pluggable algorithms can be added to this folder and be used.
* `inbuilt_algos/`: algorithms already written are in this folder.
* `utils/`: various utility methods and helpers.

# Running Synthetic Experiments

Synthetic experiments are run by generating 1000 'x' values from our data
model.

Experiments are run with the script `tools/run_expts.py`, with an invocation
of the function `run_stats_for_these_matrices()`. An example invocation of
this is:

```python
  run_stats_for_these_matrices(
      [
      "optimized_M_45_105_kirkman",
      "optimized_M_93_961_kirkman"],
      save=True
    )
```

Here `"optimized_M_45_105_kirkman"` and `"optimized_M_93_961_kirkman"` are
labels of matrices, of size 45x105 and 93x961 respectively. This function runs many
experiments for each of the matrices. We look deeper into the implementation
of `run_stats_for_these_matrices()`:

```python
def run_stats_for_these_matrices(labels, save):
  mats = [MDict[label] for label in labels]

  ts = [M.shape[0] for M in mats]

  d_ranges = [[5, 8, 12, 15, 17] for t in ts] #list(range(1, 4))
  #d_ranges = [ list(range(1, 16)) + [20, 25, 30, 35, 40] for item  in labels]
  #d_ranges = [ list(range(1, (t // 3) + 1)) for t in ts ] 
  #d_ranges = [list(range(1, 6)) for label in labels]

  num_expts = 1

  #algos = ['COMP', 'SBL', 'combined_COMP_NNOMP_random_cv',
  #    'combined_COMP_l1ls_cv']
  algos = ['SBL', 'ZEROS', 'l1ls', 'l1ls_cv']

  run_many_parallel_expts_many_matrices(mats, labels, d_ranges, algos,
      num_expts, save)

```
This function first finds the numpy matrices corresponding to the given matrix labels
(from MDict) and puts them in the list `mat`. It then finds the number of
rows `t` for each matrix. The number of infections for which each matrix
must be run is given in 

## Data Model

The data model is following. Say there are 'n' samples, out of which 'd' are
infected. Then we create a 0/1 vector of infections such that each choice of
'd' infections is equally likely. This is found in
`core/comp.py::create_infection_array_with_num_cases()`. Let's call this `bool_x`

'x', the vector of viral loads in each sample, is generated from `bool_x`. The
positive values of 'x' are drawn from the range `[1/32768, 1]`
uniformly at random. Notice that this 1) has high dynamic range, spanning many
orders of magnitude as we may expect in real viral loads 2) is bounded away
from 0.

Once we have 'x', then given a t x n matrix A, we may generate 'y', the vector
of test results. In the noiseless case, `y = Ax`. In the noisy case, `y =
Ax(1+p)^eps`. Here p = 0.95 and eps = vector of independent Gaussians with mean
0 and standard deviation 0.1. This model comes from the PCR amplification
process, where we can only observe cycle times for a given threshold
fluorescence value. The observed cycle time itself has some variability, which
eps is modelling.  Essentially, log y is Gaussian. This is found in
`core/cs.py::get_quantitative_results()`
                                                                             

# Adding Sensing Matrices

Matrices are present in the `mats/` subfolder in the form of txt files which are
loaded at startup by `core/matrices.py`. You may drop your own matrices in
`mats/extra/` folder, which will be automatically loaded and assigned a
matrix label, which is the name of the matrix file minus the extension. Matrix
txt files must have one line per row, space-separated values and no other
delimiters of any kind. This is the default format used by `numpy.loadtxt` and
`numpy.savetxt`.

All matrices are added to the matrices dictionary `MDict` in `core/matrices.py`.

## Deployed Matrices

# Adding Algorithms

New algorithms can be easily added to the `algos/` folder. Please edit
`algos_dict` to add your algorithm name and register the function
corresponding function to be called.

## Available algorithms

Inbuilt algorithms are in `inbuilt/` folder. The following algorithms are
available: `'COMP', 'SBL', 'l1ls', 'l1ls_cv', 'NNOMP', 'NNOMP_random_cv'`. You
may also run algorithms like `'combined_COMP_l1ls', 'combined_COMP_NNOMP'` etc,
which first run COMP and filter out the columns corresponding to negative
samples (as well as rows which correspond to negative tests).

Following inbuilt algorithms are deprecated: 
  `'lasso', 'lasso_cv', 'ista', 'fista'`

## Detailed instructions for adding algorithms

Detailed instructions for adding new algorithms can be found in
`algos/__init__.py` and `algos/zeros.py`

# Running Algorithms on Lab Experiments

The script `tools/run_ncbs_harvard_data.py` runs the algorithms `COMP`, `SBL`, and
`combined_COMP_NNOMP_random_cv` on these datasets. Please modify this script
and add any algorithms you may want to run.

The script `tools/experimental_data_manager.py` contains methods which return
data from experiments

## Experimental data location

# Advanced Behaviour / Details

## Statistics

Stat for each individual experiment may be saved by setting `save=True` when invoking 
`run_many_parallel_expts_many_matrices() or run_stats_for_these_matrices()`.
These are saved in a directory structure under the folder `expts_stats/`. This
folder *must not be added* to the repo, since its contents can be of the order
of several hundred MegaBytes. 

The directory structure for this is of the form
`expt_stats/[matrix_label]/[algo]/[d]/expt_stats.p.gz`. `expt_stats.p.gz`
contains a list containing the stats for each individual experiment. 

These stats can be used to perform post-facto analysis, such as finding
confidence intervals for various stats or finding why a particular algorithm fails on a
particular 'x'.

The script `tools/stats_tools.py` finds confidence intervals using these
pickled experiments.

<!---
![Python application](https://github.com/ssgosh/group-test/workflows/Python%20application/badge.svg)
-->

