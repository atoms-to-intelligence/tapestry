# Compressed Sensing Algorithms for Pooled Sampling and Testing using qPCR 
Code for the algorithms used in the following paper:
https://www.medrxiv.org/content/10.1101/2020.04.23.20077727v2

<!-- vim-markdown-toc GFM -->

* [Installation and Setup](#installation-and-setup)
* [Code Layout](#code-layout)
* [Running Synthetic Experiments](#running-synthetic-experiments)
  * [Output Statistics](#output-statistics)
  * [Saving the output](#saving-the-output)
  * [Data Model](#data-model)
* [Adding Sensing Matrices](#adding-sensing-matrices)
  * [Deployed Matrices](#deployed-matrices)
* [Adding Algorithms](#adding-algorithms)
  * [Adding an algorithm file to `algos/` folder](#adding-an-algorithm-file-to-algos-folder)
  * [Adding Algorithm to `algo_dict`](#adding-algorithm-to-algo_dict)
  * [Unit Testing your new algorithm](#unit-testing-your-new-algorithm)
  * [Running synthetic expts with new algorithm](#running-synthetic-expts-with-new-algorithm)
  * [Inbuilt algorithms](#inbuilt-algorithms)
  * [Detailed instructions for adding algorithms](#detailed-instructions-for-adding-algorithms)
  * [Combining With COMP or Definite Defects](#combining-with-comp-or-definite-defects)
  * [Performing Cross-validation for your algorithm](#performing-cross-validation-for-your-algorithm)
* [Running Algorithms on Lab Experiments](#running-algorithms-on-lab-experiments)
  * [Experimental data location](#experimental-data-location)
* [Running Tests](#running-tests)
* [Advanced Behaviour / Details](#advanced-behaviour--details)
  * [Detailed Statistics](#detailed-statistics)
  * [Other directory Layout](#other-directory-layout)
  * [Core code layout](#core-code-layout)

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
* `matrix_gen/`: code for generating sensing matrices.
* `test/`: Code which runs tests
* `matlab/`: some algorithm implementation in matlab. These are not actively maintained.

Data and other files are layed out like this:

* `mats/`: Matrices in ".txt" format.
  * `mats/extra`: Extra matrices automatically loaded.
  * `mats/kirman`: Kirkman matrices
  * `mats/sts`: STS matrices
  * `mats/unparsed`: 
* `stats/`: Various text files containing statistics for various runs.
* `data/`: Contains data from wet lab experiments - either done by us or
  elsewhere.
* `expt_stats/`: This is not part of the repo, but may be automatically
  created later. See [Detailed Statistics](#detailed-statistics) section.

# Running Synthetic Experiments

Synthetic experiments are run by generating 1000 'x' values from our data
model.

Experiments are run with the script `tools/run_expts.py` as

```
python3 tools/run_expts.py
```

Be careful to stay on the top-level directory and invoke the script from
there - else the script may not be able to find other modules, matrices etc.

This script invokes function `run_stats_for_these_matrices()`. An
example invocation of this function is:

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

  num_expts = 1000

  #algos = ['COMP', 'SBL', 'combined_COMP_NNOMP_random_cv',
  #    'combined_COMP_l1ls_cv']
  algos = ['SBL', 'ZEROS', 'l1ls', 'l1ls_cv']

  run_many_parallel_expts_many_matrices(mats, labels, d_ranges, algos,
      num_expts, save)

```

* `mat`: numpy matrices corresponding to the matrix labels.
* `algos`: which algorithms to run 
* `d_ranges`: list of lists. Specifies for each matrix which 'd' values to run experiments for.
* `num_expts`: (=1000 by default) These many independent experiments are run for each 
  (matrix, algo, d) triplet.
* `ts`: list of number of rows in each matrix. Useful for specifying 'd' in
  some cases

These are then passed to the method `run_many_parallel_expts_many_matrices()`.
The `save` flag saves each experiment. See [Detailed Statistics](#detailed-statistics)
for detail.

For each experiment, 'x' (viral load vector for each sample) is sampled from
our [data model](#data-model), `y = Ax` is computed (A being the matrix), and
noise is added to y according to our noise model.  Only `y` and `A` are passed
to the algorithm. The result `x_est`, the estimated `x`, is compared with the
ground-truth `x` and statistics are computed as below.

## Output Statistics

For each matrix label and algoithm pair, one table is printed listing various
stats. One such table is shown below:

```
t = 21, n = 70, matrix = optimized_M_21_70_kirkman

COMP
	d	Precision	Recall (Sensitivity) 	Specificity	surep	unsurep	false_pos
	1	1.000			1.000		1.000		 1.0	  0.0	    0.0
	2	1.000			1.000		1.000		 2.0	  0.0	    0.0
	3	0.760			1.000		0.986		 2.7	  1.2	    0.9
	4	0.571			1.000		0.955		 1.7	  5.3	    3.0
	5	0.454			1.000		0.908		 0.6	 10.4	    6.0
	6	0.376			1.000		0.844		 0.1	 15.9	   10.0
	7	0.334			1.000		0.778		 0.0	 21.0	   14.0
```

Each row of the above table corresponds to a different value of 'd', the
number of infections (1- 7 in this case) out of 'n' (70). For each 'd', 1000
expts have been done. From those 1000, aggregate statistics such as Precision,
Recall (Sensitivity or True Poitive Rate) and Specificity (or True Negative
Rate) are computed. `false_pos` is the average number of false positives per
expt. `surep` is the average number of positives we are sure about, and rest of
the predicted positives are classified as `unsurep`. `surep` is computed using
the Definite Defects algorithm.

If you want to modify the columns which are printed in the above table, see
`core/cs_expts.py::print_expts()`. Many more statistics are already computed.
See `__init__()`, `print_stats()` and `return_stats()` methods of class `CSExpts` defined in 
`core/cs_expts.py`.

## Saving the output

Typically, the output tables may be saved in a text file using redirection.
e.g.:

```bash
python3 tools/run_expts.py > output_MY_MATRIX_MY_ALG.txt
```

Stats for individual experiments may be saved as in [Detailed Statistics](#detailed-statistics).

Saving the aggregate stats in a dictionary for later perusal also desirable
but not currently implemented. 

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
`core/cs.py::CS::get_quantitative_results()`.

# Adding Sensing Matrices

Matrices are present in the `mats/` subfolder in the form of txt files which are
loaded at startup by `core/matrices.py`. *You may copy your own matrices in
`mats/extra/` folder*, which will be automatically loaded and assigned a
matrix label, which is the name of the matrix file minus the extension. Matrix
txt files must have one line per row, space-separated values and no other
delimiters of any kind. This is the default format used by `numpy.loadtxt()` and
`numpy.savetxt()`.

All matrices are added to the matrices dictionary `MDict` in `core/matrices.py`.
Matrix labels are also available as python variable names in `tools/run_expts.py`
and other relevant places.

Matrix labels for matrices `mats/extra/` are also found in the variable
`extra_mlabels`, and the dictionary `extra_mats_dict`. Kirkman matrices are found
in `kirkman_mlabels`, while sts matrices are found in `sts_mlabels`. 

STS matrices are generated using the Bose construction. Code is found in
`matrix_gen/sts.py`.

All matrices have their `writeable` flag set to `False`, so that they are not
accidentally written to. If you want to modify a matrix at runtime, it is best to
work with a copy of that matrix. Doing something like `A = A / some_constant` is ok
since it creates a copy and reassigns the variable `A` to that copy.

## Deployed Matrices

Currently deployed matrices in our Android app BYOM Smart Testing are found in
`core/get_test_results.py` in the dictionary `MSizeToLabelDict`. All matrices ever
deployed form the keys of the dictionary `mat_codenames` in the same file.

# Adding Algorithms

New algorithms can be easily added to the `algos/` folder. Please edit `algos_dict`
to add your algorithm name and register the corresponding function to be called.

## Adding an algorithm file to `algos/` folder

For example, to add an algorithm called `MY_ALG`, add a python file called
`my_alg.py`, with the following content: 

```python
# core/config.py contains various configs like noise model, some global
# hyperparameters like tau
from core import config

def my_alg(params):
  A = params["A"]
  y = params["y"]

  A_bool = (A > 0).astype(np.int32)
  y_bool = (y > 0).astype(np.int32)

  # Your code goes here. it produces x_est, an n-length numpy array containing the
  # estimated viral loads
  res = {'x_est' : x_est}
  return res
```

The input `params` is a dictionary containing the matrix A and the observed viral load
vector `y`. Both `A` and `y` can contain float values, so if you need only 0/1
values, you should convert and use `A_bool` and `y_bool` as shown above.

If the algorithm you're using doesn't give exact 0's for the 0-entries of `x_est`,
you will need to convert it using.

```python
x_est[x_est < tau] = 0
```
where `tau` is an appropriately chosen threshold. A good rule of thumb for `tau` is
`0.01 * np.min(y/row_sums)`, where `row_sums = np.sum(A, axis=-1)`. 

## Adding Algorithm to `algo_dict`

Add your `MY_ALG` to algos/__init__.py as follows:

```python
# import your file here
from . import my_alg

# Add to the algo_dict
algo_dict = {
    "ZEROS" : zeros.zeros,
    "ZEROS_cv" : zeros.zeros_cv,
    "MY_ALG" : my_alg.my_alg, # "ALGO_NAME" : module_name.function_name
    }
```

## Unit Testing your new algorithm

It is recommended that for every algorithm you add a test file which tests that
algorithm in isolation. For `algos/my_alg.py`, add `algos/test_my_alg.py`. For
example:

```python
import sys
# Following hack is needed for importing core/matrices.py etc in this test file. If you
don't need to import anything from the top level directory then you may skip this.
sys.path.append('.')

import my_alg

# Following to be added if you want to use any matrices defined in mats/
from core.matrices import *

def test_my_alg():
  # test code goes here

if __name__ == '__main__':
  test_my_alg()
```

Please see `inbuilt_algos/test_nnompcv.py` and `inbuilt_algos/test_sbl.py` as
examples.

Run the file as:

```bash
python3 algos/test_my_alg.py
```

Also add an import guard in `my_alg.py`. This is so that it is never run directly. 

```python
if __name__ == '__main__':
  raise ValueError('Please run algos/test_my_alg.py. This is a library file')
```

## Running synthetic expts with new algorithm

Once the algorithm is added, modify `core/cs.py::run_stats_for_these_matrices()`

```python
def run_stats_for_these_matrices(labels, save):

  .......

  algos = ['MY_ALG']

  run_many_parallel_expts_many_matrices(mats, labels, d_ranges, algos,
      num_expts, save)

```

and run

```bash
python3 tools/run_expts.py
```

## Inbuilt algorithms

Inbuilt algorithms are in `inbuilt/` folder. The following algorithms are
available: `'COMP', 'SBL', 'l1ls', 'l1ls_cv', 'NNOMP', 'NNOMP_random_cv'`. You may
also run algorithms like `'combined_COMP_l1ls', 'combined_COMP_NNOMP'` etc, which
first run COMP and filter out the columns corresponding to negative samples 
(as well as rows which correspond to negative tests).

Following inbuilt algorithms are deprecated: 
  `'lasso', 'lasso_cv', 'ista', 'fista'`

## Detailed instructions for adding algorithms

Detailed instructions for adding new algorithms can also be found in
`algos/__init__.py` and `algos/zeros.py`. The latter implements a trivial algorithm
which returns all 0's as `x_est`.

## Combining With COMP or Definite Defects

If you have registered an algorithm (say `MY_ALG`) in `algos/__init__.py`,
then the core code automatially also recognizes an algorithm by the name
`combined_COMP_MY_ALG`. This algorithm performs preprocessing using the `COMP`
algorithm and removes negative tests and columns corresponding to samples
which are definitely negative according to COMP. Typically this reduces the
size of the problem drastically and better performance is observed.

The Definite Defects algorithm *is run on top of any algorithm by
default* and the final predictions are the union of the predictions made by the
internal algorithm and the definite defects algorithm. This behaviour can be
changed by commenting out the following line:

```python
infected = (infected + infected_dd > 0).astype(np.int32)
```
from `core/cs.py::CS::decode_lasso()`.

## Performing Cross-validation for your algorithm

It is a good idea to perform cross-validation for your algorithms to determine the
ideal hyperparameter settings given the data. Currently cross validation must be
done separately for each algorithm. You should add a separate function in
`my_alg.py` called `my_alg_cv()` which performs the cross-validation using subsets
of the data. A separate algo name called `MY_ALG_cv` must be added to `algo_dict`
for the cross-validation algo to be visible to rest of the code.

# Running Algorithms on Lab Experiments

The script `tools/run_ncbs_harvard_data.py` runs the algorithms `COMP`, `SBL`, and
`combined_COMP_NNOMP_random_cv` on these datasets. Please modify this script
and add any algorithms you may want to run.

The file `tools/experimental_data_manager.py` contains methods which return
data from experiments performed at Harvard. Any new experimental data will be added
to this file. It also contains a method to generate synthetic data using our noise
model.

## Experimental data location

Data is located in the `data/` folder, with subfolders indicating the source of the
data (such as `harvard`, `ncbs`). These are in various formats, and may contain raw
flourescence values for each well and cycle or cycle time thresholds for a given
flourescence. 

# Running Tests

`test/test_all.py` performs various sanity checks for deployed matrices and also runs
each deployed matrix on various algorithms. Do add your algorithm to file so that
tests are run on it.

Some smaller tests are also present in various files, such as `inbuilt_algos/test_sbl.py`,
`core/test_get_test_results.py` etc.

# Advanced Behaviour / Details

## Detailed Statistics

Stat for each individual experiment may be saved by setting `save=True` when invoking 
`run_many_parallel_expts_many_matrices() or run_stats_for_these_matrices()`.
These are saved in a directory structure under the folder `expts_stats/`. This
folder *must not be added* to the repo, since its contents can be of the order
of several hundred MegaBytes. Please use separate backup mechanisms such as
Google Drive, Dropbox or disk backup to save experiments. The files saved in
this directory come precompressed with gzip.

The directory structure for this is of the form
`expt_stats/[matrix_label]/[algo]/[d]/expt_stats.p.gz`. `expt_stats.p.gz`
contains a list containing the stats for each individual experiment. 

These stats can be used to perform post-facto analysis, such as finding
confidence intervals for various stats or finding why a particular algorithm fails on a
particular 'x'.

The script `tools/stats_tools.py` finds confidence intervals using these
pickled experiments via bootstrapping.

## Other directory Layout

## Core code layout

TBD

<!---
![Python application](https://github.com/ssgosh/group-test/workflows/Python%20application/badge.svg)
-->

