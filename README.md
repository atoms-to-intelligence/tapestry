# Compressed Sensing Algorithms for Pooled Sampling and Testing using qPCR 
Code for the algorithms used in the following paper:
https://www.medrxiv.org/content/10.1101/2020.04.23.20077727v2

<!-- vim-markdown-toc GFM -->

* [Installation and Setup](#installation-and-setup)
* [Running Synthetic Experiments](#running-synthetic-experiments)
* [Adding Sensing Matrices](#adding-sensing-matrices)
  * [Deployed Matrices](#deployed-matrices)
* [Adding Algorithms](#adding-algorithms)
  * [Detailed instructions for adding algorithms](#detailed-instructions-for-adding-algorithms)
* [Running Algorithms on Lab Experiments](#running-algorithms-on-lab-experiments)

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

# Running Synthetic Experiments

Synthetic experiments are run by generating 1000 'x' values from our noise
model. See details of the noise model in the paper or 

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
For each matrix, all `d`

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

## Detailed instructions for adding algorithms

Detailed instructions for adding new algorithms can be found in
`algos/__init__.py` and `algos/zeros.py`

# Running Algorithms on Lab Experiments

The script `tools/run_ncbs_harvard_data.py` runs the algorithms `COMP`, `SBL`, and
`combined_COMP_NNOMP_random_cv` on these datasets. Please modify this script
and add any algorithms you may want to run.

<!---
![Python application](https://github.com/ssgosh/group-test/workflows/Python%20application/badge.svg)
-->

