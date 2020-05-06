# Compressed Sensing Algorithms for Pooled Sampling and Testing using qPCR 
<!-- vim-markdown-toc GFM -->

* [Installation and Setup](#installation-and-setup)
* [Running Synthetic Experiments](#running-synthetic-experiments)
* [Adding Sensing Matrices](#adding-sensing-matrices)
  * [Deployed Matrices](#deployed-matrices)
* [Adding Algorithms](#adding-algorithms)
  * [Detailed instructions for adding algorithms](#detailed-instructions-for-adding-algorithms)
* [Running Algorithms on Lab Experiments](#running-algorithms-on-lab-experiments)

<!-- vim-markdown-toc -->

Code for the algorithms used in the following paper:
https://www.medrxiv.org/content/10.1101/2020.04.23.20077727v2

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
experiments for each of the matrices

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
