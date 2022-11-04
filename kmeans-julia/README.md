# Running Hurst and k-means in Julia

The present README file contains the instructions to run the codes in this folder, to obtain the *k-means clustering*.

## Folders and organization

This folder contains the following subfolders:

- raw-data: containing the raw-data used for analysis;
- src: containing the Julia modified modules for K-means;
- examples: containing the codes for analysing the raw-data/ files. In this subfolder, two other folders are created when running the codes inside:
    - results-k-means;
    - results-msd;

Now, let's see how to run everything and reproduce the results.

## Julia and package instalation

All codes are written in [Julia Language](https://julialang.org/).

Hence, in order to run everything, you must [Download](https://julialang.org/downloads/) and install properly the current stable release (in our analysis, it was used v1.7.2 -- see [older releases](https://julialang.org/downloads/oldreleases/).

After installing Julia (and add to `PATH`), you can open Julia simply by typing `julia` in the terminal.
Then, the following packages must be installed:

- DelimitedFiles
- StatsBase
- CSV
- DataFrames
- CurveFit
- Clustering

To install a package, simply run `import Pkg; Pkg.add("PackageName")`

## Running k-means

To reproduce the k-means clustering, open Julia in this main folder, then run

`include("examples-k-means.jl")`

A subfolder `results-msd` will be created, with the mean squared displacement as a function of time for each plasmodium used to obtain the α exponent to be clustered.

Finally, the subfolder `results-k-means` will be created with the results of the clustering inside for invadors and non-invadors.