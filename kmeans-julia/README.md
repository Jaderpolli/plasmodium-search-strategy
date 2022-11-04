# Running k-means in Julia

The present README file contains the instructions to run the codes in this folder, to obtain the *k-means clustering*.

## System requirements

The codes were tested in Windows 10, but, after the proper packages installation - see next section - one should be able to run it either on Linux or Mac.

All codes are written in [Julia Language](https://julialang.org/).

In order to run everything, you must [Download](https://julialang.org/downloads/) and install properly the current stable release (in our analysis, it was used v1.7.2 -- see [older releases](https://julialang.org/downloads/oldreleases/).

## Julia and package instalation

After installing Julia (and add to `PATH`), you can open Julia simply by typing `julia` in the terminal.
Then, the following packages must be installed:

- DelimitedFiles
- StatsBase
- CSV
- DataFrames
- CurveFit
- Clustering

To install a package, simply run `import Pkg; Pkg.add("PackageName")`.

The installation time must not exceed 5 to 10 minutes, depending on the internet conection and how many packages must be manually installed.

## Running k-means and expected output

To reproduce the k-means clustering, open Julia in this main folder, then run

`include("examples-k-means.jl")`

A subfolder `msd-results` will be created, with the mean squared displacement as a function of time for each plasmodium used to obtain the α exponent to be clustered.

Finally, the subfolder `k-means-results` will be created with the results of the clustering inside for invadors and non-invadors.

The expected outputs are already on `expected-msd-results` and `expected-k-means-results` folders.

The expected run time for the demo should not surpasse one minute of running in any modern computer with the correct packages and versions properly installed.