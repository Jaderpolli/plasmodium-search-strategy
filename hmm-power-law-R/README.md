# Running R codes


## System requirements

All codes were run in Windows 10 OS, but with the fullfilment of the requirements bellow, one should be able to run it either on Linux or Mac too.

## Instalation guide

First, you need to install the last version of R language for your operational system from the [CRAN project web site](https://cran.r-project.org/).

Then, install the following packages:

- moveHMM
- ggplot2
- quantmod
- poweRlaw
- cowplot
- pracma
- dplyr

To install a package, run, into the RGui, or terminal, the following command

`install.packages("package_name")`

In the first run, you'll have to sellect a CRAN mirror for the section. Select the one that is closer to your location.

## Demo

To run and reproduce the results of the HMM, simply run the `hmm.r` code.
This will create the `invader-hmm-states-results.csv` and `non-invader-hmm-states-results.csv` archives, where the informations on the movement state performed by each sporozoite at each time step is recorded.

To run and reproduce the levy power law results, run the `levy-power-law.r` code.
This will generate the reproduction of figures 2(b-d).


## Reproduction instructions

In order to reproduce figure 1(b-e) you can merge data from results of HMM with the raw data on turning angles and speed and plot in a software of your preference.