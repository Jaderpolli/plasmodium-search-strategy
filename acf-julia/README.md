# Running the codes

To run the codes to reproduce the autocorrelation function for the $v \cos (\theta_{t,t-1})$ time series, the following Julia Packages must be installed:

- `DataFrames`
- `CSV`
- `Plots`
- `LinearAlgebra`
- `StatsBase`
- `DelimitedFiles`
- `Measurements`

Then the code to process the data and transform it into a collection of stretches of each state with information about the turning angles and velocity is `data_processing.jl` and to run it, just open this folder on Julia and type

`include("data_processing.jl")`

To obtain the results of autocorrelation function and Figure 1 of Supplementary material, just run the `acf_angles.jl` code by

`include("acf_angles.jl")`