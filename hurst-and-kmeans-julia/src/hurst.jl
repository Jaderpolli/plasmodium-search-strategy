module Hurst

using DelimitedFiles
using StatsBase
using CSV
using DataFrames
using CurveFit
using Clustering

    function hurst(x, y, max_t_window, t_step)
        msdt = []
        ts = []
        for t in 1:max_t_window
            dx = (x[1+t:end] .- x[1:end-t]).^2
            dy = (y[1+t:end] .- y[1:end-t]).^2
            sd = dx .+ dy
            meansd = mean(sd)
            msdt = vcat(msdt, meansd)
            ts = vcat(ts, t*t_step)
        end
        data = hcat(ts, msdt)
        log_msd = log10.(data[:,2])
        log_t = log10.(data[:,1])
        hurst = ["H" "t"]
        for l in 1:length(msdt)-1
            y = log_msd[l:l+1]
            x = log_t[l:l+1]
            h = 0.5*(y[2]-y[1])/(x[2]-x[1])
            ht = [h l*t_step]
            hurst = vcat(hurst, ht)
        end
        return(hurst)
    end
end