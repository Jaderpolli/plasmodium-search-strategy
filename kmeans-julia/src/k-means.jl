module KMeans

using DelimitedFiles
using StatsBase
using CSV
using DataFrames
using CurveFit
using Clustering

    function alpha(x, y, max_t_window, t_step)
        msdt = []
        for t in 1:max_t_window
            dx = (x[1:end-t] .- x[1+t:end]).^2
            dy = (y[1:end-t] .- y[1+t:end]).^2
            sd = dx .+ dy
            meansd = mean(sd)
            msdt = vcat(msdt, meansd)
        end
        t = range(t_step, t_step*max_t_window, length = max_t_window)
        α = power_fit(t, msdt)
        return(α[2])
    end

    function meansquare(x, y, max_t_window, t_step)
        msdt = []
        for t in 1:max_t_window
            dx = (x[1:end-t] .- x[1+t:end]).^2
            dy = (y[1:end-t] .- y[1+t:end]).^2
            sd = dx .+ dy
            meansd = mean(sd)
            msdt = vcat(msdt, meansd)
        end
        t = range(t_step, t_step*max_t_window, length = max_t_window)
        data = ["t" "MSD"]
        data1 = hcat(t, msdt)
        data = vcat(data, data1)
        return(data)
    end
end