using DelimitedFiles
using StatsBase
using CSV
using DataFrames
using CurveFit
using Clustering

include("k-means.jl")

function kmeans_example()
    mkpath("results-k-means")
    mkpath("results-msd")
    types = ["INV" "NINV"]
    movies = ["20180206","20180216", "20180801", "20180813"]
    INVmoviesFirstParasite = [1, 10, 18, 25, 28]
    NINVmoviesFirstParasite = [1, 8, 27, 69, 88]
    l = 1
    for i in 1:length(types)
        type = types[i]
        αs = Float64[]
        for j in 1:length(movies)
            movie = movies[j]
            if type == "INV"
                FirstParasite = INVmoviesFirstParasite[j]
                LastParasite = INVmoviesFirstParasite[j+1]-1
            elseif type == "NINV"
                FirstParasite = NINVmoviesFirstParasite[j]
                LastParasite = NINVmoviesFirstParasite[j+1]-1
            end
            for k in FirstParasite:LastParasite
                data = CSV.read("raw-data/position_$(type)$(k)_$(movie).csv", DataFrame)
                datax = Array(data[:,1])
                datay = Array(data[:,2])
                max_t_window = round(Int64,length(datax)/2)
                t_step = 4
                msdt = KMeans.meansquare(datax, datay, max_t_window, t_step)
                writedlm("results-msd/msd_$(type)$(k)_$(movie).csv", msdt, ',', header = true)
                α = KMeans.alpha(datax, datay, max_t_window, t_step)
                αs = vcat(αs,α)
                l += 1
            end
        end
        R = kmeans(transpose(αs), 3)
        data = ["α" "Cluster"]
        data1 = hcat(αs, R.assignments)
        data = vcat(data, data1)
        writedlm("results-k-means/k-means_$(type).csv", data, ',', header = true)
    end
end

kmeans_example()
