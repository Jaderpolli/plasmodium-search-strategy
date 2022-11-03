using DelimitedFiles
using StatsBase
using CSV
using DataFrames
using CurveFit
using Clustering

include("../src/hurst.jl")

function hurst_individual()
    mkpath("examples/results-hurst/all")
    types = ["NINV","INV"]
    movies = ["20180206","20180216", "20180801", "20180813"]
    INVmoviesFirstParasite = [1, 10, 18, 25, 28]
    NINVmoviesFirstParasite = [1, 8, 27, 69, 88]
    for i in 1:length(types)
        type = types[i]
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
                if length(datax) < 22
                    writedlm("examples/results-hurst/all/hurst_$(type)$(k)_$(movie).csv", [0])
                else
                    max_t_window = round(Int64,length(datax)/4)
                    t_step = 4
                    hurstk = Hurst.hurst(datax, datay, max_t_window, t_step)
                    writedlm("examples/results-hurst/all/hurst_$(type)$(k)_$(movie).csv", hurstk, ',', header = true)
                end
            end
        end
    end
end

function hurst_mean()
    mkpath("examples/results-hurst/mean")
    types = ["NINV","INV"]
    movies = ["20180206","20180216", "20180801", "20180813"]
    INVmoviesFirstParasite = [1, 10, 18, 25, 28]
    NINVmoviesFirstParasite = [1, 8, 27, 69, 88]
    max_step = 25
    t_step = 4.0
    H = ["H" "t"]
    for t in 1:25
        ht = []
        for i in 1:length(types)
            type = types[i]
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
                    data = CSV.read("examples/results-hurst/all/hurst_$(type)$(k)_$(movie).csv", DataFrame)
                    data = Array(data[:,:])
                    if t > length(data[:,1])
                        continue
                    else
                        h = data[t,1]
                        if h ≤ 0
                            continue
                        else
                            ht = vcat(ht, h)
                        end
                    end
                end
            end
        end
        h_mean = mean(ht)
        h_mean = [h_mean t*t_step]
        H = vcat(H, h_mean)
    end
    writedlm("examples/results-hurst/mean/hurst_mean_all.csv", H, header = true, ',')
end

function hurst_mean_type()
    mkpath("examples/results-hurst/mean-type")
    types = ["NINV","INV"]
    movies = ["20180206","20180216", "20180801", "20180813"]
    INVmoviesFirstParasite = [1, 10, 18, 25, 28]
    NINVmoviesFirstParasite = [1, 8, 27, 69, 88]
    max_step = 25
    t_step = 4.0
    for i in 1:length(types)
        H_type =  ["H" "t"]
        type = types[i]
        for t in 1:25
            ht = []
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
                    data = CSV.read("examples/results-hurst/all/hurst_$(type)$(k)_$(movie).csv", DataFrame)
                    data = Array(data[:,:])
                    if t > length(data[:,1])
                        continue
                    else
                        h = data[t,1]
                        if h < 0
                            continue
                        else
                            ht = vcat(ht, h)
                        end
                    end
                end
            end
            h_mean_type = mean(ht)
            h_mean_type = [h_mean_type t*t_step]
            H_type = vcat(H_type, h_mean_type)
        end
        writedlm("examples/results-hurst/mean-type/hurst_mean_type_$(type).csv", H_type, header = true, ',')
    end
end

hurst_individual()
hurst_mean()
hurst_mean_type()