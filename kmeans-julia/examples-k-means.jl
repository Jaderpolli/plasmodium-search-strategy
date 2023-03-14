using DelimitedFiles
using StatsBase
using CSV
using DataFrames
using CurveFit
using Clustering
using LazyStack
using Plots

include("k-means.jl")

function kmeans_example()
    #mkpath("msd-results")
    #mkpath("k-means-results")
    types = ["INV" "NINV"]
    for type in types
        αs = Float64[]
        datas = readdir("kmeans-julia/raw-data/$(type)")
        j = 0
        for data in datas
            d = CSV.read("kmeans-julia/raw-data/$(type)/$(data)", DataFrame)
            datax = Array(d[!,:x])
            datay = Array(d[!,:y])
            max_t_window = round(Int64,length(datax)/2)
            t_step = 4
            msdt = KMeans.meansquare(datax, datay, max_t_window, t_step)
            h = KMeans.hurst(msdt, t_step)
            j += 1
            writedlm("kmeans-julia/expected-msd-results/msd_$(data)", msdt, ',', header = true)
            α = KMeans.alpha(datax, datay, max_t_window, t_step)
            αs = vcat(αs,α)
        end
        R = kmeans(transpose(αs), 3)
        data = ["ID" "α" "Cluster"]
        IDS = datas
        data1 = [IDS αs R.assignments]
        data = vcat(data, data1)
        writedlm("kmeans-julia/expected-k-means-results/k-means_$(type).csv", data, ',', header = true)
    end
end

kmeans_example()

function juntar_msd()
    types = ["INV", "NINV"]
    for type in types
        datas = readdir("kmeans-julia/expected-msd-results/$(type)")
        dados = DataFrame(ID = [], t = [], MSD = [])
        for data in datas
            msd = CSV.read("kmeans-julia/expected-msd-results/$(type)/$(data)", DataFrame, header = false)
            id = repeat([chop(data, head = 4, tail = 4)], outer = length(msd.Column1))
            dat = hcat(id, msd)
            rename!(dat, :x1 => :ID, :Column1 => :t, :Column2 => :MSD)
            dados = vcat(dados, dat)
        end
        CSV.write("kmeans-julia/expected-msd-results/$(type).csv", dados)
    end
end

#juntar_msd()

function juntar_cluster()
    types = ["INV", "NINV"]
    for type in types
        data = CSV.read("kmeans-julia/expected-msd-results/$(type).csv", DataFrame)
        kmean = CSV.read("kmeans-julia/expected-k-means-results/k-means_$(type).csv", DataFrame)
        datai = groupby(data, :ID)
        for i in 1:length(datai)
            datai[i][!, :Cluster] = repeat([kmean.Cluster[i]], outer = length(datai[i].t))
        end
        final = combine(datai, :ID, :t, :MSD, :Cluster)
        CSV.write("kmeans-julia/expected-msd-results/$(type)_cluster.csv", final)
    end
end

#juntar_cluster()

function plot_msd()
    types = ["INV", "NINV"]
    for type in types
        data = CSV.read("kmeans-julia/expected-msd-results/$(type)_cluster.csv", DataFrame)
        datai = groupby(data, :Cluster)
        colors = [:firebrick4, :firebrick2, :darkorange]
        cs = ["c1", "c2", "c3"]
        for i in 1:length(datai)
            df = DataFrame(datai[i])
            plt = plot()
            dataii = groupby(df, :ID)
            for j in 1:length(dataii)
                plt = plot!(dataii[j].t, dataii[j].MSD, lc = colors[datai[i].Cluster], label = false, lw = 0.65, yscale = :log10, xscale = :log10, msw = 0, title = "$(type)",
                xlabel = "Time window (s)", ylabel = "MSD (μm²)", frame = :box, size = (300, 300), xrange = (1, 50000), yrange = (1,50000))
            end
            plt = plot!([0.1;0.1],[100001;100002], lc = colors[i], lw = 3,  label = cs[i], fg_legend = false)
            savefig(plt, "kmeans-julia/expected-msd-results/$(type)_cluster_$(i).pdf")
        end        
    end
end

#plot_msd()

#hurst_example()