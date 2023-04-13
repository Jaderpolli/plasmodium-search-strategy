using StatsBase, DelimitedFiles, DataFrames, CSV, Plots, Plots.Measures, Measurements

function procedures(D)
    minlength = [30,30,100,100]
    types = ["inv1", "inv2", "ninv1", "ninv2"]
    j = 0
    for type in types
        j += 1
        dfs = CSV.read("acf-julia/$(type)_ang.csv", DataFrame)
        dropmissing!(dfs)
        dfs = DataFrame(ID = dfs.Column1, v = dfs.vang)
        dfs_ID = groupby(dfs,:ID)
        for i in 1:length(dfs_ID)
            data = Float64.(Array(dfs_ID[i].v))
            if length(data) < minlength[j]
                nothing
            else
                println(length(data), type)
                lags = (0:length(dfs_ID[i].v)-1)
                acov = autocov(data,lags)
                acor = autocor(data,lags)
                dfs_ID[i][!, :lags] = lags
                dfs_ID[i][!, :acov] = acov
                dfs_ID[i][!, :acor] = acor
            end
        end
        mkpath("acf-julia/acov")
        mkpath("acf-julia/acor")
        dfs = combine(dfs_ID, :ID, :v, :lags, :acov, :acor)
        CSV.write("acf-julia/acor/acor_$(type)_ang_D_$(4*D).csv", dfs)
    end
end

function plotacf(D)
    
    plt_ninv = plot(xlabel = "lag (s)", ylabel = "⟨C(lag)⟩ ± σ")
    dfs = CSV.read("acf-julia/acor/acor_ninv1_ang_D_4.csv", DataFrame)
    dfs = DataFrame(ID = dfs.ID, lags = dfs.lags, acor = dfs.acor)
    dfs_lags = groupby(dfs,:lags)
    acfs = zeros(length(dfs_lags))
    stds = zeros(length(dfs_lags))
    for i in 1:length(dfs_lags)
        data = dfs_lags[i].acor
        acfs[i] = mean(skipmissing(data))
        stds[i] = std(skipmissing(data))
    end
    lag = range(0, length(acfs)-1)
    acf_std = acfs[1:20] .± stds[1:20]
    writedlm("acf-julia/acf_std_ninv1.csv", acf_std)
    plt_ninv = plot!(4*D .* lag[1:20], acf_std, #marker = :circle,# yscale = :log10, xscale = :log10, fontfamily = "Helvetica",
                    linecolor = :blue, markerstrokecolor = :blue,  frame = :box, yticks = -0.25:0.25:1,
                    linealpha = 1/(2),# markeralpha = 1/2, #markerstrokewidth = 0, 
                    size = (500, 300), label = false, dpi = 400, foreground_color_legend = nothing,
                    left_margin = 5mm, bottom_margin = 2mm)
    plt_ninv1 = plot!([0;0], [10,11], linecolor = :blue, markercolor = :blue, label = "NINV St1", yrange = (-0.25,1.1), xrange = (-1,4*10+1))

    dfs = CSV.read("acf-julia/acor/acor_ninv2_ang_D_4.csv", DataFrame)
    dfs = DataFrame(ID = dfs.ID, lags = dfs.lags, acor = dfs.acor)
    dfs_lags = groupby(dfs,:lags)
    acfs = zeros(length(dfs_lags))
    stds = zeros(length(dfs_lags))
    for i in 1:length(dfs_lags)
        data = dfs_lags[i].acor
        acfs[i] = mean(skipmissing(data))
        stds[i] = std(skipmissing(data))
    end
    lag = range(0, length(acfs)-1)
    acf_std = acfs[1:20] .± stds[1:20]
    writedlm("acf-julia/acf_std_ninv2.csv", acf_std)
    plt_ninv = plot!(4*D .* lag[1:20], acf_std, #marker = :circle,# yscale = :log10, xscale = :log10, fontfamily = "Helvetica",
                    linecolor = :red, markerstrokecolor = :red,  frame = :box,
                    linealpha = 1/(2),# markeralpha = 1/(2*k), #markerstrokewidth = 0, 
                    size = (500, 300), label = false, dpi = 400, foreground_color_legend = nothing,
                    left_margin = 5mm, bottom_margin = 2mm)
    plt_ninv = plot!([0;0], [10,11], linecolor = :red, annotationfontsize = 10, 
                markercolor = :red, label = "NINV St2", yrange = (-0.25,1.1), xrange = (-1,4*10+1), top_margin = -5mm, annotation = (20,0.95, ("b", "Helvetica Bold")))

    plt_inv = plot(ylabel = "⟨C(lag)⟩ ± σ")
    dfs = CSV.read("acf-julia/acor/acor_inv1_ang_D_4.csv", DataFrame)
    dfs = DataFrame(ID = dfs.ID, lags = dfs.lags, acor = dfs.acor)
    dfs_lags = groupby(dfs,:lags)
    acfs = zeros(length(dfs_lags))
    stds = zeros(length(dfs_lags))
    for i in 1:length(dfs_lags)
        data = dfs_lags[i].acor
        acfs[i] = mean(skipmissing(data))
        stds[i] = std(skipmissing(data))
    end
    lag = range(0, length(acfs)-1)
    acf_std = acfs[1:20] .± stds[1:20]
    writedlm("acf-julia/acf_std_inv1.csv", acf_std)
    plt_inv = plot!(4*D .* lag[1:20], acf_std, #marker = :circle,# yscale = :log10, xscale = :log10, fontfamily = "Helvetica",
                    linecolor = :blue, markerstrokecolor = :blue, frame = :box, yticks = -0.25:0.25:1,
                    linealpha = 1/(2),# markeralpha = 1/(2*k), #markerstrokewidth = 0, 
                    size = (500, 300), label = false, dpi = 400, foreground_color_legend = nothing,
                    left_margin = 5mm, bottom_margin = 2mm)
    plt_inv1 = plot!([0;0], [10,11], linecolor = :blue, label = "INV St1", yrange = (-0.5,1.1), xrange = (-1,4*10+2))


    dfs = CSV.read("acf-julia/acor/acor_inv2_ang_D_4.csv", DataFrame)
    dfs = DataFrame(ID = dfs.ID, lags = dfs.lags, acor = dfs.acor)
    dfs_lags = groupby(dfs,:lags)
    acfs = zeros(length(dfs_lags))
    stds = zeros(length(dfs_lags))
    for i in 1:length(dfs_lags)
        data = dfs_lags[i].acor
        acfs[i] = mean(skipmissing(data))
        stds[i] = std(skipmissing(data))
    end
    lag = range(0, length(acfs)-1)
    acf_std = acfs[1:20] .± stds[1:20]
    writedlm("acf-julia/acf_std_inv2.csv", acf_std)
    plt_inv = plot!(4*D .* lag[1:20], acf_std, #marker = :circle,# yscale = :log10, xscale = :log10, fontfamily = "Helvetica",
                    linecolor = :red, markerstrokecolor = :red,  frame = :box, 
                    linealpha = 1/(2),# markeralpha = 1/(2*k), #markerstrokewidth = 0, 
                    size = (500, 300), label = false, dpi = 400, foreground_color_legend = nothing,
                    left_margin = 5mm, bottom_margin = 2mm)
    plt_inv = plot!([0;0], [10,11], linecolor = :red, annotationfontsize = 10, 
                label = "INV St2", yrange = (-0.5,1.1), xrange = (-1,4*10+1),
                annotation = (20,0.95, ("a", "Helvetica Bold")),)

    plt = plot(plt_inv, plt_ninv, layout = (2,1), size = (400,300))

    savefig(plt, "acf-julia/acf_angle_D_$(4*D)_junto.pdf")    
    
end

function meanacf(D)
    types = ["ninv", "inv"]
    colors = ["blue", "red"]
    for type in types
        plt = plot(xlabel = "lag(s)", ylabel = "<C(lag)>")
        for k in reverse(1:2)
            typei = string(type,k)
            acors = readdlm("acf-julia/acor/acor_$(typei)_D_$(4*D).csv", header = false)
            acors = replace(acors, "" => missing, NaN => missing)
            meanac = zeros(length(acors[1,2:end]))
            for j in 2:length(acors[1,2:end])
                meanac[j-1] = mean(skipmissing(acors[:,j]))
            end
            lag = range(0, length(meanac)-1)
            plt = plot!(4*D .*lag, meanac, marker = :circle,# yscale = :log10, xscale = :log10,
                            linecolor = colors[k], markercolor = colors[k],  
                            linealpha = 1/(2), #markeralpha = 1/(2*k), 
                            markerstrokewidth = 0, size = (500, 300), label = false, dpi = 400,
                            left_margin = 5mm, bottom_margin = 2mm)
            plt = scatter!([10,10], markersize = 1, markerstrokewidth = 0, markercolor = colors[k], label = "$(typei)", yrange = (-0.3,0.3), title = "group each $(4*D)s")
        end
        savefig(plt, "acf-julia/mean-acf-$(type)_D_$(4*D).png")    
    end
end

procedures(1)
plotacf(1)

