using CSV, DataFrames

function final(type)
    data = CSV.read("acf-julia/$(type)_vang.csv", DataFrame)
    data = select(data, :ID, :x, :y, :t, :v, :ang, :State)
    moviesINV = ["20180206","20180216", "20180801", "20180813"]
    INVmoviesFirstParasite = [1, 10, 18, 25, 28]
    moviesNINV = ["20180206","20180216", "20180801", "20180813", "20221220", "20230116", "20230111"]
    NINVmoviesFirstParasite = [1, 8, 27, 69, 88, 109, 139, 162]
    data2 = groupby(data, :ID)
    dates = []
    if type == "NINV"
        msdk = CSV.read("kmeans-julia/expected-k-means-results/k-means_NINV.csv", DataFrame)
        for j in 1:length(moviesNINV)
            datesj = repeat([moviesNINV[j]], NINVmoviesFirstParasite[j+1]-NINVmoviesFirstParasite[j])
            dates = vcat(dates, datesj)
        end
        for i in 1:length(data2)
            n = string(i, pad = 3)
            data2[i][!,:ID] = repeat(["NINV"*n], length(data2[i][!,:ID]))
            msd = repeat([msdk[i,:α]], length(data2[i][!,:ID]))
            k = repeat([msdk[i,:Cluster]], length(data2[i][!,:ID]))
            datesk = repeat([dates[i]], length(data2[i][!,:ID]))
            data2[i][!, :α] = msd
            data2[i][!, :kmeans] = k
            data2[i][!, :date] = datesk
        end
    else
        msdk = CSV.read("kmeans-julia/expected-k-means-results/k-means_INV.csv", DataFrame)
        for j in 1:length(moviesINV)
            datesj = repeat([moviesINV[j]], INVmoviesFirstParasite[j+1]-INVmoviesFirstParasite[j])
            dates = vcat(dates, datesj)
        end
        for i in 1:length(data2)
            n = string(i, pad = 2)
            data2[i][!,:ID] = repeat(["INV"*n], length(data2[i][!,:ID]))
            msd = repeat([msdk[i,:α]], length(data2[i][!,:ID]))
            k = repeat([msdk[i,:Cluster]], length(data2[i][!,:ID]))
            datesk = repeat([dates[i]], length(data2[i][!,:ID]))
            data2[i][!, :α] = msd
            data2[i][!, :kmeans] = k
            data2[i][!, :date] = datesk
        end
    end
    data = combine(data2, :ID, :date, :x, :y, :t, :v, :ang, :α, :kmeans, :State)
    #rename!(data1, :x1 => :ID)
    CSV.write("final-data-$(type).csv", data)
    kmeanfile = DataFrame()
    IDi = []
    alpha = []
    kmean = []
    date = []
    for j in 1:length(data2)
        IDi = vcat(IDi, data2[j][1, :ID])
        date = vcat(date, data2[j][1, :date])
        alpha = vcat(alpha, data2[j][1,:α])
        kmean = vcat(kmean, data2[j][1,:kmeans])
    end
    kmeanfile = DataFrame(ID = IDi, date = date, α = alpha, kmeans = kmean)
    CSV.write("final-kmean-$(type).csv", kmeanfile)
end

types = ["NINV", "INV"]
for type in types
    final(type)
end