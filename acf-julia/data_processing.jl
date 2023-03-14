using DataFrames, CSV, Plots, LinearAlgebra

# This function takes a dataframe where column :nameColumn is a list of integer numbers
# and returns another list of the discontinuities of the list. That is, those points where
# list[j] ≠ list[j-1]+1, allowing to obtain the positions where the list is continuous in the integer set.

function separateContinuous(df, nameColumn::Symbol)
    list = df[!, nameColumn]
    disc = []
    j = 0
    for i in 2:length(list)
        if list[i] == list[i-1]+1
            j += 1
            continue
        else
            disci = [i-j-1, i-1]
            disc = vcat(disc, disci)
            j = 0
        end
    end
    disci = [length(list)-j, length(list)]
    disc = vcat(disc, disci)
    return(disc)
end

# This function separates the states 1 and 2 into lists of unique st1 and st2 with a unique ID 
# for each list of unique st1 or st2.

function separateStates(type)
    raw = CSV.read("acf-julia/$(type)_20230210_states.csv", DataFrame)
    sel = select(raw, :Column1, :ID, :x, :y, :t, :v, :State)
    sep1 = groupby(sel, :ID)
    data = DataFrame()
    for j in eachindex(sep1)
        X = sep1[j][!,:x]
        Y = sep1[j][!,:y]
        Θ = []
        for t in 2:length(X)-1
            Ri = [X[t]-X[t-1], Y[t]-Y[t-1], 0]
            Rf = [X[t+1]-X[t], Y[t+1]-Y[t], 0]
            RiRf = dot(Ri, Rf)
            absRi = norm(Ri)
            absRf = norm(Rf)
            Θt = acos((RiRf)/(absRi*absRf))
            RixRf = Ri[1]*Rf[2] - Ri[2]*Rf[1]
            if RixRf < 0
                Θt = -Θt
            else
                nothing
            end
            Θ = vcat(Θ, Θt)
        end
        vΘ = sep1[j][2:end-1, :v] .* cos.(Θ)
        dataj = DataFrame(Column1 = sep1[j][2:end-1, :Column1],
                            ID = sep1[j][2:end-1, :ID],
                            x = sep1[j][2:end-1, :x],
                            y = sep1[j][2:end-1, :y],
                            t = sep1[j][2:end-1, :t],
                            v = sep1[j][2:end-1, :v],
                            ang = Θ,
                            vang = vΘ,
                            State = sep1[j][2:end-1, :State]
                            )
        data = vcat(data, dataj)
    end
    CSV.write("acf-julia/$(type)_vang.csv", data)
    sep = groupby(data, :State)
    c = :Column1
    for i in 1:length(sep)
        disc = separateContinuous(sep[i], c)
        for j in 1:2:(length(disc)-1)
            sep[i][disc[j]:disc[j+1],c] = repeat([sep[i][disc[j], c]], disc[j+1]-disc[j]+1)
        end
        final = select(sep[i], c, :v, :vang)
        CSV.write("acf-julia/$(type)$(i)_ang.csv", final)
    end
end

# This runs this code:

types = ["ninv", "inv"]
for type in types
    separateStates(type)
end
