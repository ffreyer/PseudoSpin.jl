"""
    write_header!(
        file, N_points, TH_sweeps, TH_Temp, ME_sweeps, sys_size,
        N_nodes, K_edges, Js, h, T
    )

Writes the file header.
"""
function write_header!(
        file::IOStream,
        N_points::Int64,
        TH_sweeps::Int64,
        TH_Temp::Float64,
        ME_sweeps::Int64,
        sys_size::Int64,
        N_nodes::Int64,
        K_edges::Int64,
        Js::Vector{Tuple{Float64, Float64}},
        h::Point3{Float64},
        g::Float64,
        T::Float64
    )

    write(file, "V03")
    write(file, N_points)
    write(file, TH_sweeps)
    write(file, TH_Temp)
    write(file, ME_sweeps)
    write(file, sys_size)
    write(file, N_nodes)
    write(file, K_edges)

    write(file, length(Js))
    write(file, length(Js[1]))
    for J in Js
        for value in J
            write(file, value)
        end
    end
    for _h in h; write(file, _h) end
    write(file, g)
    write(file, T)

    nothing
end


"""
    write_BA!(file, BinningAnalysis, ID)

Adds data from a Binning Analysis to the file. The section will start with *BA*
and the ID.
"""
function write_BA!(file::IOStream, B::BinnerA, ID::String)
    #### tag
    @assert length(ID) == 5
    write(file, "BA")
    write(file, ID)

    #### metadata
    # binning levels
    write(file, length(B.compressors))
    # Number of values used on each level
    for x in B.count
        write(file, x)
    end
    write(file, Int64(B.i-1))


    #### data
    # means
    for x in B.x_sum
        write(file, x)
    end
    write(file, sum(B.output[1:B.i-1]))

    # variances
    for lvl in 0:length(B.compressors)
        write(file, var(B, lvl))
    end

    # taus
    for lvl in 0:length(B.compressors)
        write(file, tau(B, lvl))
    end

    nothing
end


"""
    write_JK!(file, mean, var, ID)

Adds data from a jackknife analysis to the file. The section will start with
*JK* and the ID.
"""
function write_JK!(file::IOStream, mean::Float64, var::Float64, ID::String)

    # tag
    @assert length(ID) == 5
    write(file, "JK")
    write(file, ID)

    # data
    write(file, mean)
    write(file, var)

    nothing
end


"""
    write_HB!(file, HistogramBinner, ID)

Adds data from a Histogram Binner to the file. The section will start with
*HB* and the ID.
"""
function write_HB!(file::IOStream, B::BinnerH, ID::String)
    # tag
    @assert length(ID) == 5
    write(file, "HB")
    write(file, ID)

    # metadata
    write(file, B.bin_width)
    write(file, length(B.data))

    # data
    for (i, N) in B.data
        write(file, (i + 0.5) * B.bin_width)
        write(file, N)
    end

    nothing
end


"""
    write_SC!(file, spins, ID)

Adds a spin configuration to the file. The section will start with *SC* and the
ID.
"""
function write_SC!(file::IOStream, spins::Vector{Point3{Float64}}, ID::String)
    # tag
    @assert length(ID) == 5
    write(file, "SC") # "Spin Configuration"
    write(file, ID)

    # metadata
    write(file, length(spins))

    # data
    for s in spins
        for x in s
            write(file, x)
        end
    end

    nothing
end
