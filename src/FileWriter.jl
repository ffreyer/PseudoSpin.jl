"""
    write_header!(
        file, N_points, TH_sweeps, TH_Temp, ME_sweeps, sys_size,
        N_nodes, K_edges, Js, h, T
    )

Writes the file header.
"""
function write_header!(
        file::IOStream;
        N_points::Int64,
        TH_sweeps::Int64,
        TH_Temp::Float64,
        ME_sweeps::Int64,
        sys_size::Int64,
        N_nodes::Int64,
        K_edges::Int64,
        parameters::Parameters,
        T::Float64,
        do_parallel_tempering::Bool,
        batch_size::Int64,
        adaptive_sample_size::Int64,
        sampler::Union{Function, AbstractLocalUpdate},
        sweep::Function,
        do_global_updates::Bool,
        global_rate::Int64,
        global_update::AbstractGlobalUpdate,
        Mhist_cutoff::Float64
    )

    write(file, "V08")
    write(file, N_points)
    write(file, TH_sweeps)
    write(file, TH_Temp)
    write(file, ME_sweeps)
    write(file, sys_size)
    write(file, N_nodes)
    write(file, K_edges)

    write(file, 5);                 write(file, 2)
    write(file, parameters.J1[1]);  write(file, parameters.J1[2])
    write(file, parameters.J2[1]);  write(file, parameters.J2[2])
    write(file, parameters.J3[1]);  write(file, parameters.J3[2])
    write(file, parameters.K);      write(file, 0.0)
    write(file, 0.0);               write(file, 1.0)
    for _h in parameters.h; write(file, _h) end
    write(file, parameters.g)
    write(file, parameters.zeta)
    write(file, T)

    write(file, do_parallel_tempering)
    write(file, batch_size)
    write(file, adaptive_sample_size)

    write(file, "&!" * string(sampler) * "&!")
    write(file, "&!" * string(sweep) * "&!")

    write(file, do_global_updates)
    write(file, global_rate)
    write(file, "&!" * string(typeof(global_update)) * "&!")

    write(file, Mhist_cutoff)

    nothing
end


"""
    write_Comp!(file, Compressor, ID)

Adds data from a Compressor to the file. The section will start with *CP*
and the ID.
"""
function write_Comp!(file::IOStream, c::Compressor, ID::String)
    #### tag
    @assert length(ID) == 5
    write(file, "CP")
    write(file, ID)

    out = output(c)
    write(file, length(out))
    for x in out
        write(file, x)
    end

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
    write_SSHB!(file, SSHBinner, ID)

Adds data from a SphereSurfaceHistogram Binner to the file. The section will
start with *SH* and the ID.
"""
function write_SSHB!(file::IOStream, B::SSHBinner, ID::String)
    # tag
    @assert length(ID) == 5
    write(file, "SH")
    write(file, ID)

    # metadata
    write(file, length(B.bins))

    # data
    for bin in B.bins
        write(file, bin)
    end

    nothing
end


"""
    write_SC!(file, spins, ID)

Adds a spin configuration to the file. The section will start with *SC* and the
ID.
"""
function write_SC!(file::IOStream, spins::Vector{SVector{3, Float64}}, ID::String)
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
