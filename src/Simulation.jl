function thermalize!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        thermalizer::AbstractThermalizationMethod,
        h::Point3{Float64},
        g::Float64,
        E_comp::Compressor
    )
    # print("correct\n")
    init_edges!(sgraph, spins)
    E_tot = totalEnergy(sgraph, spins, Js, h, g)

    if is_parallel(thermalizer)
        beta, state = initialize(thermalizer, T, sgraph, spins)
        while !done(thermalizer, state)
            E_tot = sweep(sgraph, spins, E_tot, Js, beta, h, g)
            beta, E_tot, state = next(thermalizer, state, E_tot)
            push!(E_comp, E_tot)
            yield()
        end
    else
        beta, state = initialize(thermalizer, T)
        while !done(thermalizer, state)
            E_tot = sweep(sgraph, spins, E_tot, Js, beta, h, g)
            beta, state = next(thermalizer, state)
            push!(E_comp, E_tot)
            yield()
        end
    end

    # NOTE Safety check, to be removed
    init_edges!(sgraph, spins)
    E_check = totalEnergy(sgraph, spins, Js, h, g)
    if !(E_tot ≈ E_check)
        warn("E_tot inconsistent after thermalization. $E_tot =/= $(E_check)")
        warn("On process #$(MPI.Comm_rank(MPI.COMM_WORLD))")
        MPI.Finalize()
        exit()
    end

    last(thermalizer, state)
end


function thermalize_no_paths!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        thermalizer::AbstractThermalizationMethod,
        h::Point3{Float64},
        g::Float64,
        E_comp::Compressor
    )
    init_edges!(sgraph, spins)
    E_tot = totalEnergy(sgraph, spins, Js, h, g)

    if is_parallel(thermalizer)
        beta, state = initialize(thermalizer, T, sgraph, spins)
        while !done(thermalizer, state)
            E_tot = sweep_no_paths(sgraph, spins, E_tot, Js, beta, h, g)
            beta, E_tot, state = next(thermalizer, state, E_tot)
            push!(E_comp, E_tot)
            yield()
        end
    else
        beta, state = initialize(thermalizer, T)
        while !done(thermalizer, state)
            E_tot = sweep_no_paths(sgraph, spins, E_tot, Js, beta, h, g)
            beta, state = next(thermalizer, state)
            push!(E_comp, E_tot)
            yield()
        end
    end

    # NOTE Safety check, to be removed
    init_edges!(sgraph, spins)
    E_check = totalEnergy(sgraph, spins, Js, h, g)
    if !(E_tot ≈ E_check)
        warn("E_tot inconsistent after thermalization. $E_tot =/= $(E_check)")
        warn("On process #$(MPI.Comm_rank(MPI.COMM_WORLD))")
        MPI.Finalize()
        exit()
    end

    last(thermalizer, state)
end


function thermalize!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        thermalizer::AbstractThermalizationMethod,
        h::Point3{Float64},
        E_comp::Compressor
    )
    init_edges!(sgraph, spins)
    E_tot = totalEnergy(sgraph, spins, Js, h)

    if is_parallel(thermalizer)
        beta, state = initialize(thermalizer, T, sgraph, spins)
        while !done(thermalizer, state)
            E_tot = sweep(sgraph, spins, E_tot, Js, beta, h)
            beta, E_tot, state = next(thermalizer, state, E_tot)
            push!(E_comp, E_tot)
            yield()
        end
    else
        beta, state = initialize(thermalizer, T)
        while !done(thermalizer, state)
            E_tot = sweep(sgraph, spins, E_tot, Js, beta, h)
            beta, state = next(thermalizer, state)
            push!(E_comp, E_tot)
            yield()
        end
    end

    # NOTE Safety check, to be removed
    init_edges!(sgraph, spins)
    E_check = totalEnergy(sgraph, spins, Js, h)
    if !(E_tot ≈ E_check)
        warn("E_tot inconsistent after thermalization. $E_tot =/= $(E_check)")
        warn("On process #$(MPI.Comm_rank(MPI.COMM_WORLD))")
        MPI.Finalize()
        exit()
    end

    last(thermalizer, state)
end


################################################################################
################################################################################
################################################################################



# Single Temperature simulate
"""
    simulate!(sgraph, spins, sys_size, path, filename, T, Js, TH_method, ME_sweeps, h)

Starts a simulation on a given lattice (sgraph) with given spins, including
thermalization and measurement sweeps. The results will be saved to a file
created in path.
"""
function simulate!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        sys_size::Int64,
        path::String,
        filename::String,
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        thermalizer::AbstractThermalizationMethod,
        ME_sweeps::Int64,
        h::Point3{Float64}=Point3(0.),
        g::Float64 = 0.
    )

    # println("Simulate...")
    # println("\t #therm = ", length(TH_method)) # NOTE
    # println("\t T = ", T)
    # println("\t Js = ", Js)
    # println("\t g = ", g)
    # println("\t h = ", h)
    # println("\t #messure = ", ME_sweeps)

    # Thermalization
    if T <= 0.0
        warn("T = 0.0 is not fully implemented, Assuming T = 1e-10")
        T = 1e-10
    end

    # TODO
    # Add a variable to control the compression level here?
    E_comp = Compressor(1000)

    if g == 0.0
        beta = thermalize!(
            sgraph, spins, T, Js, thermalizer, h, E_comp
        )
    elseif (Js[3][1] == Js[3][2] == 0.0) || (Js[4][1] == Js[4][2] == 0.0)
        beta = thermalize_no_paths!(
            sgraph, spins, T, Js, thermalizer, h, g, E_comp
        )
    else
        beta = thermalize!(
            sgraph, spins, T, Js, thermalizer, h, g, E_comp
        )
    end

    # Fool-proof? file creation
    if !isdir(path)
        println(
            path, " does not exist",
            ". New files will be inserted."
        )
        mkdir(path)
    end
    if isfile(path * filename * ".part")
        println("File ", filename, " exists in path ", path)
        file = open(path * filename * string(time()) * ".part", "w")
    else
        file = open(path * filename * ".part", "w")
    end

    write_header!(
        file, 1, length(thermalizer), T_max(thermalizer), ME_sweeps, sys_size,
        Int64(sgraph.N_nodes), sgraph.K_edges, Js, h, g, max(0.0, 1.0 / beta), #T,
        is_parallel(thermalizer),
        batch_size(thermalizer),
        adaptive_sample_size(thermalizer)
    )

    write_Comp!(file, E_comp, "E-th ")

    # Measurement
    if g == 0.0
        measure!(
            sgraph, spins, beta, Js, file, ME_sweeps, h,
            is_parallel(thermalizer), batch_size(thermalizer)
        )
    elseif (Js[3][1] == Js[3][2] == 0.0) || (Js[4][1] == Js[4][2] == 0.0)
        measure_no_paths!(
            sgraph, spins, beta, Js, file, ME_sweeps, h, g,
            is_parallel(thermalizer), batch_size(thermalizer)
        )
    else
        measure!(
            sgraph, spins, beta, Js, file, ME_sweeps, h, g,
            is_parallel(thermalizer), batch_size(thermalizer)
        )
    end

    close(file)

    nothing
end


function simulate!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        sys_size::Int64,
        path::String,
        filename::String,
        Ts::Vector{Float64},    # <- multiple
        Js::Vector{Tuple{Float64, Float64}},
        thermalizer::AbstractThermalizationMethod,
        ME_sweeps::Int64,
        h::Point3{Float64}=Point3(0.),
        g::Float64 = 0.
    )

    if is_parallel(thermalizer)
        try
            MPI.Init()
        catch e
            println("MPI has to be loaded before starting a simulation!")
            throw(e)
        end
        @assert MPI.Comm_size(MPI.COMM_WORLD) == length(Ts) "The number of processes has to match the number of Temperatures!"
        i = MPI.Comm_rank(MPI.COMM_WORLD)+1
        simulate!(
            sgraph, spins, sys_size,
            path, filename * string(i),
            Ts[i], Js,
            thermalizer, ME_sweeps, h, g
        )
        MPI.Finalize()
    else
        for (i, T) in enumerate(Ts)
            simulate!(
                sgraph, spins, sys_size,
                path, filename * string(i),
                T, Js,
                thermalizer, ME_sweeps, h, g
            )
        end
    end

    nothing
end


"""
    simulate!(;
        path::String = "",
        folder::String = "",
        filename::String = "",

        neighbor_search_depth::Int64 = 2,
        do_paths::Bool = true,
        L::Int64 = 6,
        spins::Union{Vector{Point3{Float64}}, Void} = rand_spin(2*L^3),

        Ts::Vector{Float64} = [1.0],
        J1::Float64 = 0.,
        J2::Float64 = 0.,
        K::Float64 = 0.,
        lambda::Float64 = 0.
        Js::Vector{Tuple{Float64, Float64}} = [
            (J1, lambda*J1),
            (J2, lambda*J2),
            (K, 0.0), (0.0, 1.0)
        ],
        h::Point3{Float64} = Point3(0.),
        g::Float64 = 0.

        TH_sweeps::Int64 = 2_000_000,
        N_switch::Int64 = div(TH_sweeps, 2),
        Freeze_temperature::Float64 = 1.5*maximum(Ts),
        thermalizer::AbstractThermalizationMethod = Freezer(
            N = TH_sweeps,
            T_max = Freeze_temperature,
            N_switch = N_switch
        ),
        ME_sweeps::Int64 = 5_000_000
    )
"""
function simulate!(;
        # files
        path::String = "",
        folder::String = "",
        filename::String = "",
        # Simulation graph
        neighbor_search_depth::Int64 = 2,
        do_paths::Bool = true,
        L::Int64 = 6,
        spins::Vector{Point3{Float64}} = rand_spin(2*L^3),
        # Simulation parameter
        T::Float64 = 1.0,
        Ts::Vector{Float64} = [T],
        J1::Float64 = 0.,
        J2::Float64 = 0.,
        K::Float64 = 0.,
        lambda::Float64 = 0.,
        Js::Vector{Tuple{Float64, Float64}} = [
            (J1, lambda*J1),
            (J2, lambda*J2),
            (K, 0.0), (0.0, 1.0)
        ],
        h::Point3{Float64} = Point3(0.),
        g::Float64 = 0.,
        # Thermalization - Temperature generator
        TH_sweeps::Int64 = 2_000_000,
        N_switch::Int64 = div(TH_sweeps, 2),
        Freeze_temperature::Float64 = 1.5*maximum(Ts),
        TGen_method::Type{<:AbstractTemperatureGenerator} = Freezer,
        # Thermalization - Parallel Tempering
        batch_size::Int64 = 10,
        adaptive_sample_size::Int64 = 100batch_size,
        skip::Int64 = div(TH_sweeps, 2),
        thermalizer_method::Type{
            <:AbstractParallelTemperingAlgorithm
        } = NoParallelTempering,
        # Thermalization - build thermalizer
        thermalizer::AbstractThermalizationMethod = thermalizer_method{
            TGen_method
        }(
            N = TH_sweeps,
            N_switch = N_switch,
            T_max = Freeze_temperature,
            batch_size = batch_size,
            adaptive_sample_size = adaptive_sample_size,
            skip = skip
        ),
        # Measurement
        ME_sweeps::Int64 = 5_000_000
    )

    @assert !is_parallel(thermalizer) || (length(Ts) > 1) "Parallel tempering only\
     works with multiple Temperatures!"

    # println("Do parallel tempering? $do_parallel_tempering")

    # Simulation graph
    rgraph = RGraph(diamond("A"), neighbor_search_depth)
    do_paths && generate_paths!(rgraph)
    sim, flat_index, lattice_indices = Basisfill(rgraph, L)

    if folder == ""
        if !endswith(path, "/")
            path = path * "/"
        end
    else
        if !endswith(folder, "/")
            folder = folder * "/"
        end
        if startswith(folder, "/")
            folder = folder[2:end]
        end
    end

    simulate!(
        sim, spins,L,
        path * folder, filename,
        length(Ts) == 1 ? T : Ts,
        Js,
        thermalizer,
        ME_sweeps,
        h, g
    )

    nothing
end
