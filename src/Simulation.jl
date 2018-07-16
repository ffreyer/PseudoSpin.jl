function thermalize!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        parameters::Parameters,
        thermalizer::AbstractThermalizationMethod,
        sweep::Function,
        E_comp::Compressor
    )
    # print("correct\n")
    init_edges!(sgraph, spins)
    E_tot = totalEnergy(sgraph, spins, parameters)

    if is_parallel(thermalizer)
        beta, state = initialize(thermalizer, T, sgraph, spins)
        while !done(thermalizer, state)
            E_tot = sweep(sgraph, spins, E_tot, beta, parameters)
            beta, E_tot, state = next(thermalizer, state, E_tot)
            push!(E_comp, E_tot)
            yield()
        end
    else
        beta, state = initialize(thermalizer, T)
        while !done(thermalizer, state)
            E_tot = sweep(sgraph, spins, E_tot, beta, parameters)
            beta, state = next(thermalizer, state)
            push!(E_comp, E_tot)
            yield()
        end
    end

    # NOTE Safety check, can be removed
    init_edges!(sgraph, spins)
    E_check = totalEnergy(sgraph, spins, parameters)
    if !(E_tot ≈ E_check)
        warn("E_tot inconsistent after thermalization. $E_tot =/= $(E_check)")
        warn("On process #$(MPI.Comm_rank(MPI.COMM_WORLD))")
        MPI.Finalize()
        exit()
    end

    last(thermalizer, state)
end


# Single Temperature simulate
# """
#     simulate!(
#         sgraph,
#         spins,
#         sys_size,
#         path,
#         filename,
#         T,
#         parameters,
#         thermalizer,
#         ME_sweeps
#     )
#
# Starts a simulation on a given lattice (sgraph) with given spins for a single
# temperature T. The results will be saved to a file filename created in path.
# """
function simulate!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        sys_size::Int64,
        path::String,
        filename::String,
        T::Float64,
        parameters::Parameters,
        thermalizer::AbstractThermalizationMethod,
        ME_sweeps::Int64,
    )

    # Thermalization
    if T <= 0.0
        warn("T = 0.0 is not fully implemented, Assuming T = 1e-10")
        T = 1e-10
    end

    # TODO
    # Add a variable to control the compression level here?
    E_comp = Compressor(1000)
    sweep = sweep_picker(parameters)
    beta = thermalize!(sgraph, spins, T, parameters, thermalizer, sweep, E_comp)

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
        Int64(sgraph.N_nodes), sgraph.K_edges,
        parameters,
        max(0.0, 1.0 / beta),
        is_parallel(thermalizer),
        batch_size(thermalizer),
        adaptive_sample_size(thermalizer)
    )

    write_Comp!(file, E_comp, "E_th ")

    # Measurement
    measure!(
        sgraph, spins, beta, parameters, file, sweep, ME_sweeps,
        is_parallel(thermalizer), batch_size(thermalizer)
    )

    close(file)

    nothing
end


# Multi temperature simulate!
# """
#     simulate!(
#         sgraph,
#         spins,
#         sys_size,
#         path,
#         filename,
#         Ts,
#         parameters,
#         thermalizer,
#         ME_sweeps
#     )
#
# Starts a simulation on a given lattice (sgraph) with given spins for a set of
# temperatures Ts. The results will be saved to a file filename1...filenameN
# created in path.
# """
function simulate!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        sys_size::Int64,
        path::String,
        filename::String,
        Ts::Vector{Float64},    # <- multiple
        parameters::Parameters,
        thermalizer::AbstractThermalizationMethod,
        ME_sweeps::Int64
    )

    if is_parallel(thermalizer)
        MPI.Init()
        @assert MPI.Comm_size(MPI.COMM_WORLD) == length(Ts) "The number of processes has to match the number of Temperatures!"
        i = MPI.Comm_rank(MPI.COMM_WORLD)+1
        simulate!(
            sgraph, spins, sys_size,
            path, filename * string(i),
            Ts[i], parameters,
            thermalizer, ME_sweeps
        )
        MPI.Finalize()
    else
        for (i, T) in enumerate(Ts)
            simulate!(
                sgraph, spins, sys_size,
                path, filename * string(i),
                T, parameters,
                thermalizer, ME_sweeps, h, g
            )
        end
    end

    nothing
end


"""
    simulate!(;
        # location/name of data files (results)
        path::String = "",
        folder::String = "",
        filename::String = "T",

        # Simulation graph setup
        neighbor_search_depth::Int64 = 2,
        do_paths::Bool = true,
        L::Int64 = 6,
        spins::Vector{Point3{Float64}} = rand_spin(2*L^3),

        # Temperatures
        T::Float64 = 1.0,
        Ts::Vector{Float64} = [T],

        # Parameters
        J1::Float64 = 0.,
        J2::Float64 = 0.,
        lambda::Float64 = 0.,
        J1s::NTuple{2, Float64} = (J1, lambda*J1),
        J2s::NTuple{2, Float64} = (J2, lambda*J2),
        K::Float64 = 0.,
        h::Point3{Float64} = Point3(0.),
        g::Float64 = 0.,
        parameters::Parameters = Parameters(
            J1s = J1s,
            J2s = J2s,
            K = K,
            g = g,
            h = h
        ),

        # Thermalization - Temperature generator
        # Number of thermalization sweeps
        TH_sweeps::Int64 = 2_000_000,
        # Number of steps after which Freezer becomes constant
        N_switch::Int64 = div(TH_sweeps, 2),
        # Starting temperature of Freezer
        Freeze_temperature::Float64 = 1.5*maximum(Ts),
        # can be: Freezer (simulated annealing) or ConstantT
        TGen_method::Type{<:AbstractTemperatureGenerator} = Freezer,

        # Thermalization - Parallel Tempering
        # Number of sweeps between parallel_tempering exchanges
        batch_size::Int64 = 10,
        # Number of sweeps between adaptive steps
        adaptive_sample_size::Int64 = 100batch_size,
        # Number of sweeps skipped in adaptive algorithm (to thermalize system)
        skip::Int64 = div(TH_sweeps, 2),
        # can be: NoParallelTempering, ParallelTempering, ProbabilityEqualizer
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

        # Number of measurement sweeps
        ME_sweeps::Int64 = 5_000_000
    )

Starts a simulation with the given keyword arguments.
"""
function simulate!(;
        # files
        path::String = "",
        folder::String = "",
        filename::String = "T",
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
        lambda::Float64 = 0.,
        J1s::NTuple{2, Float64} = (J1, lambda*J1),
        J2s::NTuple{2, Float64} = (J2, lambda*J2),
        K::Float64 = 0.,
        h::Point3{Float64} = Point3(0.),
        g::Float64 = 0.,
        zeta::Float64 = 0.,
        parameters::Parameters = Parameters(
            J1s = J1s,
            J2s = J2s,
            K = K,
            g = g,
            h = h,
            zeta = zeta
        ),
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
        parameters,
        thermalizer,
        ME_sweeps
    )

    nothing
end
