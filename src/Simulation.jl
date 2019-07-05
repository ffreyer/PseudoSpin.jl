function thermalize!(
        sgraph::SGraph,
        spins::Vector{SVector{3, Float64}},
        sampler::Function,
        T::Float64,
        parameters::Parameters,
        thermalizer::AbstractThermalizationMethod,
        sweep::Function,
        E_comp::Compressor
    )
    init_edges!(sgraph, spins)
    E_tot = totalEnergy(sgraph, spins, parameters)

    if is_parallel(thermalizer)
        beta, state = initialize(thermalizer, T, sgraph, spins)
        while !done(thermalizer, state)
            E_tot = sweep(sgraph, spins, sampler, E_tot, beta, parameters)
            beta, E_tot, state = next(thermalizer, state, E_tot)
            push!(E_comp, E_tot)
            yield()
        end
    else
        beta, state = initialize(thermalizer, T)
        while !done(thermalizer, state)
            E_tot = sweep(sgraph, spins, sampler, E_tot, beta, parameters)
            beta, state = next(thermalizer, state)
            push!(E_comp, E_tot)
            yield()
        end
    end

    # NOTE Safety check
    init_edges!(sgraph, spins)
    E_check = totalEnergy(sgraph, spins, parameters)
    if !(E_tot ≈ E_check)
        error(
            "E_tot inconsistent after thermalization. $E_tot =/= $(E_check)" *
            " on process #$(myid()) ($(corr_id[]))"
        )
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
        spins::Vector{SVector{3, Float64}},
        sys_size::Int64,
        sampler::Function,
        path::String,
        filename::String,
        T::Float64,
        parameters::Parameters,
        thermalizer::AbstractThermalizationMethod,
        ME_sweeps::Int64,
    )

    # Thermalization
    if T <= 0.0
        @warn("T = 0.0 is not fully implemented, Assuming T = 1e-10")
        T = 1e-10
    end

    # TODO
    # Add a variable to control the compression level here?
    E_comp = Compressor(1000)
    sweep = sweep_picker(parameters)
    beta = thermalize!(
        sgraph, spins, sampler,
        T, parameters, thermalizer, sweep, E_comp
    )

    # Fool-proof? file creation
    if !isdir(path)
        if corr_id[] == 1
            println("'$path' does not exist. Path will be created.")
            mkdir(path)
        else
            println(
                "'path' does not exist. Waiting for main process to create " *
                "path."
            )
        end
    end
    if isfile(path * filename * ".part")
        stats = open(path * filename * ".part") do f; stat(f) end
        if stats.size == 0
            @warn(
                "File '$filename' already exists in path '$path', but it is " *
                "empty. Overwriting '$filename'."
            )
            file = open(path * filename * ".part", "w")
        else
            println(
                "File '$filename' already exists in path '$path'. Since it " *
                "is not empty, a new file will be created."
            )
            file = open(path * filename * string(time()) * ".part", "w")
        end
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
        sgraph, spins, sampler,
        beta, parameters, file, sweep, ME_sweeps,
        is_parallel(thermalizer), batch_size(thermalizer)
    )

    close(file)

    nothing
end


const start_id = Ref{Int64}(0)
const stop_id = Ref{Int64}(0)
const corr_id = Ref{Int64}(0)
function simulate!(
        sgraph::SGraph,
        spins::Vector{SVector{3, Float64}},
        sys_size::Int64,
        sampler::Function,
        path::String,
        filename::String,
        Ts::Vector{Float64},    # <- multiple
        parameters::Parameters,
        thermalizer::AbstractThermalizationMethod,
        ME_sweeps::Int64
    )

    if is_parallel(thermalizer)
        np = nprocs()
        if np == length(Ts)
            start_id[] = 1
            stop_id[] = length(Ts)
            corr_id[] = myid()
        elseif np == length(Ts) + 1
            start_id[] = 2
            stop_id[] = length(Ts) + 1
            corr_id[] = myid() - 1
        else
            error(
                "The number of processes ($(nprocs())) has to match the " *
                "number of Temperatures ($(length(Ts)))!"
            )
        end
        simulate!(
            sgraph, spins, sys_size, sampler,
            path, filename * string(corr_id[]),
            Ts[corr_id[]], parameters,
            thermalizer, ME_sweeps
        )
    else
        for (i, T) in enumerate(Ts)
            simulate!(
                sgraph, spins, sys_size, sampler,
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
        spins::Vector{SVector{3, Float64}} = rand_spin(2*L^3),

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
        h::SVector{3, Float64} = SVector(0.0, 0.0, 0.0),
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
        neighbor_search_depth::Int64 = 3,
        do_paths::Bool = true,
        L::Int64 = 6,
        spins::Vector{SVector{3, Float64}} = rand_spin(2*L^3),
        sampler::Function = rand_spin,
        # Simulation parameter
        T::Float64 = 1.0,
        Ts::Vector{Float64} = [T],
        J1::Float64 = 0.,
        J2::Float64 = 0.,
        J3::Float64 = 0.,
        lambda::Float64 = 0.,
        J1s::NTuple{2, Float64} = (J1, lambda*J1),
        J2s::NTuple{2, Float64} = (J2, lambda*J2),
        J3s::NTuple{2, Float64} = (J3, lambda*J3),
        K::Float64 = 0.,
        h::SVector{3, Float64} = SVector(0.0, 0.0, 0.0),
        g::Float64 = 0.,
        zeta::Float64 = 0.,
        parameters::Parameters = Parameters(
            J1s = J1s,
            J2s = J2s,
            J3s = J3s,
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

    @assert(
        !is_parallel(thermalizer) || (length(Ts) > 1),
        "Parallel tempering only works with multiple Temperatures!"
    )

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
        sim, spins, L, sampler,
        path * folder, filename,
        length(Ts) == 1 ? Ts[1] : Ts,
        parameters,
        thermalizer,
        ME_sweeps
    )

    nothing
end
