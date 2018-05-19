


################################################################################


function thermalize!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        TH_method::Freezer,
        h::Point3{Float64},
        g::Float64
    )
    init_edges!(sgraph, spins)
    for beta in cool_to(TH_method, T)
        sweep(sgraph, spins, Js, beta, h, g)
        yield()
    end
    nothing
end


function thermalize_no_paths!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        TH_method::Freezer,
        h::Point3{Float64},
        g::Float64
    )
    init_edges!(sgraph, spins)
    for beta in cool_to(TH_method, T)
        sweep_no_paths(sgraph, spins, Js, beta, h, g)
        yield()
    end
    nothing
end


function thermalize!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        TH_method::Freezer,
        h::Point3{Float64}
    )
    init_edges!(sgraph, spins)
    for beta in cool_to(TH_method, T)
        sweep(sgraph, spins, Js, beta, h)
        yield()
    end
    nothing
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
        TH_method::Freezer,
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


    # Fool-proof? file creation that was actually not fool-proof
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
        file, 1, length(TH_method), TH_method.T_max, ME_sweeps, sys_size,
        Int64(sgraph.N_nodes), sgraph.K_edges, Js, h, g, T #TODO: order: g, T
    )

    # Thermalization
    init_edges!(sgraph, spins)
    # on_g_branch = g != 0.
    if T > 0.0
        if g == 0.0
            thermalize!(sgraph, spins, T, Js, TH_method, h)
        elseif (Js[3][1] == Js[3][2] == 0.0) || (Js[4][1] == Js[4][2] == 0.0)
            thermalize_no_paths!(sgraph, spins, T, Js, TH_method, h, g)
        else
            thermalize!(sgraph, spins, T, Js, TH_method, h, g)
        end

        beta = 1. / T
    else
        warn("Thermalization skipped, measurement extended.")
        ME_sweeps += length(TH_method)
        beta = -1.0
    end

    # Measurement
    if g == 0.0
        measure!(sgraph, spins, beta, Js, file, ME_sweeps, h)
    elseif (Js[3][1] == Js[3][2] == 0.0) || (Js[4][1] == Js[4][2] == 0.0)
        measure_no_paths!(sgraph, spins, beta, Js, file, ME_sweeps, h, g)
    else
        measure!(sgraph, spins, beta, Js, file, ME_sweeps, h, g)
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
        TH_method::Union{ConstantT, Freezer},
        ME_sweeps::Int64,
        h::Point3{Float64}=Point3(0.),
        g::Float64 = 0.
    )


    for (i, T) in enumerate(Ts)
        simulate!(
            sgraph, spins, sys_size,
            path, filename * string(i),
            T, Js,
            TH_method, ME_sweeps, h, g
        )
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
        ME_sweeps::Int64 = 5_000_000,

        do_parallel_tempering::Bool = false,
        batch_size::Int64 = 1000
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
        # Thermalization/Measurement parameters
        TH_sweeps::Int64 = 2_000_000,
        N_switch::Int64 = div(TH_sweeps, 2),
        Freeze_temperature::Float64 = 1.5*maximum(Ts),
        ME_sweeps::Int64 = 5_000_000,
        # Parallel Tempering
        do_parallel_tempering::Bool = false,
        batch_size::Int64 = 1000
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

    if do_parallel_tempering
        simulate_PT!(
            sim,
            spins,
            L,
            path * folder,
            filename,
            Ts,
            # Ts,
            Js,
            Freezer(TH_sweeps, Freeze_temperature, N_switch=N_switch),
            ME_sweeps,
            h,
            g,
            batch_size
        )
    else
        simulate!(
            sim,
            spins,
            L,
            path * folder,
            filename,
            length(Ts) == 1 ? T : Ts,
            # Ts,
            Js,
            Freezer(TH_sweeps, Freeze_temperature, N_switch=N_switch),
            ME_sweeps,
            h,
            g
        )

    nothing
end