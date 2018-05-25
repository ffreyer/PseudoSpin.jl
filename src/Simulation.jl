function thermalize!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        TH_method::AbstractTGen,
        h::Point3{Float64},
        g::Float64,
        do_pt::Bool,
        do_adaptive::Bool,
        batch_size::Int64,
        adaptive_sample_size::Int64 = 100batch_size
    )
    print("correct\n")
    init_edges!(sgraph, spins)

    if do_adaptive
        i = 0       # count against batch_size
        switch = 0  # switch between forward and backward propagation
        E_tot = totalEnergy(sgraph, spins, Js, h, g)
        beta = 1.0 / T

        MPI.Barrier(MPI.COMM_WORLD)
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        comm_size = MPI.Comm_size(MPI.COMM_WORLD)

        temperatures = MPI.Allgather(T, MPI.COMM_WORLD)
        T_min, T_max = extrema(temperatures)
        cum_prob = 1.0
        N_prob = 0

        for i in 1:length(TH_method)
            E_tot = sweep(sgraph, spins, E_tot, Js, beta, h, g)

            if i % batch_size == 0
                E_tot, p = parallel_tempering_adaptive!(
                    sgraph, spins, E_tot, beta, switch
                )
                if p >= 0.0
                    # print("[$rank] p = $p\n")
                    cum_prob += p
                    N_prob += 1
                end
                switch = 1 - switch
            end

            if i % adaptive_sample_size == 0
                # print("\n")
                MPI.Barrier(MPI.COMM_WORLD)
                p = N_prob == 0 ? 0.0 : cum_prob / N_prob
                probabilities = MPI.Allgather(p, MPI.COMM_WORLD)
                # rank == 0 && print("proabilities = $probabilities\n")
                if (rank != 0) && (rank != comm_size-1)
                    T_steps = (T_max - T_min) * cumsum(probabilities) / sum(probabilities)
                    # Tsteps[1] describes how far Ts[2] should be from Ts[1]
                    # therefore skip rank == 0
                    # print("[$rank] calculated Tsteps = $T_steps\n")
                    new_T = (T_min + T_steps[rank])
                    beta = 1.0 / new_T
                    print("[$rank] T = $T   ->  $(round(new_T, 3)) \t\t proabilities = $(probabilities[1:end-1])\n")
                    # print("[$rank] T = $T  ->  $new_T\n")
                end
                cum_prob = 1.0
                N_prob = 0
            end
            yield()
        end

        E_check = totalEnergy(sgraph, spins, Js, h, g)
        if !(E_tot ≈ E_check)
            warn(
                "E_tot inconsistent on $(MPI.Comm_rank(MPI.COMM_WORLD)) after $i\
                $E_tot =/= $(E_check)\
                in full thermalize."
            )
            MPI.Finalize()
            exit()
        end

    elseif do_pt
        i = 0       # count against batch_size
        switch = 0  # switch between forward and backward propagation
        E_tot = totalEnergy(sgraph, spins, Js, h, g)

        # blocked_time = 0.0
        # total_time = -time()

        for beta in cool_to(TH_method, T)
            E_tot = sweep(sgraph, spins, E_tot, Js, beta, h, g)
            i += 1

            # parallel tempering step
            if i % batch_size == 0
                E_tot = parallel_tempering!(sgraph, spins, E_tot, beta, switch)
                # E_tot, bt = parallel_tempering_time!(spins, E_tot, beta, switch)
                # blocked_time += bt
                # init_edges!(sgraph, spins)
                switch = 1 - switch
            end
            yield()
        end

        E_check = totalEnergy(sgraph, spins, Js, h, g)
        if !(E_tot ≈ E_check)
            warn(
                "E_tot inconsistent on $(MPI.Comm_rank(MPI.COMM_WORLD)) after $i\
                $E_tot =/= $(E_check)\
                in full thermalize."
            )
            # MPI.Barrier()
            MPI.Finalize()
            exit()
        end
        # total_time += time()
        # print("[$(MPI.Comm_rank(MPI.COMM_WORLD))] was blocked for $blocked_time / $total_time = $(round(blocked_time/total_time*100, 1))%. \n")
    else
        tic()
        for beta in cool_to(TH_method, T)
            sweep(sgraph, spins, Js, beta, h, g)
            yield()
        end
        print("[T = $T] $(toq())\n")
    end

    1. / T
end


function thermalize_no_paths!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        TH_method::AbstractTGen,
        h::Point3{Float64},
        g::Float64,
        do_pt::Bool,
        do_adaptive::Bool,
        batch_size::Int64
    )
    init_edges!(sgraph, spins)
    if do_pt
        i = 0       # count against batch_size
        switch = 0  # switch between forward and backward propagation
        E_tot = totalEnergy(sgraph, spins, Js, h, g)

        for beta in cool_to(TH_method, T)
            E_tot = sweep_no_paths(sgraph, spins, E_tot, Js, beta, h, g)
            i += 1

            # parallel tempering step
            if i % batch_size == 0
                E_tot = parallel_tempering!(sgraph, spins, E_tot, beta, switch)
                # init_edges!(sgraph, spins)
                switch = 1 - switch
            end
            yield()
        end

        E_check = totalEnergy(sgraph, spins, Js, h, g)
        if !(E_tot ≈ E_check)
            warn(
                "E_tot inconsistent on $(MPI.Comm_rank(MPI.COMM_WORLD)) after $i\
                $E_tot =/= $(E_check)\
                in full thermalize."
            )
            # MPI.Barrier()
            MPI.Finalize()
            exit()
        end
    else
        for beta in cool_to(TH_method, T)
            sweep_no_paths(sgraph, spins, Js, beta, h, g)
            yield()
        end
    end

    1. / T
end


function thermalize!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        T::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        TH_method::AbstractTGen,
        h::Point3{Float64},
        do_pt::Bool,
        do_adaptive::Bool,
        batch_size::Int64
    )
    init_edges!(sgraph, spins)
    if do_pt
        i = 0       # count against batch_size
        switch = 0  # switch between forward and backward propagation
        E_tot = totalEnergy(sgraph, spins, Js, h)

        for beta in cool_to(TH_method, T)
            E_tot = sweep(sgraph, spins, E_tot, Js, beta, h)
            i += 1

            # parallel tempering step
            if i % batch_size == 0
                E_tot = parallel_tempering!(sgraph, spins, E_tot, beta, switch)
                # E_tot, _ = parallel_tempering_time!(spins, E_tot, beta, switch)
                # init_edges!(sgraph, spins)
                switch = 1 - switch
            end
            yield()
        end

        init_edges!(sgraph, spins)
        E_check = totalEnergy(sgraph, spins, Js, h)
        if !(E_tot ≈ E_check)
            warn(
                "E_tot inconsistent on $(MPI.Comm_rank(MPI.COMM_WORLD)) after $i\
                $E_tot =/= $(E_check)\
                in full thermalize."
            )
            # MPI.Barrier()
            MPI.Finalize()
            exit()
        end
    else
        for beta in cool_to(TH_method, T)
            sweep(sgraph, spins, Js, beta, h)
            yield()
        end
    end

    1. / T
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
        TH_method::AbstractTGen,
        ME_sweeps::Int64,
        h::Point3{Float64}=Point3(0.),
        g::Float64 = 0.,
        do_pt::Bool = false,
        do_adaptive::Bool = false,
        batch_size::Int64 = 1000
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
        file, 1, length(TH_method), last(TH_method), ME_sweeps, sys_size,
        Int64(sgraph.N_nodes), sgraph.K_edges, Js, h, g, T,
        do_pt, batch_size
    )

    # Thermalization
    # init_edges!(sgraph, spins)
    # on_g_branch = g != 0.
    if T > 0.0
        if g == 0.0
            beta = thermalize!(
                sgraph, spins, T, Js, TH_method, h,
                do_pt, do_adaptive, batch_size
            )
        elseif (Js[3][1] == Js[3][2] == 0.0) || (Js[4][1] == Js[4][2] == 0.0)
            beta = thermalize_no_paths!(
                sgraph, spins, T, Js, TH_method, h, g,
                do_pt, do_adaptive, batch_size
            )
        else
            beta = thermalize!(
                sgraph, spins, T, Js, TH_method, h, g,
                do_pt, do_adaptive, batch_size
            )
        end
    else
        warn("Thermalization skipped, measurement extended.")
        ME_sweeps += length(TH_method)
        beta = -1.0
    end

    # Measurement
    if g == 0.0
        measure!(sgraph, spins, beta, Js, file, ME_sweeps, h, do_pt, batch_size)
    elseif (Js[3][1] == Js[3][2] == 0.0) || (Js[4][1] == Js[4][2] == 0.0)
        measure_no_paths!(
            sgraph, spins, beta, Js, file, ME_sweeps, h, g, do_pt, batch_size
        )
    else
        measure!(
            sgraph, spins, beta, Js, file, ME_sweeps, h, g, do_pt, batch_size
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
        TH_method::AbstractTGen,
        ME_sweeps::Int64,
        h::Point3{Float64}=Point3(0.),
        g::Float64 = 0.,
        do_parallel_tempering::Bool = false,
        do_adaptive::Bool = false,
        batch_size::Int64 = 1000
    )

    if do_parallel_tempering || do_adaptive
        MPI.Init()
        @assert MPI.Comm_size(MPI.COMM_WORLD) == length(Ts) "The number of processes has to match the number of Temperatures!"
        i = MPI.Comm_rank(MPI.COMM_WORLD)+1
        simulate!(
            sgraph, spins, sys_size,
            path, filename * string(i),
            Ts[i], Js,
            TH_method, ME_sweeps, h, g,
            do_parallel_tempering || do_adaptive, do_adaptive, batch_size
        )
        MPI.Finalize()
    else
        for (i, T) in enumerate(Ts)
            simulate!(
                sgraph, spins, sys_size,
                path, filename * string(i),
                T, Js,
                TH_method, ME_sweeps, h, g,
                do_parallel_tempering, batch_size
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
        ME_sweeps::Int64 = 5_000_000,
        TH_Method::AbstractTGen = if do_parallel_tempering
            ConstantT(TH_sweeps)
        else
            Freezer(TH_sweeps, Freeze_temperature, N_switch=N_switch)
        end,

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
        # Parallel Tempering
        do_parallel_tempering::Bool = false,
        do_adaptive::Bool = false,
        batch_size::Int64 = 1000,
        # Thermalization/Measurement parameters
        TH_sweeps::Int64 = 2_000_000,
        N_switch::Int64 = div(TH_sweeps, 2),
        Freeze_temperature::Float64 = 1.5*maximum(Ts),
        TH_method::DataType = if do_parallel_tempering || do_adaptive
            ConstantT
        else
            Freezer
        end,
        ME_sweeps::Int64 = 5_000_000
    )

    @assert !do_parallel_tempering || (length(Ts) > 1) "Parallel tempering only\
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
        TH_method(TH_sweeps, Freeze_temperature, N_switch=N_switch),
        ME_sweeps,
        h, g,
        do_parallel_tempering,
        do_adaptive,
        batch_size
    )

    nothing
end
