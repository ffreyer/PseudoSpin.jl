# function forward_T2(p1::Float64, m1::Float64, T1::Float64, b1::Float64, m2::Float64, b2::Float64)
#     p = -T1 - m1 * T1 / m2 + T1 / m2 * log(p1) + (b2 - b1) / m2
#     q = m1 * T1^2 / m2 - (b2 - b1) / m2 * T1
#     # print("FORWARD \t p = $p; \t q = $q\n")
#     root = sqrt(0.25p^2-q)
#     -0.5p + root, -0.5p - root
# end
#
# function backward_T2(p2::Float64, m2::Float64, T2::Float64, b2::Float64, m1::Float64, b1::Float64)
#     p = -T2 - m2 * T2 / m1 + T2 / m1 * log(p2) - (b2 - b1) / m1
#     q = m2 * T2^2 / m1 + (b2 - b1) / m1 * T2
#     # print("BACKWARD \t p = $p; \t q = $q\n")
#     root = sqrt(0.25p^2-q)
#     -0.5p + root, -0.5p - root
# end

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
        k = 0       # count against batch_size
        M = div(length(TH_method), adaptive_sample_size)
        switch = 0  # switch between forward and backward propagation
        E_tot = totalEnergy(sgraph, spins, Js, h, g)
        beta = 1.0 / T

        MPI.Barrier(MPI.COMM_WORLD)
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        comm_size = MPI.Comm_size(MPI.COMM_WORLD)

        initial_temperatures = MPI.Allgather(T, MPI.COMM_WORLD)
        betas = 1. ./ initial_temperatures
        T_min, T_max = extrema(initial_temperatures)
        # dT = (T_max - T_min)
        # dT = 0.1minimum(temperatures[2:end] - temperatures[1:end-1])
        cum_prob = 0.0
        N_prob = 0
        cum_E = 0.0

        for i in 1:length(TH_method)
            E_tot = sweep(sgraph, spins, E_tot, Js, beta, h, g)
            cum_E += E_tot

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
                if div(i, adaptive_sample_size) == 1
                    cum_prob = 0.0 #prob_sum
                    N_prob = 0
                    cum_E = 0.
                    print("[$rank] Skip first\n")
                    continue
                end
                # print("\n")
                MPI.Barrier(MPI.COMM_WORLD)
                p = N_prob == 0 ? 0.0 : cum_prob / N_prob
                # probabilities = MPI.Allgather(p, MPI.COMM_WORLD)
                probs = MPI.Allgather(p, MPI.COMM_WORLD)
                # rank == 0 && print("proabilities = $probabilities\n")

                # energies = MPI.Allgather(cum_E, MPI.COMM_WORLD)
                # energies = MPI.Allgather(cum_E / adaptive_sample_size, MPI.COMM_WORLD)
                temperatures = MPI.Allgather(1.0/beta, MPI.COMM_WORLD)
                if (rank != 0) && (rank != comm_size-1)
                    # _probs = probs[1:end-1] + sum(probs[1:end-1])
                    # prob_sum = sum(probs)
                    # T_steps = (T_max - T_min) * cumsum(probs) / prob_sum
                    # Tsteps[1] describes how far Ts[2] should be from Ts[1]
                    # therefore skip rank == 0
                    # print("[$rank] calculated Tsteps = $T_steps\n")
                    # new_T = 0.5(1.0 / beta + T_min + T_steps[rank])

                    # new_T = T_min + dT * (probabilities[rank] - probabilities[rank+1]) / prob_sum

                    # es = energies / sum(energies)
                    # # distance has to scale ~ 1/dE?
                    # T_steps = (T_max - T_min) ./ (cumsum(es[2:end] - es[1:end-1]) + 0.000001)
                    # new_T = T_min + T_steps[rank]

                    # flow = proabilities[rank] - probabilities[rank+1]
                    # Fs = probs[1:end-2] ./ (probs[1:end-2] .+ probs[2:end-1])
                    # push!(Fs, 0.5)
                    # Fs2 = (1-w)Fs + w(1 - k/M)
                    # T_step = (T_max - T_min) .* sum(Fs2[1:rank]) / sum(Fs2)
                    # new_T = T_min + T_step

                    # new_T = 1.0 / beta - (log(probs[rank]) - log(probs[rank+1])) * dT
                    # new_T = 1.0 / beta + (probs[rank] - probs[rank+1]) * dT

                    # equalize p from linear energy fits?
                    # p = mean(exp.(
                    #     (energies[2:end] - energies[1:end-1]) ./
                    #     (1.0 ./ temperatures[2:end] - 1.0 ./ temperatures[1:end-1])
                    # ))
                    # ms = (energies[2:end] - energies[1:end-1]) ./
                    # (temperatures[2:end] - temperatures[1:end-1])
                    # bs = energies[2:end] - ms .* temperatures[2:end]
                    # forward_Ts = deepcopy(temperatures)
                    # for j in 2:comm_size-1
                    #     forward_Ts[j] = forward_T2(
                    #         p, ms[j-1], forward_Ts[j-1], bs[j-1], ms[j], bs[j]
                    #     )[1]
                    # end
                    # # backward_Ts = deepcopy(temperatures)
                    # # for j in comm_size-1:-1:2
                    # #     backward_Ts[j] = backward_T2(
                    # #         p, ms[j], backward_Ts[j], bs[j], ms[j-1], bs[j-1]
                    # #     )[1]
                    # # end
                    # new_T = 0.5(forward_Ts[rank+1] + temperatures[rank+1])
                    # # new_T = forward_Ts[rank+1]

                    # # flatten energies?
                    # E_min, E_max = extrema(energies)
                    # # T = mE + b
                    # ms = (temperatures[2:end] - temperatures[1:end-1]) ./
                    # (energies[2:end] - energies[1:end-1])
                    # bs = temperatures[2:end] - ms .* energies[2:end]
                    #
                    # wanted_E = linspace(E_min, E_max, comm_size)[rank+1]
                    # if cum_E <= wanted_E
                    #     new_T = ms[rank+1] * wanted_E + bs[rank+1]
                    # else
                    #     new_T = ms[rank] * wanted_E + bs[rank]
                    # end

                    prob_sum = sum(probs)
                    norm_prob = probs / prob_sum
                    mean_prob = mean(norm_prob)
                    down = (norm_prob[rank] - mean_prob) *
                        (temperatures[rank+1] - temperatures[rank])
                    up = (norm_prob[rank+1] - mean_prob) *
                        (temperatures[rank+1] - temperatures[rank+2])
                    # down = 0.01(log(norm_prob[rank]) - log(mean_prob)) *
                    #     (temperatures[rank+1] - temperatures[rank])
                    # up = 0.01(log(norm_prob[rank+1]) - log(mean_prob)) *
                    #     (temperatures[rank+1] - temperatures[rank+2])
                    new_T = max(
                        temperatures[rank+1] + (1.0 - k/M) * (down + up),
                        temperatures[rank] + 0.01
                    )
                    # new_T = temperatures[rank+1] + (down + up) * k
                    # unfix: norm_prob[rank], temperatures[rank+1] (from rank+1, rank+2)


                    # T_step = (T_max - T_min) * cumsum(probs) / sum(probs)
                    # new_T = T_min + T_step[rank]
                    beta = 1.0 / new_T
                    sleep(0.01*rank)
                    print("[$rank] T = $(round(temperatures[rank+1], 3))   ->  $(round(new_T, 3)) \t\t down = $down \t up = $up \t\t proabilities = $(probs)\n")
                    # print("[$rank] T = $T   ->  $(round(new_T, 3)) \t\t proabilities = $(probs)\n")
                    # print("[$rank] T = $T   ->  $(round(new_T, 3)) \t\t Energies = $es\n")
                    # print("[$rank] T = $T   ->  $(round(new_T, 3)) \t\t Energies = $Fs\n")
                    # print("[$rank] T = $T   ->  $(round(new_T, 3)) \t\t Energies = $energies \t\t probabilities = $probs\n")
                    # print("[$rank] T = $T   ->  $(round(new_T, 3)) \t\t forward: $forward_Ts \t\t backward $backward_Ts\n")
                    # print("[$rank] T = $T   ->  $(round(new_T, 3))  \t\t forward: $forward_Ts \t\t last $temperatures\n")
                    # print("[$rank] T = $T   ->  $(round(new_T, 3))  \t\t energies = $energies \t\t probs = $probs\n")
                    # print("[$rank] T = $T  ->  $new_T\n")
                end
                k += 1
                cum_prob = 0.0 #prob_sum
                N_prob = 0
                # dT *= 0.75
                cum_E = 0.0
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
