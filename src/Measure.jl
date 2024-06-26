@inline function flip(i::Int64, s::SVector{3, Float64}, Nhalf::Int64)
    if i <= Nhalf
        return s
    else
        return SVector{3, Float64}(-s[1], -s[2], -s[3])
    end
end



"""
    measure!(
        sgraph,
        spins,
        beta,
        parameters,
        file,
        sweep,
        N_sweeps,
        do_pt,
        batch_size
    )

Starts a measurement (over N_sweeps) on an existing lattice (sgraph) with a
given set of spins and the simulation parameters. The results will be
saved to file. (Note: This does not append a file header)
"""
function measure!(
        sgraph::SGraph,
        spins::Vector{SVector{3, Float64}},
        sampler::Union{Function, AbstractLocalUpdate},
        beta::Float64,
        parameters::Parameters,
        file::IOStream,
        sweep::Function,
        N_sweeps::Int64,
        do_pt::Bool,
        batch_size::Int64,
        do_global_updates::Bool,
        global_rate::Int64,
        global_update::AbstractGlobalUpdate,
        Mhist_cutoff::Float64
    )

    invN = 1. / sgraph.N_nodes
    Nhalf = div(sgraph.N_nodes, 2)

    # Initialize binning analysis
    E_BA = BinnerA(200)
    Es = Array{Float64}(undef, N_sweeps)

    Mx_BA = BinnerA(200);   My_BA = BinnerA(200);   Mz_BA = BinnerA(200)

    M2xy_BA = BinnerA(200);             M2z_BA = BinnerA(200)
    M2xys = Array{Float64}(undef, N_sweeps)
    M2zs = Array{Float64}(undef, N_sweeps)
    Mquad_BA = BinnerA(200);            Moct_BA = BinnerA(200)
    Dimer_xy = BinnerA(200);            Dimer_xy_var = BinnerA(200)
    Dimer_z = BinnerA(200);             Dimer_z_var = BinnerA(200)

    srMx = 0.0;     srMy = 0.0;     srMz = 0.0
    srdMx = 0.0;    srdMy = 0.0;    srdMz = 0.0
    srMxabs = 0.0;  srMyabs = 0.0;  srMzabs = 0.0
    srdMxabs = 0.0; srdMyabs = 0.0; srdMzabs = 0.0

    srM2xy = 0.0;   srM2z = 0.0
    srdM2xy = 0.0;  srdM2z = 0.0

    ssh_binner = SSHBinner(10_000)
    # ssh_binner2 = SSHBinner(10_000)
    circ_hist = CircularHistogram(1080)

    E_tot = totalEnergy(sgraph, spins, parameters)
    spin_sum = sum(spins)
    if do_global_updates
        M1 = normalize(spin_sum)
        new_spins = apply(global_update, spins)
        M2 = normalize(sum(new_spins))
    end

    for i in 1:N_sweeps
        E_tot, spin_sum = sweep(sgraph, spins, spin_sum, sampler, E_tot, beta, parameters)

        if do_global_updates
            if i % global_rate == 0
                new_spins = apply(global_update, spins)
                new_E_tot = totalEnergy(sgraph, new_spins, parameters)
                if (new_E_tot < E_tot) || rand() < exp(-beta * (new_E_tot - E_tot))
                    accept(global_update)
                    E_tot = new_E_tot
                    spins .= new_spins
                    spin_sum = sum(new_spins)
                    if typeof(sampler) in (self_balancing_update, self_balancing_update2)
                        sampler = typeof(sampler)(spins)
                        # sampler.M = sum(spins)
                        # sampler.eM = normalize(sampler.M)
                        # sampler.eM_perp = cross(SVector(0., 0., 1.), sampler.eM)
                    end
                end
            end
        end

        if parameters.dual_rot && (i % 1000 == 0)
            @inbounds for j in eachindex(spins)
                n = norm(spins[j])
                n ≈ 1.0 || begin
                    @warn "Normalization actually necessary"
                    spins[j] = spins[j] / n
                end
            end
            spin_sum = sum(spins)
            m = normalize(spin_sum)
            if (!(m ≈ M1)) && (!(m ≈ M2))
                @warn "magnetization changed!"
                printstyled("\tM1 = $M1\n")
                printstyled("\tm  = $m\n", color=:red)
                printstyled("\tM2 = $M2\n")
                M = m
            end
        end

        @inbounds Es[i] = E_tot * invN
        push!(E_BA, E_tot * invN)

        # We're only updating spin_sum when the kappa term is used
        parameters.kappa == 0.0 && (spin_sum = sum(spins))
        S = spin_sum * invN
        # _norm = sum(S.^2)
        push!(circ_hist, normalize(spin_sum))
        # if true #_norm > Mhist_cutoff # norm > 0.32
        #     push!(ssh_binner2, S ./ _norm)
        # end

        @inbounds push!(Mx_BA, S[1])
        @inbounds push!(My_BA, S[2])
        @inbounds push!(Mz_BA, S[3])

        @inbounds M2xys[i] = sqrt(S[1]^2 + S[2]^2)
        @inbounds M2zs[i] = abs(S[3])
        @inbounds push!(M2xy_BA, M2xys[i])
        @inbounds push!(M2z_BA, abs(S[3]))

        @inbounds temp = mapreduce(
            v -> SVector{3, Float64}(0., sqrt(v[1] * v[1] + v[2] * v[2]), abs(v[3])),
            +,
            spins
        ) * invN
        @inbounds push!(Mquad_BA, temp[2])
        @inbounds push!(Moct_BA, temp[3])


        # 1/N ∑_sites ∑_{e ∈ NN} e.xy
        Dxy = map(n -> mapreduce(e -> e.xy, +, n.first), sgraph.nodes)
        Dxy_mean = sum(Dxy) * invN
        Dxy_var = sum(Dxy.^2) * invN - Dxy_mean^2
        Dz = map(n -> mapreduce(e -> e.z, +, n.first), sgraph.nodes)
        Dz_mean = sum(Dz) * invN
        Dz_var = sum(Dz.^2) * invN - Dz_mean^2

        @inbounds push!(Dimer_xy, Dxy_mean)
        @inbounds push!(Dimer_xy_var, Dxy_var)
        @inbounds push!(Dimer_z, Dz_mean)
        @inbounds push!(Dimer_z_var, Dz_var)


        _spins = map(t -> flip(t..., Nhalf), enumerate(spins))
        S = reduce(+, _spins) * invN
        srMx += S[1]
        srMy += S[2]
        srMz += S[3]

        # NOTE
        # map, ... used with a function containing a constant causes that
        # constant to be boxed and lose its type.
        # This is a workaround for this
        temp2 = similar(_spins)
        for j in eachindex(temp2); temp2[j] = _spins[j] - S; end
        vars = mapreduce(s -> s.^2, +, temp2) * invN
        srdMx += vars[1]
        srdMy += vars[2]
        srdMz += vars[3]

        srM2xy += sqrt(S[1]^2 + S[2]^2)
        srdM2xy += sqrt(vars[1]^2 + vars[2]^2)

        S = mapreduce(v -> abs.(v), +, _spins) * invN
        srMxabs += S[1]
        srMyabs += S[2]
        srMzabs += S[3]

        temp2 = similar(_spins)
        for j in eachindex(temp2); temp2[j] = abs.(_spins[j]) - S; end
        vars = mapreduce(s -> (s).^2, +, temp2) * invN
        srdMxabs += vars[1]
        srdMyabs += vars[2]
        srdMzabs += vars[3]

        append!(ssh_binner, spins)

        if do_pt && (i % batch_size == 0)
            E_tot = _parallel_tempering!(sgraph, spins, E_tot, beta)
            __switch__[] = 1 - __switch__[]
        end

        yield()
    end

    init_edges!(sgraph, spins)
    E_check = totalEnergy(sgraph, spins, parameters)
    if !(E_tot ≈ E_check)
        @warn(
            "E_tot diverged by $(round(100(E_check - E_tot) / E_check, 3))%" *
            " (E_check = $E_check, E_tot = $E_tot)."
        )
    end

    # TODO - This still has type instability due to constants in lambda function
    # This formula is wrong for k =/= 1
    # dc/dT?
    f(x, x2, x3) = (
        (x3 - 3*x*x2 + 2*x^3) * beta^4 * sgraph.N_nodes -
        2 * (x2 - x^2) * beta^3
    ) * sgraph.N_nodes

    cv, dcv = jackknife((x, x2) -> (x2 - x^2) * beta^2 * sgraph.N_nodes, Es, Es.^2)
    dcvdT, ddcvdT = jackknife(f, Es, Es.^2, Es.^3)
    # TODO - This still has type instability due to constants in lambda function

    dM2xy, ddM2xy = jackknife((x, x2) -> (x2 - x^2), M2xys, M2xys.^2)
    dM2z, ddM2z = jackknife((x, x2) -> (x2 - x^2), M2zs, M2zs.^2)

    # maybe do this for Ms too?
    Es_HB = BinnerH(0.0001)
    for E in Es; push!(Es_HB, E) end

    # saving
    write_BA!(file, E_BA, "Energ")
    write_BA!(file, Mx_BA, "M1x  ")
    write_BA!(file, My_BA, "M1y  ")
    write_BA!(file, Mz_BA, "M1z  ")
    write_BA!(file, M2xy_BA, "M2xy ")
    write_BA!(file, M2z_BA,  "M2z  ")
    write_BA!(file, Mquad_BA, "Mquad")
    write_BA!(file, Moct_BA, "Moct ")

    write_BA!(file, Dimer_xy, "Dxy  ")
    write_BA!(file, Dimer_xy_var, "DxyV ")
    write_BA!(file, Dimer_z, "Dz   ")
    write_BA!(file, Dimer_z_var, "DzV  ")

    K = 1. / N_sweeps
    write_JK!(file, srMx * K, srdMx * K, "rMx  ")
    write_JK!(file, srMy * K, srdMy * K, "rMy  ")
    write_JK!(file, srMz * K, srdMz * K, "rMz  ")
    write_JK!(file, srMxabs * K, srdMxabs * K, "rMxa ")
    write_JK!(file, srMyabs * K, srdMyabs * K, "rMya ")
    write_JK!(file, srMzabs * K, srdMzabs * K, "rMza ")
    write_JK!(file, srM2xy * K, srdM2xy * K, "rMxy ")

    write_JK!(file, cv, dcv, "cV   ")
    write_JK!(file, dcvdT, ddcvdT, "dcVdT")
    write_JK!(file, dM2xy, ddM2xy, "dM2xy")
    write_JK!(file, dM2z, ddM2z, "dM2z ")

    write_HB!(file, Es_HB, "Energ")
    write_SSHB!(file, ssh_binner, "spins")
    # write_SSHB!(file, ssh_binner2, "Magn ")
    write_SSHB!(file, circ_hist, "Magn ")
    write_SC!(file, spins, "spins")

    nothing
end
