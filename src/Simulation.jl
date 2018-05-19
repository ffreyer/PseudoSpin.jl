# Just measure! and simulate!

"""
    measure!(sgraph, spins, beta, Js, file, N_sweeps, h)

Starts a measurement (over N_sweeps) on an existing lattice (sgraph) with a
given set of spins and the simulation parameters Js and h. The results will be
saved to file. (Note: This does not append a file header)
"""
function measure!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        beta::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        file::IOStream,
        N_sweeps::Int64,
        h::Point3{Float64},
        g::Float64
    )

    # println("measure all...")  # NOTE
    # println("\t T = ", 1./beta)
    # println("\t Js = ", Js)
    # println("\t h = ", h)
    # println("\t g = ", g)
    # println("\t #messure = ", N_sweeps)

    const invN = 1. / sgraph.N_nodes

    # Initialize binning analysis
    E_BA = BinnerA(200)
    Es = Array{Float64}(N_sweeps)

    Mx_BA = BinnerA(200)
    My_BA = BinnerA(200)
    Mz_BA = BinnerA(200)

    M2xy_BA = BinnerA(200)
    M2z_BA = BinnerA(200)
    M2xys = Array{Float64}(N_sweeps)
    M2zs = Array{Float64}(N_sweeps)

    Mquad_BA = BinnerA(200)
    Moct_BA = BinnerA(200)

    Dimer_xy = BinnerA(200)
    Dimer_xy_var = BinnerA(200)
    Dimer_z = BinnerA(200)
    Dimer_z_var = BinnerA(200)

    # "Staggered" magnetizatio
    # Flips the xy and/or z direction if it maps the current system to a ferro-
    # magnetic one. Mostly use for automatic data evaluation
    Nhalf = div(sgraph.N_nodes, 2)
    @inline @inbounds flip1(i::Int64, s::Point3{Float64}) = i <= Nhalf ? s : [-s[1], -s[2], s[3]]
    @inline flip2(i::Int64, s::Point3{Float64}) = i <= Nhalf ? s : -s
    @inline flip3(i::Int64, s::Point3{Float64}) = s
    @inline @inbounds flip4(i::Int64, s::Point3{Float64}) = i <= Nhalf ? s : [s[1], s[2], -s[3]]

    if sign(Js[1][1]) >= 0.0
        if sign(Js[3][1]) >= 0.0
            flip = flip1
        else
            flip = flip2
        end
    elseif sign(Js[3][1]) >= 0.0
        flip = flip3
    else
        flip = flip4
    end

    srMx = 0.0
    srMy = 0.0
    srMz = 0.0

    srdMx = 0.0
    srdMy = 0.0
    srdMz = 0.0

    srMxabs = 0.0
    srMyabs = 0.0
    srMzabs = 0.0

    srdMxabs = 0.0
    srdMyabs = 0.0
    srdMzabs = 0.0

    srM2xy = 0.0
    srM2z = 0.0

    srdM2xy = 0.0
    srdM2z = 0.0

    E_tot = totalEnergy(sgraph, spins, Js, h, g)

    for i in 1:N_sweeps
        E_tot = sweep(sgraph, spins, E_tot, Js, beta, h, g)
        @inbounds Es[i] = E_tot * invN
        push!(E_BA, E_tot * invN)

        S = sum(spins) * invN
        @inbounds push!(Mx_BA, S[1])
        @inbounds push!(My_BA, S[2])
        @inbounds push!(Mz_BA, S[3])

        @inbounds M2xys[i] = sqrt(S[1]^2 + S[2]^2)
        @inbounds M2zs[i] = abs(S[3])
        @inbounds push!(M2xy_BA, M2xys[i])
        @inbounds push!(M2z_BA, abs(S[3]))

        @inbounds temp = mapreduce(v -> Point3{Float64}(0., sqrt(v[1] * v[1] + v[2] * v[2]), abs(v[3])), +, spins) * invN
        @inbounds push!(Mquad_BA, temp[2])
        @inbounds push!(Moct_BA, temp[3])

        # @inbounds map(push!, dimer, get_dimer_parameter(sgraph, spins))

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


        # if additional_observables
        _spins = map(t -> flip(t...), enumerate(spins))
        S = sum(_spins) * invN
        srMx += S[1]
        srMy += S[2]
        srMz += S[3]

        vars = mapreduce(s -> (s - S).^2, +, _spins) * invN
        srdMx += vars[1]
        srdMy += vars[2]
        srdMz += vars[3]

        srM2xy += sqrt(S[1]^2 + S[2]^2)
        srdM2xy += sqrt(vars[1]^2 + vars[2]^2)

        S = mapreduce(v -> abs.(v), +, _spins) * invN
        srMxabs += S[1]
        srMyabs += S[2]
        srMzabs += S[3]

        vars = mapreduce(s -> (abs.(s) - S).^2, +, _spins) * invN
        srdMxabs += vars[1]
        srdMyabs += vars[2]
        srdMzabs += vars[3]
        # end
        yield()
    end
    E_check = totalEnergy(sgraph, spins, Js, h, g)
    if !(E_tot ≈ E_check)
        warn("E_tot diverged by ", (E_check - E_tot) / E_check)
    end

    # This formula is wrong for k =/= 1
    f(x, x2, x3) = (
        (x3 - 3*x*x2 + 2*x^3) * beta^4 * sgraph.N_nodes -
        2 * (x2 - x^2) * beta^3
    ) * sgraph.N_nodes

    cv, dcv = jackknife((x, x2) -> (x2 - x^2) * beta^2 * sgraph.N_nodes, Es, Es.^2)
    dcvdT, ddcvdT = jackknife(f, Es, Es.^2, Es.^3)
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

    #if additional_observables
    K = 1. / N_sweeps
    write_JK!(file, srMx * K, srdMx * K, "rMx  ")
    write_JK!(file, srMy * K, srdMy * K, "rMy  ")
    write_JK!(file, srMz * K, srdMz * K, "rMz  ")
    write_JK!(file, srMxabs * K, srdMxabs * K, "rMxa ")
    write_JK!(file, srMyabs * K, srdMyabs * K, "rMya ")
    write_JK!(file, srMzabs * K, srdMzabs * K, "rMza ")
    write_JK!(file, srM2xy * K, srdM2xy * K, "rMxy ")
    # end

    write_JK!(file, cv, dcv, "cV   ")
    write_JK!(file, dcvdT, ddcvdT, "dcVdT")
    write_JK!(file, dM2xy, ddM2xy, "dM2xy")
    write_JK!(file, dM2z, ddM2z, "dM2z ")

    write_HB!(file, Es_HB, "Energ")
    write_SC!(file, spins, "spins")

    nothing
end


function measure_no_paths!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        beta::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        file::IOStream,
        N_sweeps::Int64,
        h::Point3{Float64},
        g::Float64
    )

    # println("measure no paths...")  # NOTE
    # println("\t T = ", 1./beta)
    # println("\t Js = ", Js)
    # println("\t h = ", h)
    # println("\t g = ", g)
    # println("\t #messure = ", N_sweeps)

    const invN = 1. / sgraph.N_nodes

    E_BA = BinnerA(200)
    Es = Array{Float64}(N_sweeps)

    Mx_BA = BinnerA(200)
    My_BA = BinnerA(200)
    Mz_BA = BinnerA(200)

    M2xy_BA = BinnerA(200)
    M2z_BA = BinnerA(200)
    M2xys = Array{Float64}(N_sweeps)
    M2zs = Array{Float64}(N_sweeps)

    Mquad_BA = BinnerA(200)
    Moct_BA = BinnerA(200)

    # dimer = [BinnerA(200) for _ in eachindex(spins)]
    Dimer_xy = BinnerA(200)
    Dimer_xy_var = BinnerA(200)
    Dimer_z = BinnerA(200)
    Dimer_z_var = BinnerA(200)

    # additional_observables = true
    Nhalf = div(sgraph.N_nodes, 2)
    @inline @inbounds flip1(i::Int64, s::Point3{Float64}) = i <= Nhalf ? s : [-s[1], -s[2], s[3]]
    @inline flip2(i::Int64, s::Point3{Float64}) = i <= Nhalf ? s : -s
    @inline flip3(i::Int64, s::Point3{Float64}) = s
    @inline @inbounds flip4(i::Int64, s::Point3{Float64}) = i <= Nhalf ? s : [s[1], s[2], -s[3]]

    if sign(Js[1][1]) >= 0.0
        if sign(Js[3][1]) >= 0.0
            flip = flip1
        else
            flip = flip2
        end
    elseif sign(Js[3][1]) >= 0.0
        flip = flip3
    else
        flip = flip4
    end

    srMx = 0.0
    srMy = 0.0
    srMz = 0.0

    srdMx = 0.0
    srdMy = 0.0
    srdMz = 0.0

    srMxabs = 0.0
    srMyabs = 0.0
    srMzabs = 0.0

    srdMxabs = 0.0
    srdMyabs = 0.0
    srdMzabs = 0.0

    srM2xy = 0.0
    srM2z = 0.0

    srdM2xy = 0.0
    srdM2z = 0.0

    E_tot = totalEnergy(sgraph, spins, Js, h, g)

    for i in 1:N_sweeps
        E_tot = sweep_no_paths(sgraph, spins, E_tot, Js, beta, h, g)

        @inbounds Es[i] = E_tot * invN
        push!(E_BA, E_tot * invN)

        S = sum(spins) * invN
        @inbounds push!(Mx_BA, S[1])
        @inbounds push!(My_BA, S[2])
        @inbounds push!(Mz_BA, S[3])

        @inbounds M2xys[i] = sqrt(S[1]^2 + S[2]^2)
        @inbounds M2zs[i] = abs(S[3])
        @inbounds push!(M2xy_BA, M2xys[i])
        @inbounds push!(M2z_BA, abs(S[3]))

        @inbounds temp = mapreduce(v -> Point3{Float64}(0., sqrt(v[1] * v[1] + v[2] * v[2]), abs(v[3])), +, spins) * invN
        @inbounds push!(Mquad_BA, temp[2])
        @inbounds push!(Moct_BA, temp[3])

        # @inbounds map(push!, dimer, get_dimer_parameter(sgraph, spins))

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


        # if additional_observables
        _spins = map(t -> flip(t...), enumerate(spins))
        S = sum(_spins) * invN
        srMx += S[1]
        srMy += S[2]
        srMz += S[3]

        vars = mapreduce(s -> (s - S).^2, +, _spins) * invN
        srdMx += vars[1]
        srdMy += vars[2]
        srdMz += vars[3]

        srM2xy += sqrt(S[1]^2 + S[2]^2)
        srdM2xy += sqrt(vars[1]^2 + vars[2]^2)

        S = mapreduce(v -> abs.(v), +, _spins) * invN
        srMxabs += S[1]
        srMyabs += S[2]
        srMzabs += S[3]

        vars = mapreduce(s -> (abs.(s) - S).^2, +, _spins) * invN
        srdMxabs += vars[1]
        srdMyabs += vars[2]
        srdMzabs += vars[3]
        # end
        yield()
    end
    E_check = totalEnergy(sgraph, spins, Js, h, g)
    if !(E_tot ≈ E_check)
        warn("E_tot diverged by ", (E_check - E_tot) / E_check)
    end

    # This formula is wrong for k =/= 1
    f(x, x2, x3) = (
        (x3 - 3*x*x2 + 2*x^3) * beta^4 * sgraph.N_nodes -
        2 * (x2 - x^2) * beta^3
    ) * sgraph.N_nodes

    cv, dcv = jackknife((x, x2) -> (x2 - x^2) * beta^2 * sgraph.N_nodes, Es, Es.^2)
    dcvdT, ddcvdT = jackknife(f, Es, Es.^2, Es.^3)
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

    #if additional_observables
    K = 1. / N_sweeps
    write_JK!(file, srMx * K, srdMx * K, "rMx  ")
    write_JK!(file, srMy * K, srdMy * K, "rMy  ")
    write_JK!(file, srMz * K, srdMz * K, "rMz  ")
    write_JK!(file, srMxabs * K, srdMxabs * K, "rMxa ")
    write_JK!(file, srMyabs * K, srdMyabs * K, "rMya ")
    write_JK!(file, srMzabs * K, srdMzabs * K, "rMza ")
    write_JK!(file, srM2xy * K, srdM2xy * K, "rMxy ")
    # end

    write_JK!(file, cv, dcv, "cV   ")
    write_JK!(file, dcvdT, ddcvdT, "dcVdT")
    write_JK!(file, dM2xy, ddM2xy, "dM2xy")
    write_JK!(file, dM2z, ddM2z, "dM2z ")

    write_HB!(file, Es_HB, "Energ")
    write_SC!(file, spins, "spins")

    nothing
end


function measure!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        beta::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        file::IOStream,
        N_sweeps::Int64,
        h::Point3{Float64}
    )

    # println("measure no g...")  # NOTE
    # println("\t T = ", 1./beta)
    # println("\t Js = ", Js)
    # println("\t h = ", h)
    # println("\t #messure = ", N_sweeps)

    zeroT = beta < 0.0
    const invN = 1. / sgraph.N_nodes

    E_BA = BinnerA(200)
    Es = Array{Float64}(N_sweeps)

    Mx_BA = BinnerA(200)
    My_BA = BinnerA(200)
    Mz_BA = BinnerA(200)

    M2xy_BA = BinnerA(200)
    M2z_BA = BinnerA(200)
    M2xys = Array{Float64}(N_sweeps)
    M2zs = Array{Float64}(N_sweeps)

    Mquad_BA = BinnerA(200)
    Moct_BA = BinnerA(200)

    # dimer = [BinnerA(200) for _ in eachindex(spins)]
    Dimer_xy = BinnerA(200)
    Dimer_xy_var = BinnerA(200)
    Dimer_z = BinnerA(200)
    Dimer_z_var = BinnerA(200)

    # additional_observables = true
    Nhalf = div(sgraph.N_nodes, 2)
    @inline @inbounds flip1(i::Int64, s::Point3{Float64}) = i <= Nhalf ? s : [-s[1], -s[2], s[3]]
    @inline flip2(i::Int64, s::Point3{Float64}) = i <= Nhalf ? s : -s
    @inline flip3(i::Int64, s::Point3{Float64}) = s
    @inline @inbounds flip4(i::Int64, s::Point3{Float64}) = i <= Nhalf ? s : [s[1], s[2], -s[3]]

    if sign(Js[1][1]) >= 0.0
        if sign(Js[3][1]) >= 0.0
            flip = flip1
        else
            flip = flip2
        end
    elseif sign(Js[3][1]) >= 0.0
        flip = flip3
    else
        flip = flip4
    end

    srMx = 0.0
    srMy = 0.0
    srMz = 0.0

    srdMx = 0.0
    srdMy = 0.0
    srdMz = 0.0

    srMxabs = 0.0
    srMyabs = 0.0
    srMzabs = 0.0

    srdMxabs = 0.0
    srdMyabs = 0.0
    srdMzabs = 0.0

    srM2xy = 0.0
    srM2z = 0.0

    srdM2xy = 0.0
    srdM2z = 0.0

    E_tot = totalEnergy(sgraph, spins, Js, h)

    for i in 1:N_sweeps
        if zeroT
            E_tot = sweep(sgraph, spins, E_tot, Js, h)
        else
            E_tot = sweep(sgraph, spins, E_tot, Js, beta, h)
        end
        @inbounds Es[i] = E_tot * invN
        push!(E_BA, E_tot * invN)

        S = sum(spins) * invN
        @inbounds push!(Mx_BA, S[1])
        @inbounds push!(My_BA, S[2])
        @inbounds push!(Mz_BA, S[3])

        @inbounds M2xys[i] = sqrt(S[1]^2 + S[2]^2)
        @inbounds M2zs[i] = abs(S[3])
        @inbounds push!(M2xy_BA, M2xys[i])
        @inbounds push!(M2z_BA, abs(S[3]))

        @inbounds temp = mapreduce(v -> Point3{Float64}(0., sqrt(v[1] * v[1] + v[2] * v[2]), abs(v[3])), +, spins) * invN
        @inbounds push!(Mquad_BA, temp[2])
        @inbounds push!(Moct_BA, temp[3])

        # @inbounds map(push!, dimer, get_dimer_parameter(sgraph, spins))

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


        # if additional_observables
        _spins = map(t -> flip(t...), enumerate(spins))
        S = sum(_spins) * invN
        srMx += S[1]
        srMy += S[2]
        srMz += S[3]

        vars = mapreduce(s -> (s - S).^2, +, _spins) * invN
        srdMx += vars[1]
        srdMy += vars[2]
        srdMz += vars[3]

        srM2xy += sqrt(S[1]^2 + S[2]^2)
        srdM2xy += sqrt(vars[1]^2 + vars[2]^2)

        S = mapreduce(v -> abs.(v), +, _spins) * invN
        srMxabs += S[1]
        srMyabs += S[2]
        srMzabs += S[3]

        vars = mapreduce(s -> (abs.(s) - S).^2, +, _spins) * invN
        srdMxabs += vars[1]
        srdMyabs += vars[2]
        srdMzabs += vars[3]
        # end
        yield()
    end
    E_check = totalEnergy(sgraph, spins, Js, h)
    if !(E_tot ≈ E_check)
        warn("E_tot diverged by ", (E_check - E_tot) / E_check)
    end

    # This formula is wrong for k =/= 1
    f(x, x2, x3) = (
        (x3 - 3*x*x2 + 2*x^3) * beta^4 * sgraph.N_nodes -
        2 * (x2 - x^2) * beta^3
    ) * sgraph.N_nodes

    cv, dcv = jackknife((x, x2) -> (x2 - x^2) * beta^2 * sgraph.N_nodes, Es, Es.^2)
    dcvdT, ddcvdT = jackknife(f, Es, Es.^2, Es.^3)
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

    #if additional_observables
    K = 1. / N_sweeps
    write_JK!(file, srMx * K, srdMx * K, "rMx  ")
    write_JK!(file, srMy * K, srdMy * K, "rMy  ")
    write_JK!(file, srMz * K, srdMz * K, "rMz  ")
    write_JK!(file, srMxabs * K, srdMxabs * K, "rMxa ")
    write_JK!(file, srMyabs * K, srdMyabs * K, "rMya ")
    write_JK!(file, srMzabs * K, srdMzabs * K, "rMza ")
    write_JK!(file, srM2xy * K, srdM2xy * K, "rMxy ")
    # end

    write_JK!(file, cv, dcv, "cV   ")
    write_JK!(file, dcvdT, ddcvdT, "dcVdT")
    write_JK!(file, dM2xy, ddM2xy, "dM2xy")
    write_JK!(file, dM2z, ddM2z, "dM2z ")

    write_HB!(file, Es_HB, "Energ")
    write_SC!(file, spins, "spins")

    nothing
end



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

    println("Simulate...")
    println("\t #therm = ", length(TH_method)) # NOTE
    println("\t T = ", T)
    println("\t Js = ", Js)
    println("\t g = ", g)
    println("\t h = ", h)
    println("\t #messure = ", ME_sweeps)


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
        # Thermalization/Measurement parameters
        TH_sweeps::Int64 = 2_000_000,
        N_switch::Int64 = div(TH_sweeps, 2),
        Freeze_temperature::Float64 = 1.5*maximum(Ts),
        ME_sweeps::Int64 = 5_000_000
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
