
################################################################################
#### New Simulation stack
################################################################################

# change: N/2 -> var N_switch
# Freezer is a Generator for slow cooling iterators. For the first half the
# temperatures follow an exponential decay from T_max to T with sin-noise. In
# the second half the iterator yields T. Usuage:

# Freezer(Number_of_points, T_max;
#       N_exp_points = number of points to sample exp from,
#       exp_strength = result of exp that represents T_max (higher -> faster decay),
#       N_sin_points = number of sin-points used for a period. (more -> slower sin),
#       sin_percentage = sin is multiplied with this factor)
#
# Freezer will yield a type which acts as a Generator (and later Iterator).

# cool_to(Freezer, T)
# will set up Freezer as an Iterator for final temperature T.
# You can do:
#   for T in cool_to(Freezer, T); ... end
# to get your Temperatures

type Freezer
    sin_values::Vector{Float64}
    exp_values::Vector{Float64}
    exp_deltas::Vector{Float64}

    N::Int64
    # N_half::Int64
    N_switch::Int64
    j_step::Float64
    T_max::Float64

    j::Float64
    delta_T::Float64
    T::Float64
    beta::Float64
end

# Constructor, with basic setup
function Freezer(N::Int64, T_max::Float64;
        N_switch::Int64=-1,
        N_exp_points::Int64=10, exp_strength::Float64=100.,
        N_sin_points::Int64=10_000, sin_percentage::Float64=0.2)

    (N_switch == -1) && (N_switch = div(N, 2))
    sin_values = 1. - sin_percentage .* sin(linspace(0., 2*pi, N_sin_points))
    exp_values = (exp(linspace(log(exp_strength + 1), 0., N_exp_points)) - 1.) ./ exp_strength
    exp_deltas = exp_values[2:end] - exp_values[1:end-1]
    # N_half = div(N, 2)

    Freezer(
        sin_values,
        exp_values,
        exp_deltas,
        N,
        # N_half,
        N_switch,
        (N_exp_points - 1) / N_switch,
        T_max,
        0.,
        0.,
        0.,
        0.
    )
end

# Sets up T as the final temperature
function cool_to(F::Freezer, T::Float64)
    F.delta_T = F.T_max - T
    F.T = T
    F.beta = 1./T
    F.j = 1. - F.j_step
    return F
end

# Iterator functions
start(F::Freezer) = 0
done(F::Freezer, i::Int64) = i >= F.N
length(F::Freezer) = F.N
eltype(::Freezer) = Float64

@inline function next(F::Freezer, i::Int64)
    if i < F.N_switch # F.N_half
        i += 1
        F.j += F.j_step
        return 1. / (
            F.delta_T * (
                F.exp_values[floor(Int64, F.j)] +
                (F.j - floor(F.j)) * F.exp_deltas[floor(Int64, F.j)]
            ) * F.sin_values[(i - 1) % length(F.sin_values) + 1] + F.T
        ), i

    else
        i += 1
        return F.beta, i
    end
end

# "Shock"-freezing at constant temeprature
# init: F = ConstantT(N); cool_to(F, T)
immutable ConstantT
    beta::Float64
    N::Int64

    ConstantT(T::Float64, N::Int64) = new(1./T, N)
    ConstantT(N::Int64) = new(0., N)
end

cool_to(F::ConstantT, T::Float64) = ConstantT(T, F.N)

start(F::ConstantT) = 0
next(F::ConstantT, i::Int64) = (F.beta, i+1)
done(F::ConstantT, i::Int64) = i >= F.N
eltype(::ConstantT) = Float64
length(F::ConstantT) = F.N

# Iterator end


################################################################################
#### Binning Anaylsis
################################################################################


# This is a Node in the Binning Analysis tree, that averages two values. There
# is one of these for each binning level. Since when two values should be
# compressed, this is done immediately, so that only one value needs to be saved.
# switch indicates whether value should be written to or averaging should happen.

# TODO:     value -> Type T?
type Compressor
    value::Float64
    switch::UInt8
end


type BinnerA
    # list of Compressors, one per level
    compressors::Vector{Compressor}

    # currently binned values
    output::Vector{Float64}
    # output array has a constant size. i indicates the first index that can
    # be (over-)written
    i::UInt32
    # max length of output
    max_output_size::UInt32

    # sum(x) for all values on a given lvl (excluding the upper most level)
    x_sum::Vector{Float64}
    # sum(x.^2) for all values on a given lvl (excluding the upper most level)
    x2_sum::Vector{Float64}
    # number of values that are summed on a given lvl (excluding the upper most level)
    count::Vector{Int64}
end


# Constructor
# min_output_size is the minimum number of points used for the upper most level
function BinnerA(min_output_size::Integer)
    BinnerA(
        Compressor[],
        Array(Float64, 2 * min_output_size),
        UInt32(1), UInt32(2 * min_output_size),
        Int64[], Float64[], Float64[]
    )
end


function push!(B::BinnerA, value::Float64)
    # first set of values
    if isempty(B.compressors)
        # fill output if possible
        if B.i < B.max_output_size
            B.output[B.i] = value
            B.i += UInt32(1)
            return nothing

        # once output is full, add compressor and drop one level
        else
            push!(B.compressors, Compressor(0., UInt8(0)))
            push!(B.x_sum, 0.)
            push!(B.x2_sum, 0.)
            push!(B.count, 0)

            B.output[B.i] = value
            B.i = 1
            for x in B.output
                push!(B, 1, x)
            end

            return nothing
        end

    else
        # push value in compressor stack
        if B.i <= B.max_output_size
            push!(B, 1, value)
            return nothing

        # if output is full, add new level and drop down
        else
            push!(B.compressors, Compressor(0., UInt8(0)))
            push!(B.x_sum, 0.)
            push!(B.x2_sum, 0.)
            push!(B.count, 0)

            B.i = 1
            for x in B.output
                push!(B, length(B.compressors), x)
            end
            push!(B, 1, value)

            return nothing
        end
    end

    println("How did you get here? user push!")
    return nothing
end


# recursion
function push!(B::BinnerA, lvl::Int64, value::Float64)
    C = B.compressors[lvl]

    # any value propagating through this function is new to lvl. Therefore we
    # add it to the sums. Note that values pushed to the output arrays are not
    # added here until the array drops to the next level. (New compressors are
    # added)
    B.x_sum[lvl] += value
    B.x2_sum[lvl] += value^2
    B.count[lvl] += 1

    if C.switch == 0
        # Compressor has space -> save value
        C.value = value
        C.switch = 1
        return nothing
    else
        # Do averaging
        if lvl == length(B.compressors)
            # No more propagation possible -> push to output array
            B.output[B.i] = 0.5 * (C.value + value)
            B.i += 1
            C.switch = 0
            return nothing
        else
            # propagate to next lvl
            C.switch = 0
            push!(B, lvl+1, 0.5 * (C.value + value))
            return nothing
        end
    end

    println("How did you get here? recursive push!")
    return nothing
end


# NOTE:
# This is actually Var(x)/N
# This definition is more practical for the binning analysis and error
# calculation since we need the standard error sqrt(Var(x) / N)
function var(B::BinnerA, lvl::Int64)
    # The upper most lvl needs to be calculated explicitly
    if lvl < length(B.count)
        return ((B.x2_sum[lvl+1] / B.count[lvl+1]) - (B.x_sum[lvl+1] / B.count[lvl+1])^2) / (B.count[lvl+1] - 1)
    else
        return ((sum(B.output[1:B.i-1].^2) / (B.i-1)) - (sum(B.output[1:B.i-1]) / (B.i-1))^2) / (B.i - 2)
    end
end


function tau(B::BinnerA, lvl::Int64)
    var_0 = var(B, 0)
    var_l = var(B, lvl)
    0.5 * (var_l / var_0 - 1)
end


################################################################################
#### Histogram binning
################################################################################

# This is kept very simplistic. The width of bins is fixed and bins are created
# when needed. bin width should be somewhere in the range of 0.001 - 0.0001.
type BinnerH
    data::Dict{Int64, Int64}
    bin_width::Float64

    BinnerH(bin_width::Float64) = new(Dict{Int64, UInt32}(), bin_width)
end

function push!(B::BinnerH, x::Float64)
    id = floor(Int64, x / B.bin_width)#Int64(div(x, B.bin_width))

    if haskey(B.data, id)
        B.data[id] += 1
    else
        push!(B.data, id => 1)
    end
end


################################################################################
#### Jackknife
################################################################################


# This requires the full xs arrays :(
function jackknife(f::Function, args...)
    N = length(args[1])
    x0s = map(sum, args) # (xs, ...) -> (sum xs, ...)
    x_avs_s = map((xs, x0) -> [(x0 - xs[i]) / (N-1) for i in 1:N], args, x0s)
    x0s = map(x -> x / N, x0s)

    y0 = map(f, x0s...)
    y_avs = map(f, x_avs_s...)

    # outputs mean and standard deviation of the distribution of means (which
    # is the standard error of ys)
    y0, sqrt((N-1) / N * sum((y_avs - y0).^2))
end


################################################################################
#### File handling
################################################################################


# What do we have?
#=
# metadata
N_points
TH_sweeps
TH_Temp
ME_sweeps
sys_size
N_nodes
K_edges
Js

# Binning Analysis E

=#

# Do other things in a modular fashion
#=
tag
required metadata (to determine start and end)
data
=#

# file writers
function write_header!(file::IOStream,
        N_points::Int64, TH_sweeps::Int64, TH_Temp::Float64,
        ME_sweeps::Int64, sys_size::Int64, N_nodes::Int64, K_edges::Int64,
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        h::Point3{Float64}, T::Float64)

    write(file, "V02")
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
    write(file, T)

    nothing
end


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
    write(file, sum(B.output))

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

function get_dimer_parameter(sgraph::SGraph, spins::Vector{Point3{Float64}})
    return map(eachindex(sgraph.nodes)) do i
        mapreduce(e -> e.xy + e.z, +, sgraph.nodes[i].first)
    end
end


function measure!(sgraph::SGraph, spins::Vector{Point3{Float64}}, beta::Float64,
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        file::IOStream, N_sweeps::Int64=1000, h::Point3{Float64}=Point3(0.))

    zeroT = beta < 0.0

    E_BA = BinnerA(200)
    Es = Array(Float64, N_sweeps)

    Mx_BA = BinnerA(200)
    My_BA = BinnerA(200)
    Mz_BA = BinnerA(200)

    M2xy_BA = BinnerA(200)
    M2z_BA = BinnerA(200)
    M2xys = Array(Float64, N_sweeps)
    M2zs = Array(Float64, N_sweeps)

    Mquad_BA = BinnerA(200)
    Moct_BA = BinnerA(200)

    dimer = [BinnerA(200) for _ in eachindex(spins)]

    additional_observables = true
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
        b2 = Base.gc_bytes()
        @inbounds Es[i] = E_tot / sgraph.N_nodes
        push!(E_BA, E_tot / sim.N_nodes)

        S = sum(spins) / sgraph.N_nodes
        @inbounds push!(Mx_BA, S[1])
        @inbounds push!(My_BA, S[2])
        @inbounds push!(Mz_BA, S[3])

        @inbounds M2xys[i] = sqrt(S[1]^2 + S[2]^2)
        @inbounds M2zs[i] = abs(S[3])
        @inbounds push!(M2xy_BA, M2xys[i])
        @inbounds push!(M2z_BA, abs(S[3]))

        @inbounds temp = mapreduce(v -> Point3{Float64}(0., sqrt(v[1] * v[1] + v[2] * v[2]), abs(v[3])), +, spins) / sgraph.N_nodes
        @inbounds push!(Mquad_BA, temp[2])
        @inbounds push!(Moct_BA, temp[3])

        @inbounds map(push!, dimer, get_dimer_parameter(sgraph, spins))

        if additional_observables
            _spins = map(t -> flip(t...), enumerate(spins))
            S = sum(_spins) / sgraph.N_nodes
            srMx += S[1]
            srMy += S[2]
            srMz += S[3]

            vars = mapreduce(s -> (s - S).^2, +, _spins) / sgraph.N_nodes
            srdMx += vars[1]
            srdMy += vars[2]
            srdMz += vars[3]

            srM2xy += sqrt(S[1]^2 + S[2]^2)
            srdM2xy += sqrt(vars[1]^2 + vars[2]^2)

            S = mapreduce(abs, +, _spins) / sgraph.N_nodes
            srMxabs += S[1]
            srMyabs += S[2]
            srMzabs += S[3]

            vars = mapreduce(s -> (abs(s) - S).^2, +, _spins) / sgraph.N_nodes
            srdMxabs += vars[1]
            srdMyabs += vars[2]
            srdMzabs += vars[3]
        end
        yield()
    end
    # This formula is wrong for k =/= 1
    f(x, x2, x3) = (
        (x3 - 3*x*x2 + 2*x^3) * beta^4 * sgraph.N_nodes -
        2 * (x2 - x^2) * beta^3
        # 2 * (x2 - x^2) * beta^3 +
        # (x3 - 3*x*x2 + 2*x^3) * beta^4 * sgraph.N_nodes
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
    for i in eachindex(dimer)
        write_BA!(file, dimer[i], rpad(string(i), 5))
    end

    if additional_observables
        K = 1 / N_sweeps
        write_JK!(file, srMx * K, srdMx * K, "rMx  ")
        write_JK!(file, srMy * K, srdMy * K, "rMy  ")
        write_JK!(file, srMz * K, srdMz * K, "rMz  ")
        write_JK!(file, srMxabs * K, srdMxabs * K, "rMxa ")
        write_JK!(file, srMyabs * K, srdMyabs * K, "rMya ")
        write_JK!(file, srMzabs * K, srdMzabs * K, "rMza ")
        write_JK!(file, srM2xy * K, srdM2xy * K, "rMxy ")
    end

    write_JK!(file, cv, dcv, "cV   ")
    write_JK!(file, dcvdT, ddcvdT, "dcVdT")
    write_JK!(file, dM2xy, ddM2xy, "dM2xy")
    write_JK!(file, dM2z, ddM2z, "dM2z ")

    write_HB!(file, Es_HB, "Energ")
    write_SC!(file, spins, "spins")

    nothing
end


# Single Temperature simulate
function simulate!(sgraph::SGraph, spins::Vector{Point3{Float64}}, sys_size::Int64,
        path::String, filename::String,
        T::Float64, Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        TH_method::Freezer, ME_sweeps::Int64, h::Point3{Float64}=Point3(0.))


    # Fool-proof file creation that was actually not fool-proof
    if isdir(path)
        println(
            path, " already exists",
            ". New files will be inserted."
        )
        #folder_name = folder_name * string(time())
    else
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
        Int64(sgraph.N_nodes), sgraph.K_edges, Js, h, T
    )

    # tic()
    # Thermalization
    # spins = rand_spin(sgraph.N_nodes)
    init_edges!(sgraph, spins)

    if T > 0.0
        db_exp = typemax(Int64)
        for beta in cool_to(TH_method, T) #1:TH_sweeps
            b0 = Base.gc_bytes()
            sweep(sgraph, spins, Js, beta, h)
            db = Base.gc_bytes() - b0
            if db < db_exp
                db_exp = db
            elseif db > db_exp
                println("Therm. ", db)
           end
        end

        beta = 1. / T
    else
        ME_sweeps += TH_sweeps
        beta = -1.0
    end
    # Measurement
    # Es = Array(Float64, ME_sweeps)
    # Ms, Ms2, Ms_abs, Ms_abs2 =
    measure!(sgraph, spins, beta, Js, file, ME_sweeps, h)

    # # saving
    # dump_Floats!(file, Es)
    # map(x -> write(file, x), Ms)
    # map(x -> write(file, x), Ms_abs)
    # dump_Points!(file, spins)
    # map(x -> write(file, x), Ms2)
    # map(x -> write(file, x), Ms_abs2)

    close(file)
    # println("T = ", T, " done. Time: ", toq())

    nothing
end



function simulate!(sgraph::SGraph, spins::Vector{Point3{Float64}}, sys_size::Int64,
        path::String, filename::String,
        Ts::Vector{Float64},
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        TH_method::Union{ConstantT, Freezer}, ME_sweeps::Int64,
        h::Point3{Float64}=Point3(0.)) #,
        #do_parallel::Bool=true)


    # if do_parallel
    #     pmap(
    #         (i, T) -> simulate!(
    #             deepcopy(sgraph), sys_size,
    #             path * folder_name * "/", filename * string(i),
    #             T, deepcopy(Js),
    #             deepcopy(TH_method), ME_sweeps
    #         ),
    #         enumerate(Ts)
    #     )
    # else
    for (i, T) in enumerate(Ts)
        simulate!(
            sgraph, sys_size,
            path, filename * string(i),
            T, Js,
            TH_method, ME_sweeps, h
        )
    end
    # end

    nothing
end