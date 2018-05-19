"""
    parallel_tempering!(
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        beta::Float64,
        switch::Int64
    )

Neat?
"""
function parallel_tempering!(
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        beta::Float64,
        switch::Int64
    )

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    comm_size = MPI.Comm_size(comm)

    if (switch + rank) % 2 == 0
        comm_with = rank + 1

        if comm_with < comm_size
            E_tot1 = [E_tot]
            E_tot2 = [0.0]
            MPI.Recv!(E_tot2, comm_with, 0, comm)
            beta1 = [beta]
            beta2 = [0.0]
            MPI.Recv!(beta2, comm_with, 1, comm)

            @fastmath @inbounds dEdT = (E_tot2[1] - E_tot1[1]) *
                (beta2[1] - beta1[1])
            if dEdT > 0.0
                do_swap = true
            elseif exp(dEdT) > rand()
                do_swap = true
            else
                do_swap = false
            end
            MPI.Send([do_swap], comm_with, 2, comm)

            if do_swap
                new_spins = typeof(spins)(length(spins))
                MPI.Recv!(new_spins, comm_with, 3, comm)
                MPI.Send(spins, comm_with, 4, comm)
                MPI.Send(E_tot, comm_with, 5, comm)
                spins = new_spins
                @inbounds E_tot = E_tot2[1]
            end
        end
    else
        comm_with = rank - 1

        if comm_with >= 1
            MPI.Send([E_tot], comm_with, 0, comm)
            MPI.Send([beta], comm_with, 1, comm)

            do_swap = [false]
            MPI.Recv!(do_swap, comm_with, 2, comm)

            if do_swap[1]
                new_spins = typeof(spins)(length(spins))
                E_tot2 = [0.0]
                MPI.Send(spins, comm_with, 3, comm)
                MPI.Recv!(new_spins, comm_with, 4, comm)
                MPI.Recv!(E_tot2, comm_with, 5, comm)
                @inbounds E_tot = E_tot2[1]
                spins = new_spins
            end
        end
    end

    MPI.Barrier(comm)
    return E_tot
end


type Freezer
    sin_values::Vector{Float64}
    exp_values::Vector{Float64}
    exp_deltas::Vector{Float64}

    N::Int64
    N_switch::Int64
    j_step::Float64
    T_max::Float64

    j::Float64
    delta_T::Float64
    T::Float64
    beta::Float64
end

# Constructor, with basic setup
"""
    Freezer(N, T_max[;
        N_switch,
        N_exp_points,
        exp_strength,
        N_sin_points,
        sin_percentage
    ])

A Freezer is an iterator which returns inverse temperatures oscilalting around an
exponentially decaying function starting at some temperature T_max. The final
temperature is given in cool_to(...). The iterator returns a total of N
temperatures, but is constant after N_switch temperatures.

To avoid calling sin() and exp() often, a set number of sin and exp points are
calculated in advance. The number of points can be adjusted through N_exp_points
and N_sin_points. The rate of the exponential decay can be adjusted through
exp_strength (higher = faster decay). The strength of the oscillation depends on
sin_percentage, while the rate of the oscillation is controlled by N_sin_points.
(One period over N_sin_points)
"""
function Freezer(
        N::Int64,
        T_max::Float64;
        N_switch::Int64 = -1,
        N_exp_points::Int64 = 10,
        exp_strength::Float64 = 100.,
        N_sin_points::Int64 = 10_000,
        sin_percentage::Float64 = 0.2
    )

    (N_switch == -1) && (N_switch = div(N, 2))
    sin_values = 1. - sin_percentage .* sin.(linspace(0., 2*pi, N_sin_points))
    exp_values = (exp.(linspace(log(exp_strength + 1), 0., N_exp_points)) - 1.) ./ exp_strength
    exp_deltas = exp_values[2:end] - exp_values[1:end-1]

    Freezer(
        sin_values,
        exp_values,
        exp_deltas,
        N,
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
"""
    cool_to(freezer, T)

Starts the Freezer with a final temperature T.

Example:
    f = Freezer(1_000_000, 1.0)
    for beta in cool_to(f, 0.1)
        ...
    end
"""
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

# next temperature
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


"""
    BinnerA(min_output_size)

Creates a new Binning Analysis which keeps min_output_size values around after
binning. Returns a Binning Analysis object. Use push! to add values.
"""
function BinnerA(min_output_size::Integer)
    BinnerA(
        Compressor[],
        Array{Float64}(2 * min_output_size),
        UInt32(1),
        UInt32(2 * min_output_size),
        Int64[],
        Float64[],
        Float64[]
    )
end


"""
    push!(BinningAnalysis, value)

Pushes a new value into the Binning Analysis.
"""
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

    return nothing
end


# recursion, back-end function
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
"""
    var(BinningAnalysis[, lvl])

Calculates the **variance/N** of a given level in the Binning Analysis.
"""
function var(B::BinnerA, lvl::Int64=1_000_000)
    # The upper most lvl needs to be calculated explicitly
    if lvl < length(B.count)
        return ((B.x2_sum[lvl+1] / B.count[lvl+1]) - (B.x_sum[lvl+1] / B.count[lvl+1])^2) / (B.count[lvl+1] - 1)
    else
        return ((sum(B.output[1:B.i-1].^2) / (B.i-1)) - (sum(B.output[1:B.i-1]) / (B.i-1))^2) / (B.i - 2)
    end
end


"""
    mean(BinningAnalysis[, lvl])

Calculates the mean for a given level in the Binning Analysis
"""
function mean(B::BinnerA, lvl::Int64=1_000_000)
    if lvl < length(B.count)
        return B.x_sum[lvl] / B.count[lvl]
    else
        return mean(B.output[1:B.i-1])
    end
end


"""
    tau(BinningAnalysis, lvl)

Calculates tau...
"""
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

    """
        BinnerH(bin_width)

    Creates a new Histogram Binner with a given bin width.
    """
    BinnerH(bin_width::Float64) = new(Dict{Int64, UInt32}(), bin_width)
end

"""
    push!(HistogramBinner, x)

Pushes x into the Histogram Binner.
"""
function push!(B::BinnerH, x::Float64)
    id = floor(Int64, x / B.bin_width)

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
"""
    jackknife(function, args...)

Performs a jackknife analysis given a function and its arguments. Returns the
mean and standard error.

Example:
    xs = rand(1000)
    ys = rand(1000)
    value, error = jackknife((x, y) -> x*y, xs, ys)
"""
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
