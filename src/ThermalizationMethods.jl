#=
            AbstractThermalizationMethod
AbstractTemperatureGenerator        AbstractParallelTemperingAlgorithm
- ConstantT                         - ParallelTempering
- Freezer                           - ProbabilityEqualizer
                                    - TimeEqualizer
                                    - EnergyRebalancer
.

each one should work as

    init_edges!(...)
    E_tot = totalEnergy(...)
    state = initialize(thermalizer, ...)            # or start()
    while !done(thermalizer, state)
        state = next(thermalizer, state, E_tot)
        E_tot = sweep(...)
    end
    beta = last(thermalizer, state)
=#


# TODO
#=
Every adaptive algorithm has to thermalize thoroughly first
- Preferably PT + simulated annealing
- lot's of sweeps (same order as no PT!?)
- explicitly checking/recording energy equillibration might be useful

Energy balancing
1) long thermalization, PT + SA
2) measure E_tot quickly
3) shift T down slightly
4) thermalize quickly
5) goto 2) until enough points
6) calculate dE/dT, apply some smoothness constraint
7) concentrate Ts where |dE/dT| is large
8) medium thermalization
    if T_new - T_old > 0 do constant T
    else do SA

Ts:
|   |   |   |   |   ||
(long thermalization)
|   |   |   |   |   ||  quick measure
|  |   |   |   |   | |  quick thermalize + measure
| |   |   |   |   |  |  quick thermalize + measure
||   |   |   |   |   |  quick thermalize + measure
|  |  | | |   |      |  concetrate Ts
(medium thermalization)
|  |  | | |   |      |  output

=#

#=
# TODO
Where should the thermalization method be created?

1) in main file, then given to simulate!
- can be built in mainf
- would clean up stuff

How should Tgen be created for pt methods?
=#

"""
    (<: AbstractThermalizationMethod)(; kwargs...)

Creates a thermalization method.

# Keyword Arguments
- `N::Int64`: Number of thermalization sweeps
- `T_max::Float64`: The initial/maximum temperature
- `N_switch::Int64`: Number of sweeps before Freezer switches to constant T
- `batch_size::Int64`: Number of sweeps between each parallel tempering attempt
- `adaptive_sample_size::Int64`: Number of sweeps to consider for adaptive step

# Methods
- `ConstantT(N)`
- `Freezer(N, T_max, N_switch, ...)`: Implements simulated annealing
- `ParallelTempering()`
- `ProbabilityEqualizer()`
"""
abstract type AbstractThermalizationMethod end
abstract type AbstractTemperatureGenerator <: AbstractThermalizationMethod end
abstract type AbstractParallelTemperingAlgorithm <: AbstractThermalizationMethod end

"""
    is_parallel(<:AbstractThermalizationMethod)

Returns true if the thermalization method uses parallelism. In this case
parallel tempering should be used for the measurement.
"""
is_parallel(::AbstractTemperatureGenerator) = false
is_parallel(::AbstractParallelTemperingAlgorithm) = true

"""
    initialize(<:AbstractThermalizationMethod, T::Float64)

Returns the initial state of a thermalization method.
"""
function initialize() end

"""
    next(<:AbstractThermalizationMethod, state)

Returns the next set of variables as well as the next state.
"""
function next() end

"""
    current_index(<:AbstractThermalizationMethod, state)

Returns the current index of the iteration.
"""
function current_index() end

"""
    done(<:AbstractThermalizationMethod, state)

Returns true if the thermalization method has reached its last output.
"""
function done() end

"""
    last(<:AbstractThermalizationMethod, state)

Returns the last set of variables of a thermalization process.
"""
function last() end

"""
    lenght(<:AbstractThermalizationMethod)

Returns the length of a thermalization process.
"""
function length() end


batch_size(::AbstractThermalizationMethod) = -1
adaptive_sample_size(::AbstractThermalizationMethod) = -1

# Some methods, like ConstantT, don't have a "T_max"
T_max(::AbstractThermalizationMethod) = -1.0

#######################


struct ConstantT <: AbstractTemperatureGenerator
    N::Int64
end
ConstantT(; N::Int64=0, kwargs...) = ConstantT(N)

initialize(th::ConstantT, T::Float64) = 1.0 / T, (1.0 / T, 0)
next(::ConstantT, state::Tuple{Float64, Int64}) = (state[1], (state[1], state[2]+1))
current_index(::ConstantT, state::Tuple{Float64, Int64}) = state[2]
done(th::ConstantT, state::Tuple{Float64, Int64}) = state[2] >= th.N
last(::ConstantT, state::Tuple{Float64, Int64}) = state[1]
length(th::ConstantT) = th.N


#######################


struct Freezer <: AbstractTemperatureGenerator
    N::Int64
    N_switch::Int64
    N_per_exp::Int64
    T_max::Float64

    sin_values::Vector{Float64}
    exp_values::Vector{Float64}
    exp_deltas::Vector{Float64}
end

function Freezer(;
        N::Int64 = 0,
        T_max::Float64 = 10.0,
        N_switch::Int64 = div(N, 2),
        N_exp_points::Int64 = 10,
        exp_strength::Float64 = 100.,
        N_sin_points::Int64 = 10_000,
        sin_percentage::Float64 = 0.2,
        kwargs...
    )
    sin_values = 1. .- sin_percentage * sin.(
        range(0., stop=2.0pi, length=N_sin_points)
    )
    exp_values = (exp.(
        range(log(exp_strength + 1.), stop=0., length=N_exp_points)
    ) .- 1) / exp_strength
    exp_deltas = exp_values[2:end] - exp_values[1:end-1]

    Freezer(
        N,
        N_switch,
        ceil((N_switch+1) / length(exp_deltas)),
        T_max,
        sin_values,
        exp_values,
        exp_deltas
    )
end

function initialize(th::Freezer, T_min::Float64)
    1.0 / th.T_max, (
        0,
        1,
        1,
        th.T_max - T_min,
        T_min,
        1. / T_min
    )
end

@inline function next(th::Freezer, state)
    i, j, k, dT, T_final, beta = state
    if i < th.N_switch
        # L = length(th.exp_values)
        factor = (i - th.N_per_exp * (j-1)) / th.N_per_exp
        exp_value = th.exp_values[j] + factor * th.exp_deltas[j]
        T = T_final + dT * exp_value * th.sin_values[k]
        return 1. / T, (
            i+1,
            i < j*th.N_per_exp ? j : j+1,
            k < length(th.sin_values) ? k+1 : 1,
            dT,
            T_final,
            beta
        )
    else
        return beta, (i+1, j, k, dT, T_final, beta)
    end
end

current_index(::Freezer, state::Tuple) = state[1]
done(th::Freezer, state::Tuple) = state[1] >= th.N
last(::Freezer, state::Tuple) = state[end]
length(th::Freezer) = th.N
T_max(th::Freezer) = th.T_max


##############################################
#=
Stacking parallel tempering and simulated annealing should be
fine, so let's allow this.

All adaptive methods share code with the basic pt algorithm,
so they should just spawm a pt...
- Nope, sometimes I need probabilities, sometimes not
=#
##############################################


# This is a dummy to construct TGen without a ParallelTemperingAlgirthm
struct NoParallelTempering{
        TGen <: AbstractTemperatureGenerator
    } <: AbstractParallelTemperingAlgorithm
    x::Int64
end
function (::Type{NoParallelTempering{TGen}})(;
        kwargs...
    ) where TGen <: AbstractTemperatureGenerator
    TGen(; kwargs...)
end


#######################


struct ParallelTempering{
        ATG <: AbstractTemperatureGenerator
    }  <: AbstractParallelTemperingAlgorithm

    Tgen::ATG
    batch_size::Int64
end

function (::Type{ParallelTempering{TGen}})(;
        batch_size::Int64 = 10,
        kwargs...
    )  where TGen <: AbstractTemperatureGenerator
    ParallelTempering(
        TGen(; kwargs...),
        batch_size
    )
end

function initialize(
        th::ParallelTempering,
        T::Float64,
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        # E_tot::Float64
    )
    beta, Tgen_state = initialize(th.Tgen, T)
    return beta, (Tgen_state, sgraph, spins, 0)
end

function next(
        th::ParallelTempering,
        state::Tuple,
        E_tot::Float64
    )
    Tgen_state, sgraph, spins = state
    beta, Tgen_state = next(th.Tgen, Tgen_state)
    i = current_index(th.Tgen, Tgen_state)

    if i % th.batch_size == 0
        E_tot = _parallel_tempering!(
            sgraph, spins, E_tot, beta
        )
        __switch__[] = 1 - __switch__[]
    end

    return beta, E_tot, (
        Tgen_state, sgraph, spins
    )
end

current_index(th::ParallelTempering, state::Tuple) = current_index(th.Tgen, state[1])
done(th::ParallelTempering, state::Tuple) = done(th.Tgen, state[1])
last(th::ParallelTempering, state::Tuple) = last(th.Tgen, state[1])
length(th::ParallelTempering) = length(th.Tgen)
batch_size(th::ParallelTempering) = th.batch_size
T_max(th::ParallelTempering) = T_max(th.Tgen)


#######################


const __p__ = Channel{Float64}(1)
set_p!(value) = put!(__p__, value)
const __beta__ = Channel{Float64}(1)
set_beta!(value) = put!(__beta__, value)


struct ProbabilityEqualizer{
        ATG <: AbstractTemperatureGenerator
    } <: AbstractParallelTemperingAlgorithm
    # PT::ParallelTempering{<:AbstractTemperatureGenerator}
    Tgen::ATG
    batch_size::Int64
    skip::Int64
    M::Int64
    adaptive_sample_size::Int64
end

# function ProbabilityEqualizer(
#         Tgen::ATG,
#         batch_size::Int64,
#         adaptive_sample_size::Int64
#     ) where {ATG <: AbstractTemperatureGenerator}
#
#     return ProbabilityEqualizer(
#         Tgen,
#         batch_size,
#         div(length(Tgen), adaptive_sample_size),
#         adaptive_sample_size
#     )
# end

function (::Type{ProbabilityEqualizer{TGen}})(;
        skip::Int64 = -1,
        batch_size::Int64 = 10,
        adaptive_sample_size::Int64 = 100batch_size,
        kwargs...
    )  where TGen <: AbstractTemperatureGenerator
    Tgen = TGen(; kwargs...)
    if skip == -1
        skip = div(length(Tgen), 2)
    end
    skip = div(skip, adaptive_sample_size)
    ProbabilityEqualizer(
        Tgen,
        batch_size,
        skip,
        div(length(Tgen), adaptive_sample_size) - skip,
        adaptive_sample_size
    )
end

#
function initialize(
        th::ProbabilityEqualizer,
        T::Float64,
        sgraph::SGraph,
        spins::Vector{Point3{Float64}}
    )
    cum_prob = 0.0
    N_prob = 0
    beta, Tgen_state = initialize(th.Tgen, T)

    return beta, (Tgen_state, cum_prob, N_prob, sgraph, spins)
end

function next(th::ProbabilityEqualizer, state::Tuple, E_tot::Float64)
    # comm = MPI.COMM_WORLD
    # comm_size = MPI.Comm_size(comm)
    # comm_rank = MPI.Comm_rank(comm)
    comm_size = nprocs()
    comm_rank = myid()

    Tgen_state, cum_prob, N_prob, sgraph, spins = state
    beta, Tgen_state = next(th.Tgen, Tgen_state)
    i = current_index(th.Tgen, Tgen_state)

    if i % th.batch_size == 0
        E_tot, p = _parallel_tempering_adaptive!(
            sgraph, spins, E_tot, beta
        )
        if p >= 0.0
            cum_prob += p
            N_prob += 1
        end
        __switch__[] = 1 - __switch__[]
    end

    if i % th.adaptive_sample_size == 0
        k = div(i, th.adaptive_sample_size) - th.skip

        if k > 0
            # MPI.Barrier(comm)
            p = N_prob == 0 ? 0.0 : cum_prob / N_prob
            # TODO
            # probs = MPI.Allgather(p, comm)
            # temperatures = MPI.Allgather(1.0 / beta, comm)

            if (comm_rank != 0) && (comm_rank != comm_size-1)
                prob_sum = sum(probs)
                norm_prob = probs / prob_sum
                mean_prob = mean(norm_prob)
                down = (norm_prob[comm_rank] - mean_prob) *
                    (temperatures[comm_rank+1] - temperatures[comm_rank])
                up = (norm_prob[comm_rank+1] - mean_prob) *
                    (temperatures[comm_rank+1] - temperatures[comm_rank+2])
                # down = 0.01(log(norm_prob[rank]) - log(mean_prob)) *
                #     (temperatures[rank+1] - temperatures[rank])
                # up = 0.01(log(norm_prob[rank+1]) - log(mean_prob)) *
                #     (temperatures[rank+1] - temperatures[rank+2])
                new_T = max(
                    temperatures[comm_rank+1] + (1.0 - k / th.M) * (down + up),
                    temperatures[comm_rank] + 0.01
                )

                beta = 1.0 / new_T

                # NOTE
                # sleep(0.01 * comm_rank)
                # print("[$comm_rank] T = $(round(temperatures[comm_rank+1], 3))   ->  $(round(new_T, 3)) \t\t down = $down \t up = $up \t\t proabilities = $(probs)\n")
            end
        end

        cum_prob = 0.0
        N_prob = 0
    end

    return beta, E_tot, (Tgen_state, cum_prob, N_prob, sgraph, spins)
end

current_index(th::ProbabilityEqualizer, state) = current_index(th.Tgen, state[1])
done(th::ProbabilityEqualizer, state) = done(th.Tgen, state[1])
last(th::ProbabilityEqualizer, state) = last(th.Tgen, state[1])
length(th::ProbabilityEqualizer) = length(th.Tgen)
batch_size(th::ProbabilityEqualizer) = th.batch_size
adaptive_sample_size(th::ProbabilityEqualizer) = th.adaptive_sample_size
T_max(th::ProbabilityEqualizer) = T_max(th.Tgen)


################################################################################

# MPI-ish style interprocess communication
# Non-blocking send via remotecall(set..., comm_with, data)
# Blocking recv via take!(__...__)

# to sync up processes
# left and right communicators should fix race conditions...
const __left_ready__ = Channel{Bool}(1)
const __right_ready__ = Channel{Bool}(1)
left_is_ready() = put!(__left_ready__, true)
right_is_ready() = put!(__right_ready__, true)

# data transfer
const __E_tot__ = Channel{Float64}(1)
set_E_tot!(value) = put!(__E_tot__, value) # println("$(myid()):  $value")
const __do_swap__ = Channel{Bool}(1)
set_do_swap!(value) = put!(__do_swap__, value)
const __spins__ = Channel{Vector{Point3{Float64}}}(1)
set_spins!(value) = put!(__spins__, value)

# global switch allows continuation of parallel tempering after thermalization
const __switch__ = Ref{Int64}(0)


"""
    parallel_tempering!(
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        beta::Float64,
        switch::Int64
    )

Neat?
"""
function _parallel_tempering!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        beta::Float64
    )
    rank = myid()
    comm_size = nprocs()

    if (__switch__[] + rank) % 2 == 0
        # me -> comm_with
        comm_with = rank + 1

        if comm_with <= comm_size
            # Tell comm_with "I am ready"
            remotecall(left_is_ready, comm_with)
            # Wait for comm_with to be ready
            take!(__right_ready__)

            # Receive data from comm_with
            remote_E_tot = take!(__E_tot__)
            remote_beta = take!(__beta__)

            # Determine whether to swap and notify comm_with
            @fastmath dEdT = (remote_E_tot - E_tot) * (remote_beta - beta)
            do_swap = (dEdT > 0.0) || (exp(dEdT) > rand())
            remotecall(set_do_swap!, comm_with, do_swap)

            # Perform data swap
            if do_swap
                old_spins = deepcopy(spins)
                remotecall(set_spins!, comm_with, old_spins)
                remotecall(set_E_tot!, comm_with, E_tot)
                spins .= take!(__spins__)
                init_edges!(sgraph, spins)
                return remote_E_tot
            end
        end
        return E_tot
    else
        # comm_with <- me
        comm_with = rank - 1

        if comm_with > 0
            # Tell comm_with "I am ready"
            remotecall(right_is_ready, comm_with)
            # Wait for comm_with to be ready
            take!(__left_ready__)

            # Send data to comm_with
            remotecall(set_E_tot!, comm_with, E_tot)
            remotecall(set_beta!, comm_with, beta)

            # Wait for response (whether to swap)
            do_swap = take!(__do_swap__)

            if do_swap
                old_spins = deepcopy(spins)
                remotecall(set_spins!, comm_with, old_spins)
                spins .= take!(__spins__)
                remote_E_tot = take!(__E_tot__)
                init_edges!(sgraph, spins)
                return remote_E_tot
            end
        end
        return E_tot
    end
    nothing
end


# TODO
# p needs to be transfered
# this may interfere with something in ProbabilityEqualizer
function _parallel_tempering_adaptive!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        beta::Float64
    )
    rank = myid()
    comm_size = nprocs()
    p = -1.0

    if (__switch__[] + rank) % 2 == 0
        comm_with = rank + 1

        if comm_with <= comm_size
            remotecall(left_is_ready, comm_with)
            take!(__right_ready__)

            remote_E_tot = take!(__E_tot__)
            remote_beta = take!(__beta__)

            @fastmath dEdT = (remote_E_tot - E_tot) * (remote_beta - beta)
            p = min(1.0, exp(dEdT))
            do_swap = (p == 1) || (p > rand())
            remotecall(set_do_swap!, comm_with, do_swap)

            if do_swap
                old_spins = deepcopy(spins)
                remotecall(set_spins!, comm_with, old_spins)
                remotecall(set_E_tot!, comm_with, E_tot)
                spins .= take!(__spins__)
                init_edges!(sgraph, spins)
                return remote_E_tot, p
            end
        end
        return E_tot, p
    else
        comm_with = rank - 1

        if comm_with > 0
            remotecall(right_is_ready, comm_with)
            take!(__left_ready__)

            remotecall(set_E_tot!, comm_with, E_tot)
            remotecall(set_beta!, comm_with, beta)

            do_swap = take!(__do_swap__)

            if do_swap
                old_spins = deepcopy(spins)
                remotecall(set_spins!, comm_with, old_spins)
                spins = take!(__spins__)
                remote_E_tot = take!(__E_tot__)
                init_edges!(sgraph, spins)
                return remote_E_tot, p
            end
        end
        return E_tot, p
    end
    nothing
end


# NOTE for reference, to be deleted later
#=
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
        # print("[T = $T] $(toq())\n")
    end

    1. / T
end
=#
