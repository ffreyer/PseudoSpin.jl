#=
            AbstractThermalizationMethod
AbstractTemperatureGenerator        AbstractParallelTemperingAlgorithm
- ConstantT                         - ParallelTempering
- Freezer                           - ProbabilityEqualizer
                                    TODO
                                    - TimeEqualizer
                                    - EnergyRebalancer
.
=#

# TODO
#=
Every adaptive algorithm has to thermalize thoroughly first
- Preferably PT + simulated annealing
- lot's of sweeps (same order as no PT!?)
- explicitly checking/recording energy equillibration might be useful (added)

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


# Interface methods

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
initialize

"""
    next(<:AbstractThermalizationMethod, state)

Returns the next set of variables as well as the next state.
"""
next

"""
    current_index(<:AbstractThermalizationMethod, state)

Returns the current index of the iteration.
"""
current_index

"""
    done(<:AbstractThermalizationMethod, state)

Returns true if the thermalization method has reached its last output.
"""
done

"""
    last(<:AbstractThermalizationMethod, state)

Returns the last set of variables of a thermalization process.
"""
last

"""
    lenght(<:AbstractThermalizationMethod)

Returns the length of a thermalization process.
"""
length


batch_size(::AbstractThermalizationMethod) = -1
adaptive_sample_size(::AbstractThermalizationMethod) = -1

# Some methods, like ConstantT, don't have a "T_max"
T_max(::AbstractThermalizationMethod) = -1.0


################################################################################
### ConstantT Methods
################################################################################


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
set_beta!(::ConstantT, state::Tuple{Float64, Int64}, beta::Float64) = (beta, state[2])


################################################################################
### Freezer Methods
################################################################################


struct Freezer <: AbstractTemperatureGenerator
    N::Int64
    N_switch::Int64
    N_per_exp::Int64
    T_max::Float64

    sin_values::Vector{Float64}
    exp_values::Vector{Float64}
    exp_deltas::Vector{Float64}
end

"""
    Freezer()

Creates a temperature generator with simulated annealing.
"""
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
        range(log(exp_strength + 1.), stop=0.0, length=N_exp_points)
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
set_beta!(::Freezer, state::Tuple, beta::Float64) = (state[1:end-2]..., 1.0/beta, beta)


################################################################################
### NoParallelTempering
################################################################################


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


################################################################################
### ParallelTempering (Generation)
################################################################################

##############################################
#=
Stacking parallel tempering and simulated annealing should be
fine, so let's allow this.

All adaptive methods share code with the basic pt algorithm,
so they should just spawm a pt...
- Nope, sometimes I need probabilities, sometimes not
=#
##############################################

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


struct ProbabilityEqualizer{
        ATG <: AbstractTemperatureGenerator
    } <: AbstractParallelTemperingAlgorithm
    Tgen::ATG
    batch_size::Int64
    skip::Int64
    M::Int64
    adaptive_sample_size::Int64
end

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


# This also acts as a barrier to all
const __all_gather_p_beta__ = Channel{Tuple{Int64, Float64, Float64}}(256)
set_p_beta!(value) = put!(__all_gather_p_beta__, value)
function all_gather_p_beta(p, beta)
    rank = myid()
    for comm_with in start_id[]:stop_id[]
        remotecall(set_p_beta!, comm_with, (rank, p, beta))
    end
    probs = Vector{Float64}(stop_id[])
    Ts = Vector{Float64}(stop_id[])
    for _ in start_id[]:stop_id[]
        i, p, beta = take!(__all_gather_p_beta__)
        probs[i] = p
        Ts[i] = 1.0 / beta
    end
    probs, Ts
end


function next(th::ProbabilityEqualizer, state::Tuple, E_tot::Float64)
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
            p = N_prob == 0.0 ? 0.0 : cum_prob / N_prob
            probs, Ts = all_gather_p_beta(p, beta)
            # comm_rank == 1 && println("Probs - $probs \t \t Ts - $Ts")

            if (comm_rank != start_id[]) && (comm_rank != stop_id[])
                prob_sum = sum(probs)
                norm_prob = probs ./ prob_sum
                mean_prob = mean(norm_prob)

                dT21 = (Ts[corr_id[]] - Ts[corr_id[]-1])
                dT23 = (Ts[corr_id[]] - Ts[corr_id[]+1])
                down = (norm_prob[corr_id[]-1] - mean_prob) * dT21
                up = (norm_prob[corr_id[]] - mean_prob) * dT23
                # down = 0.01(log(norm_prob[rank]) - log(mean_prob)) *
                #     (temperatures[rank+1] - temperatures[rank])
                # up = 0.01(log(norm_prob[rank+1]) - log(mean_prob)) *
                #     (temperatures[rank+1] - temperatures[rank+2])

                # Let's guarantee that we don't jumble the temperatures
                # maximum change locked to ±0.4ΔT
                new_T = min(
                    Ts[corr_id[]] - 0.4dT23, # dT23 < 0.0
                    max(
                        Ts[corr_id[]-1] + 0.6dT21,
                        Ts[corr_id[]] + (1.0 - k/th.M) * (down + up)
                    )
                )
                # println("$(Ts[comm_rank]) -> $new_T")

                beta = 1.0 / new_T
                Tgen_state = set_beta!(th.Tgen, Tgen_state, beta)
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
### parallel tempering implementation
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
# const __E_tot__ = Channel{Float64}(1)
# set_E_tot!(value) = put!(__E_tot__, value)
# const __beta__ = Channel{Float64}(1)
# set_beta!(value) = put!(__beta__, value)
# const __do_swap__ = Channel{Bool}(1)
# set_do_swap!(value) = put!(__do_swap__, value)
# const __spins__ = Channel{Vector{Point3{Float64}}}(1)
# set_spins!(value) = put!(__spins__, value)

# transfer (beta, E_tot, (rand), spins)
const __left_channel__ = Channel{Tuple{
    Float64, Float64, Float64, Vector{Point3{Float64}}
}}(1)
put_left!(value) = put!(__left_channel__, value)
const __right_channel__ = Channel{Tuple{
    Float64, Float64, Vector{Point3{Float64}}
}}(1)
put_right!(value) = put!(__right_channel__, value)

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

    if (__switch__[] + rank) % 2 == 0
        # me -> comm_with
        comm_with = rank + 1

        if comm_with <= stop_id[]
            p = rand()
            remotecall(put_left!, comm_with, (beta, E_tot, p, spins))
            remote_beta, remote_E_tot, remote_spins = take!(__right_channel__)

            @fastmath dEdT = (remote_E_tot - E_tot) * (remote_beta - beta)
            if (dEdT > 0.0) || (exp(dEdT) > p)
                spins .= remote_spins
                init_edges!(sgraph, spins)
                return remote_E_tot
            end
        end
    else
        # comm_with <- me
        comm_with = rank - 1

        if comm_with >= start_id[]
            remotecall(put_right!, comm_with, (beta, E_tot, spins))
            remote_beta, remote_E_tot, p, remote_spins = take!(__left_channel__)

            @fastmath dEdT = (remote_E_tot - E_tot) * (remote_beta - beta)
            if (dEdT > 0.0) || (exp(dEdT) > p)
                spins .= remote_spins
                init_edges!(sgraph, spins)
                return remote_E_tot
            end
        end
    end
    return E_tot
end


function _parallel_tempering_adaptive!(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        beta::Float64
    )
    error("BROKEN!")
    rank = myid()
    p = -1.0

    if (__switch__[] + rank) % 2 == 0
        comm_with = rank + 1

        if comm_with <= stop_id[]
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
    else
        comm_with = rank - 1

        if comm_with >= start_id[]
            remotecall(right_is_ready, comm_with)
            take!(__left_ready__)

            remotecall(set_E_tot!, comm_with, E_tot)
            remotecall(set_beta!, comm_with, beta)

            do_swap = take!(__do_swap__)

            if do_swap
                old_spins = deepcopy(spins)
                remotecall(set_spins!, comm_with, old_spins)
                spins .= take!(__spins__)
                remote_E_tot = take!(__E_tot__)
                init_edges!(sgraph, spins)
                return remote_E_tot, p
            end
        end
    end
    return E_tot, p
end
