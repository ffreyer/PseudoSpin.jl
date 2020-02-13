abstract type AbstractUpdate end


################################################################################
### Local updates
################################################################################

abstract type AbstractLocalUpdate <: AbstractUpdate end

function apply(g::AbstractLocalUpdate, ::Vector{SVector{3, Float64}})
    throw(ErrorException("Missing apply($(typeof(g)))"))
end

function accept(g::AbstractLocalUpdate) end


"""
    rand_spin([, N])

Returns a randomly oriented spin. (A point randomly picked from a unit sphere.)
If N is given, returns a *Vector* with N random spins. This is generally more
efficient in practise.
"""
function rand_spin()
    phi = 2 * pi * rand(Float64)
    ct = 2 * rand(Float64) - 1
    st = sqrt(1. - ct*ct) #sin(acos(ct)) # max(0., )

    SVector{3, Float64}(
        st * cos(phi),
        st * sin(phi),
        ct
    )
end

function rand_spin(N::Int64)
    phis = 2 * pi * rand(Float64, N)
    cts = 2 * rand(Float64, N) .- 1
    sts = sqrt.(1. .- cts .* cts) # sin(acos(cts)) # max(0., )

    [SVector{3, Float64}(
        sts[i] .* cos.(phis[i]),
        sts[i] .* sin.(phis[i]),
        cts[i]
    ) for i in 1:N]
end


"""
    rand_XY_spin([, N])

Returns a (an array of N) random spin(s) in the XY plane.
"""
function rand_XY_spin()
    phi = 2 * pi * rand(Float64)
    SVector{3, Float64}(cos(phi), sin(phi), 0.0)
end

function rand_XY_spin(N::Int64)
    phis = 2 * pi * rand(Float64, N)
    [SVector{3, Float64}(cos.(phis[i]), sin.(phis[i]), 0.0) for i in 1:N]
end


"""
    rand_XYm_spin([, N])

Returns a (an array of N) random spin(s) in a reduced XY plane with phi ∈
[3, 6.45] ~ [172°, 370°] (with 128 bins roughly 3 extra bins).
"""
function rand_red_XY_spin()
    phi = 3.0 + 3.45rand(Float64)
    SVector{3, Float64}(cos(phi), sin(phi), 0.0)
end

function rand_red_XY_spin(N::Int64)
    phis = 3.0 .+ 3.45 .* rand(Float64, N)
    [SVector{3, Float64}(cos.(phis[i]), sin.(phis[i]), 0.0) for i in 1:N]
end


# """
#     rand_XY_rot_matrix()
#
# Returns a random XY rotation matrix.
# """
# @inline function rand_XY_rot_matrix()
#     phi = 2pi * rand()
#     c = cos(phi)
#     s = sin(phi)
#
#     R = @SMatrix [
#         c  -s  0;
#         s   c  0;
#         0   0  1;
#     ]
# end


mutable struct self_balancing_update <: AbstractLocalUpdate
    M::SVector{3, Float64}
    eM::SVector{3, Float64}
    eM_perp::SVector{3, Float64}
    N::Int64
end

"""
    self_balancing_update(spins)

Initializes a `self_balancing_update` which keeps the diretion of magnetization
constant.

This update is biased. On top of a slight preference of spin pointing in the
direction of magnetization, there is also a strong preference of S = ±e_M⟂. 
"""
function self_balancing_update(spins)
    # constants
    M = sum(spins)
    eM = normalize(M)
    eM_perp = cross(SVector(0., 0., 1.), eM)
    N = length(spins)
    self_balancing_update(M, eM, eM_perp, N)
end


function apply(U::self_balancing_update, spins::Vector{SVector{3, Float64}})


    # algorithm
    new_spins = eltype(spins)[]
    idxs = Int64[]
    counter = 1

    while true
        # SSF
        new_S = rand_XY_spin()
        idxs = [trunc(Int64, 1 + U.N * rand())]

        new_spins = [new_S]
        dS = new_S - spins[idxs[1]]
        # component of dS in e_M⟂ direction
        x = dot(dS, U.eM_perp)

        # Somehow always converges
        @inbounds @fastmath while true
            i = trunc(Int64, 1 + U.N * rand())
            while i in idxs
                i = trunc(Int64, 1 + U.N * rand())
            end
            push!(idxs, i)

            # Maximum compensation
            # -sign(x) * 1.0 - sign(x) dot(spins[i], U.eM_perp)
            # ^- flip to ±e_M⟂ direction
            # previous component in e_M⟂ -^
            # sign(x) required to differentiate compensation & enhancement
            x -= dot(spins[i], U.eM_perp)
            if abs(x) > 1.0
                # Rotating to ±e_M⟂ not enough
                # do it, try again with another
                push!(new_spins, - sign(x) * U.eM_perp)
                x -= sign(x)
            else
                # can be compensated fully, calculate compensating rotation
                phi = acos(x)
                s = sin(phi)
                R = @SMatrix [
                    x -s  0.;
                    s  x  0.;
                    0. 0. 1.
                ]
                push!(new_spins, -R * U.eM_perp)
                break
            end
        end

        # if dot(sum(spins) - sum(spins[idxs]) + sum(new_spins), U.eM) < 0.0
        #     _new_spins = - copy(spins)
        #     _new_spins[idxs] = -new_spins
        #     return collect(eachindex(_new_spins)), _new_spins
        # end

        # Check if sign of e_M changed
        if dot(sum(spins) - sum(spins[idxs]) + sum(new_spins), U.eM) > 0.0
            break
        end

        # I feel like this is in the wrong place, but is never reached anyway
        if counter == U.N
            @warn "Proposing flip of all spins"
            _new_spins = - copy(spins)
            _new_spins[idxs] = -new_spins
            return collect(eachindex(_new_spins)), _new_spins
        end

        # Retry update if it failed
        counter += 1
    end

    idxs, new_spins
end


################################################################################
### Global updates
################################################################################


abstract type AbstractGlobalUpdate <: AbstractUpdate end


function apply(g::AbstractGlobalUpdate, ::Vector{SVector{3, Float64}})
    throw(ErrorException("Missing apply($(typeof(g)))"))
end

function accept(g::AbstractGlobalUpdate) end



"""
    rand_3fold_XY_rotation

Returns a random 3-fold rotation matrix, excluding identity rotations. So the
result will either be a rotation by +2pi/3 or -2pi/3 around the z/axis.
"""
struct rand_3fold_XY_rotation <: AbstractGlobalUpdate
    R1::SArray{Tuple{3,3},Float64,2,9}
    R2::SArray{Tuple{3,3},Float64,2,9}
    function rand_3fold_XY_rotation()
        R1 = @SMatrix[
            -0.5                 0.8660254037844387     0.0;
            -0.8660254037844387 -0.5                    0.0;
            0.0                  0.0                    1.0
        ]
        R2 = @SMatrix [
            -0.5                -0.8660254037844387     0.0;
            0.8660254037844387  -0.5                    0.0;
            0.0                  0.0                    1.0
        ]
        new(R1, R2)
    end
end


@inline function apply(
        x::rand_3fold_XY_rotation,
        spins::Vector{SVector{3, Float64}}
    )
    rot = rand(Bool) ? x.R1 : x.R2
    out = similar(spins)
    @inbounds @fastmath for i in eachindex(spins)
        out[i] = rot * spins[i]
    end
    out
end



"""
    yaxis_mirror(spins)

Returns an array of spins mirror at the y-(yz-)axis.
"""
struct yaxis_mirror <: AbstractGlobalUpdate end

@inline function apply(x::yaxis_mirror, spins::Vector{SVector{3, Float64}})
    [SVector(-S[1], S[2], S[3]) for S in spins]
end



"""
    flipflop_rot(spins)

Alternatingly flips by +/-Δϕ
"""
mutable struct flipflop_rot <: AbstractGlobalUpdate
    idx::Int64
    rots::Tuple{
        SArray{Tuple{3,3},Float64,2,9},
        SArray{Tuple{3,3},Float64,2,9}
    }
end
function flipflop_rot(delta_phi)
    c = cos(delta_phi)
    s = sin(delta_phi)

    R = @SMatrix [
        c  -s  0;
        s   c  0;
        0   0  1;
    ]
    Rinv = @SMatrix [
        c  s  0;
       -s  c  0;
        0  0  1;
    ]

    flipflop_rot(1, (R, Rinv))
end

function apply(x::flipflop_rot, spins::Vector{SVector{3, Float64}})
    R = x.rots[x.idx]
    [R * S for S in spins]
end

function accept(x::flipflop_rot)
    x.idx = 3 - x.idx
    nothing
end
