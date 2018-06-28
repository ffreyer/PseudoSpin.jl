"""
    init_edges!(sgraph, spins)

Initialize the pre-calculated values of nearest neighbor edges. This has to
happen before the first sweep (and any other operation requiring scalar products)
"""
function init_edges!(sgraph::SGraph, spins::Vector{Point3{Float64}})
    for e in sgraph.first
        @fastmath @inbounds e.xy = spins[e.n1][1] * spins[e.n2][1] + spins[e.n1][2] * spins[e.n2][2]
        @fastmath @inbounds e.z = spins[e.n1][3] * spins[e.n2][3]
    end
    nothing
end


################################################################################
#### total Energy functions
################################################################################


# Note:
# The order of calculations is not optimized here. This is fine because the
# function only gets called once and helps verifying both totalEnergy and
# deltaEnergy.
"""
    totalEnergy(sgraph, spins, Js, h)

Calculates the total energy of the current system.
"""
function totalEnergy(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Vector{Tuple{Float64, Float64}},
        h::Point3{Float64}=Point3(0.)
    )
    # println("totalEnergy w/o g...")  # NOTE
    # println("\t Js = ", Js)
    # println("\t h = ", h)


    E = 0.
    for e in sgraph.second
        E += Js[2][1] * (
            spins[e.n1][1] * spins[e.n2][1] +
            spins[e.n1][2] * spins[e.n2][2]
        ) + Js[2][2] * spins[e.n1][3] * spins[e.n2][3]
    end

    for i in eachindex(sgraph.first)
        e = sgraph.first[i]
        E += Js[1][1] * e.xy + Js[1][2] * e.z
        for p in sgraph.paths[i]
            # I'll leave this as it is, because it's much easier to see whether
            # this is right. With this being right one can easily check whether
            # dE is right because E_tot2 == E_tot1 + dE has to be true
            E += (Js[3][1] * e.xy + Js[3][2] * e.z) * (Js[4][1] * p.xy + Js[4][2] * p.z)
            E += (Js[4][1] * e.xy + Js[4][2] * e.z) * (Js[3][1] * p.xy + Js[3][2] * p.z)
        end
    end

    if !(Tuple(h) == (0., 0., 0.))
        for S in spins
            E -= dot(h, S)
        end
    end

    E
end


# total Energy for (Js, h, g)
function totalEnergy(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Vector{Tuple{Float64, Float64}},
        h::Point3{Float64},
        g::Float64
    )
    # println("totalEnergy w/ g...")  # NOTE
    # println("\t Js = ", Js)
    # println("\t h = ", h)
    # println("\t g = ", g)


    E = 0.
    for e in sgraph.second
        E += Js[2][1] * (
            spins[e.n1][1] * spins[e.n2][1] +
            spins[e.n1][2] * spins[e.n2][2]
        ) + Js[2][2] * spins[e.n1][3] * spins[e.n2][3]
    end

    for i in eachindex(sgraph.first)
        e = sgraph.first[i]
        E += Js[1][1] * e.xy + Js[1][2] * e.z
        for p in sgraph.paths[i]
            E += (Js[3][1] * e.xy + Js[3][2] * e.z) * (Js[4][1] * p.xy + Js[4][2] * p.z)
            E += (Js[4][1] * e.xy + Js[4][2] * e.z) * (Js[3][1] * p.xy + Js[3][2] * p.z)
        end
    end

    # Checked:
    # factor 0.5 (overcounting edges)
    # number of edges (12x a-b-c, 2x (overcounting) 6x c-a-b)
    # edges/paths for first node correct
    for ei in eachindex(sgraph.first)
        e12 = sgraph.first[ei]
        for e23 in sgraph.nodes[e12.n2].first   # LID: 1 -> 2 -> 1
            e12 == e23 && continue
            n3 = e23.n1 != e12.n2 ? e23.n1 : e23.n2
            E += 0.5g * (
                spins[e12.n1][1] * spins[e12.n2][1] * spins[n3][2] +
                spins[e12.n1][1] * spins[e12.n2][2] * spins[n3][1] +
                spins[e12.n1][2] * spins[e12.n2][1] * spins[n3][1] -
                spins[e12.n1][2] * spins[e12.n2][2] * spins[n3][2]
            )
        end
        for e23 in sgraph.nodes[e12.n1].first   # LID: 2 -> 1 -> 2
            e12 == e23 && continue              # overcounting cause 2 <- 1 <- 2
            n3 = e23.n1 != e12.n1 ? e23.n1 : e23.n2
            E += 0.5g * (
                spins[e12.n2][1] * spins[e12.n1][1] * spins[n3][2] +
                spins[e12.n2][1] * spins[e12.n1][2] * spins[n3][1] +
                spins[e12.n2][2] * spins[e12.n1][1] * spins[n3][1] -
                spins[e12.n2][2] * spins[e12.n1][2] * spins[n3][2]
            )
        end
    end

    if !(Tuple(h) == (0., 0., 0.))
        for S in spins
            E -= dot(h, S)
        end
    end

    E
end


################################################################################
### scalar product stuff
################################################################################


"""
    scalar_prod(edge, index, new_spin, spins)

Calculates new values to be allocated as scalar product on Edge e.
"""
function scalar_prod(
        e::SEdge1,
        i::Int64,
        new_spin::Point3{Float64},
        spins::Vector{Point3{Float64}}
    )

    j = i == e.n1 ? e.n2 : e.n1
    @fastmath @inbounds return new_spin[1] * spins[j][1] + new_spin[2] * spins[j][2], new_spin[3] * spins[j][3]
end


"""
    generate_scalar_products(sgraph, spins, i, new_spin)

Calculates all new scalar products for a node at index i with spin new_spin.
"""
function generate_scalar_products(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64}
    )

    xys = Array{Float64}(4)
    zs = Array{Float64}(4)

    for j in eachindex(sgraph.nodes[i].first)
        @inbounds xys[j], zs[j] = scalar_prod(
            sgraph.nodes[i].first[j], i, new_spin, spins
        )
    end

    xys, zs
end


"""
    update_edges!(node, xys, zs)

Pushes new scalar products (in lists xys and zs) to the relevant edges of node.
"""
function update_edges!(n::SNode, xys::Vector{Float64}, zs::Vector{Float64})
    for j in eachindex(n.first)
        @inbounds n.first[j].xy = xys[j]
        @inbounds n.first[j].z = zs[j]
    end
    nothing
end


################################################################################
#### new stack
################################################################################

# PLAN
# - generate methods such as sweepJ1J2g(..., J1, J2, K, g, h) specialized for
#   arguments J1, J2, g (silently assuming that K, h are irrelevant/0.0)
#   *
# - generate/write sweep_picker(; J1 = 0.0, ...) which return the correct sweep
#   * call this early so that later calls of sweep() are efficient
#       (see pass_function_perf.jl in random_things for performance analysis)
#   * update all sweep/spin_flip (kernel)/deltaEnergy calls to have some
#       normal ordering. Maybe this would be easier with a Container struct?
#       Efficiency is fine
struct Parameters
    J1::Tuple{Float64, Float64}
    J2::Tuple{Float64, Float64}
    K::Float64
    g::Float64
    h::Point3{Float64}
end
function Parameters(;
        J1::Float64 = 0.0,
        J2::Float64 = 0.0,
        lambda::Float64 = 1.0,
        J1s::Tuple{Float64, Float64} = (J1, lambda*J1),
        J2s::Tuple{Float64, Float64} = (J2, lambda*J2),
        K::Float64 = 0.0,
        g::Float64 = 0.0,
        h::Point3{Float64} = Point3(0.0)
    )
    Parameters(J1s, J2s, K, g, h)
end

# - The first entry defines "normal order".
# - Every entry should follow "normal order"
# - Every sweep/spin_flip/deltaEnergy function is named according to the order
#   given here
# - sweep_picker has to follow this convention
const param_groups = [
    [:J1, :J2, :K, :g, :h],
    [:J1, :J2, :K, :g],
    [:J1, :J2, :K, :h],
    [:J1, :J2, :K],
    [:J1, :g, :h],
    [:J1, :g],
    # ...
]

function sweep_picker(param::Parameters)
    doJ1 = param.J1 != (0.0, 0.0)
    doJ2 = param.J2 != (0.0, 0.0)
    doK = param.K != 0.0
    dog = param.g != 0.0
    doh = param.h != Point3(0.0)

    param_group = Symbol[]
    # Normal order!
    doJ1 && push!(param_group, :J1)
    doJ2 && push!(param_group, :J2)
    doK && push!(param_group, :K)
    dog && push!(param_group, :g)
    doh && push!(param_group, :h)

    if param_group in param_groups
        return eval(:($(Symbol(:sweep_, param_group...))))
    else
        warn(
            "No method generated for (" *
            mapreduce(string, (a, b) -> a * ", " * b, param_group) *
            "). Using the most general method instead. Consider implementing" *
            " a specialized method by adding the parameters to param_groups!"
        )
        param_group = [:J1, :J2, :K, :g, :h]
        return eval(:($(Symbol(:sweep_, param_group...))))
    end
end

# Code Generation
for param_group in param_groups
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #<<< generate sweep functions
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    @eval begin
        # Function names such as sweep_J1J2g
        function $(Symbol(:sweep_, param_group...))(
                sgraph::SGraph,
                spin::Vector{Point3{Float64}},
                E_tot::Float64,
                beta::Float64,
                param::Parameters
            )
            for (i, new_spin) in zip(
                    rand(1:sgraph.N_nodes, sgraph.N_nodes),
                    rand_spin(sgraph.N_nodes)
                )
                # calls spin_flip_J1J2g
                E_tot = $(Symbol(:spin_flip_, param_group...))(
                    sgraph, spins, i, new_spin, E_tot, beta, param
                )
            end

            E_tot
        end
    end

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #<<< generate spin_flip functions
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    @eval begin
        function $(Symbol(:spin_flip_, param_group...))(
                sgraph::SGraph,
                spins::Vector{Point3{Float64}},
                i::Int64,
                new_spin::Point3{Float64},
                E_tot::Float64,
                beta::Float64,
                param::Parameters
            )
            @inbounds n = sgraph.nodes[i]
            xys, zs = generate_scalar_products(sgraph, spins, i, new_spin)
            dE = $(Symbol(:deltaEnergy_, param_group...))(
                n, spins, i, new_spin, xys, zs, param
            )

            if dE < 0.
                @inbounds spins[i] = new_spin
                update_edges!(n, xys, zs)
                return E_tot + dE
            elseif rand() < exp(-dE * beta)
                @inbounds spins[i] = new_spin
                update_edges!(n, xys, zs)
                return E_tot + dE
            end

            E_tot
        end
    end

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #<<< generate deltaEnergy functions
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    doJ1 = :J1 in param_group
    doJ2 = :J2 in param_group
    doK = :K in param_group
    dog = :g in param_group
    doh = :h in param_group

    @eval begin
        function $(Symbol(:deltaEnergy_, param_group...))(
                n::SNode,
                spins::Vector{Point3{Float64}},
                i::Int64,
                new_spin::Point3{Float64},
                xys::Vector{Float64},
                zs::Vector{Float64},
                param::Parameters
            )

            dE = 0.

            # J1/NN accumulator variables
            $(doJ1 && quote #---------------------------------------------------
                xy = 0.
                z = 0.
            end) #--------------------------------------------------------------

            # K/4-spin accumulator variable
            $(doK && quote #----------------------------------------------------
                xyz = 0.
            end) #--------------------------------------------------------------

            # J1, K loop start
            $((doJ1 || doK) && quote #------------------------------------------
                for j in 1:4
                    @inbounds xy = xys[j]
                    @inbounds z = zs[j]

                    @inbounds e = n.first[j]
                    @fastmath dxy = xy - e.xy
                    @fastmath dz = z - e.z

                    # J1 only part - accumulate xy, z
                    $(doJ1 && quote
                        @fastmath xy += dxy        # xy - e.xy
                        @fastmath z += dz          # z - e.z
                    end)

                    # K only
                    $(doK && quote
                        @fastmath temp_xy = 0.
                        @fastmath temp_z = 0.
                        for p in n.paths[j]
                            @fastmath temp_xy += p.xy
                            @fastmath temp_z += p.z
                        end
                        @fastmath xyz += dxy * temp_z + dz * temp_xy  # (K-term)
                        # @fastmath K_xy += dxy * temp_xy       # (general 4-spin)
                        # @fastmath K_z += dz * temp_z          # (general 4-spin)
                    end)
                end
            end) #--------------------------------------------------------------

            # Add to dE
            $(doJ1 && quote #---------------------------------------------------
                @fastmath @inbounds dE += param.J1[1] * xy + param.J1[2] * z
            end) #--------------------------------------------------------------
            $(doK && quote
                @fastmath @inbounds dE += 2 * param.K * xyz
            end) #--------------------------------------------------------------
                    # Js[1][1] * xy1 +                    # J1[1]
                    # Js[1][2] * z1 +                     # J1[2]
                    # 2 * (                               # K here
                    #     Js[3][1] * Js[4][1] * xy3 +     # K * 0     (x)
                    #     Js[3][2] * Js[4][2] * z3        # 0 * 1     (x)
                    # ) + (
                    #     Js[3][1] * Js[4][2] +           # k * 1     (keep)
                    #     Js[3][2] * Js[4][1]             # 0 * 1     (x)
                    # ) * xyz3



            # for g, h, J2
            $((doJ2 || dog || doh) && quote #-----------------------------------
                @fastmath @inbounds delta_s = new_spin .- spins[i]
            end) #--------------------------------------------------------------


            # J2/NNN
            $(doJ2 && quote #---------------------------------------------------
                # xy and z can be overwritten
                xy = 0.
                z = 0.
                for j in n.second
                    @fastmath @inbounds xy += delta_s[1] * spins[j][1] +
                                              delta_s[2] * spins[j][2]
                    @fastmath @inbounds z += delta_s[3] * spins[j][3]
                end
                @fastmath @inbounds dE += param.J2[1] * xy + param.J2[2] * z
            end) #--------------------------------------------------------------

            #=
            g-paths
            original: ∑∑∑ xxy + xyx + yxx - yyy
            affected spins: o-a-b, a-o-b, where o is flipped
            Notation: dx is flipped, sum over paths (o-a-b, a-o-b) implied
            o-a-b: dx(xy)   + dx(yx)   + dy(xx)   - dy(yy)
            a-o-b: (x)dx(y) + (x)dy(x) + (y)dx(x) - y(dy)y
            notation: 1 means o-a-b term, 2 means a-o-b term
            dx * [(xy)1 + (yx)1 + (xy)2 + (yx)2]
            dy * [(xx)1 - (yy)1 + (xx)2 - (yy)2]
            Note symmetry - we can loop over (o-a-b) and (a-o-b) terms here
            further bracketing:
            dx * [x1 * (y1 + y2) + y1 * (x1 + x2)]
            dy * [x1 * (x1 + x2) - y1 * (y1 + y2)]
            Note dublicate terms - we only need to compute two inner brackets
            =#
            # g terms
            $(dog && quote #----------------------------------------------------
                g_x__ = 0.  # [x1 * g___y + y1 * g___x]
                g_y__ = 0.  # [x1 * g___x - y1 * g___y]
                for gedge in n.gpaths   # 4 elements
                    g___x = 0.  # (x1 + x2)
                    g___y = 0.  # (y1 + y2)
                    for bi in gedge.bs  # 6 5 4 3 elements
                        @inbounds s2 = spins[bi]
                        @fastmath @inbounds g___x += s2[1]
                        @fastmath @inbounds g___y += s2[2]
                    end
                    @inbounds s1 = spins[gedge.a]
                    @fastmath @inbounds g_x__ += s1[1] * g___y + s1[2] * g___x
                    @fastmath @inbounds g_y__ += s1[1] * g___x - s1[2] * g___y
                end
                @fastmath @inbounds dE += param.g *
                    (delta_s[1] * g_x__ + delta_s[2] * g_y__)
            end) #--------------------------------------------------------------

            # h
            $(doh && quote #----------------------------------------------------
                @fastmath dE -= dot(param.h, delta_s)
            end) #--------------------------------------------------------------
        end
    end
end


# For reference: g-paths
# TODO: find paths in sgraph; remove searching here
# @inbounds sj = n.first[j].n1 != i ? n.first[j].n1 : n.first[j].n2
# @fastmath @inbounds dxxyy = delta_s[1] * spins[sj][1] -
#                             delta_s[2] * spins[sj][2]
# @fastmath @inbounds dxyyx = delta_s[1] * spins[sj][2] +
#                             delta_s[2] * spins[sj][1]
#
# x = 0.
# y = 0.
#
# # second part (x - a - b)
# for k in 1:4
#     @inbounds _n = sgraph.nodes[sj]
#     @inbounds n.first[j] == _n.first[k] && continue
#     @inbounds sk = if _n.first[k].n1 != sj
#         _n.first[k].n1
#     else
#         _n.first[k].n2
#     end
#     @fastmath @inbounds x += spins[sk][1]
#     @fastmath @inbounds y += spins[sk][2]
# end
#
# # second part (b - x - a)
# for k in j+1:4
#     @inbounds sk = if n.first[k].n1 != i
#         n.first[k].n1
#     else
#         n.first[k].n2
#     end
#     @fastmath @inbounds x += spins[sk][1]
#     @fastmath @inbounds y += spins[sk][2]
# end
#
# @fastmath g_temp += dxxyy * y + dxyyx * x

################################################################################
#### single spin flip kernels
################################################################################




#=
################################################################################
#### delta Energy Functions
################################################################################
# Note:
# I tried to optimize this stuff as much as I could...


# Energy shift for anisotropic Js
# Checked vs totalEnergy() for anisotropic J3, J4 (aka K) in reordered 144-term
# form. Error ~1e-12% (1e-14 as factor)
"""
    deltaEnergy(node, spins, index, new_spin, Js, xys, zs, h)

Calculates the energy difference for a spin at some index changing to new_spin.
"""
function deltaEnergy(
        n::SNode,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        Js::Vector{Tuple{Float64, Float64}},
        xys::Vector{Float64},
        zs::Vector{Float64},
        h::Point3{Float64}=Point3(0.)
    )

    # allocations for nearest neighbour
    xy1 = 0.
    z1 = 0.

    # allocations for paths
    xyz3 = 0.
    xy3 = 0.
    z3 = 0.

    # Calculate nearest neighbours and paths
    for (xy, z, j) in zip(xys, zs, eachindex(n.first))
        @inbounds e = n.first[j]
        @fastmath dxy = xy - e.xy
        @fastmath dz = z - e.z
        @fastmath xy1 += dxy        # xy - e.xy
        @fastmath z1 += dz          # z - e.z
        @fastmath temp_xy = 0.
        @fastmath temp_z = 0.

        for p in n.paths[j]
            @fastmath temp_z += p.z
            @fastmath temp_xy += p.xy
        end

        @fastmath xyz3 += dxy * temp_z + dz * temp_xy
        @fastmath xy3 += dxy * temp_xy
        @fastmath z3 += dz * temp_z
    end
    @fastmath @inbounds dE = begin
        Js[1][1] * xy1 +
        Js[1][2] * z1 +
        2 * (
            Js[3][1] * Js[4][1] * xy3 +
            Js[3][2] * Js[4][2] * z3
        ) + (
            Js[3][1] * Js[4][2] +
            Js[3][2] * Js[4][1]
        ) * xyz3
    end

    # Calculate Next Nearest Neighbor terms
    @inbounds delta_s = new_spin - spins[i]
    xy1 = 0.
    z1 = 0.
    for j in n.second
        @fastmath @inbounds xy1 += delta_s[1] * spins[j][1] + delta_s[2] * spins[j][2]
        @fastmath @inbounds z1 += delta_s[3] * spins[j][3]
    end
    @fastmath @inbounds dE += Js[2][1] * xy1 + Js[2][2] * z1

    # Calculate field term
    if !(Tuple(h) == (0., 0., 0.))
        @fastmath @inbounds dE -= dot(h, (new_spin - spins[i]))
    end

    dE
end


# Energy shift for (Js, h, g)
function deltaEnergy(
        #n::SNode,
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        Js::Vector{Tuple{Float64, Float64}},
        xys::Vector{Float64},
        zs::Vector{Float64},
        h::Point3{Float64},
        g::Float64
    )

    @inbounds n = sgraph.nodes[i]
    @inbounds delta_s = new_spin - spins[i]

    # allocations for nearest neighbour
    xy1 = 0.
    z1 = 0.

    # allocations for paths
    xyz3 = 0.
    xy3 = 0.
    z3 = 0.

    # Allocation for 3-spin
    temp = 0.

    # Calculate nearest neighbours and paths
    for (xy, z, j) in zip(xys, zs, eachindex(n.first))
        @inbounds e = n.first[j]    # NN edge
        @fastmath dxy = xy - e.xy   # ΔS_xy
        @fastmath dz = z - e.z      # ΔS_z
        @fastmath xy1 += dxy        # xy - e.xy
        @fastmath z1 += dz          # z - e.z
        @fastmath temp_xy = 0.
        @fastmath temp_z = 0.

        for p in n.paths[j]
            @fastmath temp_xy += p.xy
            @fastmath temp_z += p.z
        end

        @fastmath xyz3 += dxy * temp_z + dz * temp_xy
        @fastmath xy3 += dxy * temp_xy
        @fastmath z3 += dz * temp_z

        # g/3-spin stuff
        # sj = n.first[j].n1 != i ? n.first[j].n1 : n.first[j].n2
        # dxxyy = delta_s[1] * spins[sj][1] - delta_s[2] * spins[sj][2]
        # dxyyx = delta_s[1] * spins[sj][2] + delta_s[2] * spins[sj][1]
        #
        # x = 0.
        # y = 0.
        #
        # # second part (x - a - b)
        # for k in 1:4
        #     _n = sgraph.nodes[sj]
        #     n.first[j] == _n.first[k] && continue
        #     sk = _n.first[k].n1 != sj ? _n.first[k].n1 : _n.first[k].n2
        #     x += spins[sk][1]
        #     y += spins[sk][2]
        # end
        #
        # # second part (b - x - a)
        # for k in j+1:4
        #     sk = n.first[k].n1 != i ? n.first[k].n1 : n.first[k].n2
        #     x += spins[sk][1]
        #     y += spins[sk][2]
        # end
        #
        # temp += dxxyy * y + dxyyx * x
    end
    @fastmath @inbounds dE = begin
        Js[1][1] * xy1 +
        Js[1][2] * z1 +
        2 * (
            Js[3][1] * Js[4][1] * xy3 +
            Js[3][2] * Js[4][2] * z3
        ) + (
            Js[3][1] * Js[4][2] +
            Js[3][2] * Js[4][1]
        ) * xyz3 #+
        # g * temp
    end


    g_x__ = 0.  # [x1 * g___y + y1 * g___x]
    g_y__ = 0.  # [x1 * g___x - y1 * g___y]
    for gedge in n.gpaths   # 4
        g___x = 0.  # (x1 + x2)
        g___y = 0.  # (y1 + y2)
        for bi in gedge.bs  # 6 5 4 3
            @inbounds s2 = spins[bi]
            @fastmath @inbounds g___x += s2[1]
            @fastmath @inbounds g___y += s2[2]
        end
        @inbounds s1 = spins[gedge.a]
        @fastmath @inbounds g_x__ += s1[1] * g___y + s1[2] * g___x
        @fastmath @inbounds g_y__ += s1[1] * g___x - s1[2] * g___y
    end
    @fastmath @inbounds dE += g * (delta_s[1] * g_x__ + delta_s[2] * g_y__)


    # Calculate Next Nearest Neighbor terms
    xy1 = 0.
    z1 = 0.
    for j in n.second
        @fastmath @inbounds xy1 += delta_s[1] * spins[j][1] + delta_s[2] * spins[j][2]
        @fastmath @inbounds z1 += delta_s[3] * spins[j][3]
    end
    @fastmath @inbounds dE += Js[2][1] * xy1 + Js[2][2] * z1

    # Calculate field term
    if !(Tuple(h) == (0., 0., 0.))
        @fastmath @inbounds dE -= dot(h, (new_spin - spins[i]))
    end

    dE
end


# Energy shift for (Js, h, g)
function deltaEnergy_no_paths(
        #n::SNode,
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        Js::Vector{Tuple{Float64, Float64}},
        xys::Vector{Float64},
        zs::Vector{Float64},
        h::Point3{Float64},
        g::Float64
    )

    @inbounds n = sgraph.nodes[i]
    @inbounds delta_s = new_spin - spins[i]

    # allocations for nearest neighbour
    xy1 = 0.
    z1 = 0.

    # Allocation for 3-spin
    temp = 0.

    # Calculate nearest neighbours and paths
    for (xy, z, j) in zip(xys, zs, eachindex(n.first))
        @inbounds e = n.first[j]    # NN edge
        @fastmath dxy = xy - e.xy   # ΔS_xy
        @fastmath dz = z - e.z      # ΔS_z
        @fastmath xy1 += dxy        # xy - e.xy
        @fastmath z1 += dz          # z - e.z

        # g/3-spin stuff
        sj = n.first[j].n1 != i ? n.first[j].n1 : n.first[j].n2
        dxxyy = delta_s[1] * spins[sj][1] - delta_s[2] * spins[sj][2]
        dxyyx = delta_s[1] * spins[sj][2] + delta_s[2] * spins[sj][1]

        x = 0.
        y = 0.

        # second part (x - a - b)
        for k in 1:4
            _n = sgraph.nodes[sj]
            n.first[j] == _n.first[k] && continue
            sk = _n.first[k].n1 != sj ? _n.first[k].n1 : _n.first[k].n2
            x += spins[sk][1]
            y += spins[sk][2]
        end

        # second part (b - x - a)
        for k in j+1:4
            sk = n.first[k].n1 != i ? n.first[k].n1 : n.first[k].n2
            x += spins[sk][1]
            y += spins[sk][2]
        end

        temp += dxxyy * y + dxyyx * x
    end
    @fastmath @inbounds dE = begin
        Js[1][1] * xy1 +
        Js[1][2] * z1 +
        g * temp
    end


    # Calculate Next Nearest Neighbor terms
    xy1 = 0.
    z1 = 0.
    for j in n.second
        @fastmath @inbounds xy1 += delta_s[1] * spins[j][1] + delta_s[2] * spins[j][2]
        @fastmath @inbounds z1 += delta_s[3] * spins[j][3]
    end
    @fastmath @inbounds dE += Js[2][1] * xy1 + Js[2][2] * z1

    # Calculate field term
    if !(Tuple(h) == (0., 0., 0.))
        @fastmath @inbounds dE -= dot(h, (new_spin - spins[i]))
    end

    dE
end



################################################################################
#### single spin flip kernels
################################################################################


# This attempts a single spin flip for a given new_spin and array position i.
# Anisotropic version, no update to total Energy
"""
    kernel(sgraph, spins, index, new_spin[, E_tot], Js, beta, h)

Attempts a single spin flip. The current total energy E_tot can be given
optionally. If done so, the function will return the updated total energy.
"""
function kernel(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64}=Point3(0.)
    )

    # energy difference implemented for anisotropy
    # calculate new scalar products
    @inbounds n = sgraph.nodes[i]
    xys, zs = generate_scalar_products(sgraph, spins, i, new_spin)
    dE = deltaEnergy(n, spins, i, new_spin, Js, xys, zs, h)

    if dE < 0.
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return nothing
    elseif rand() < exp(-dE * beta)
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return nothing
    end

    nothing
end


# anisotropic Js, returns updated E_tot
function kernel(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64}=Point3(0.)
    )

    # energy difference implemented for anisotropy
    # calculate new scalar products
    @inbounds n = sgraph.nodes[i]
    xys, zs = generate_scalar_products(sgraph, spins, i, new_spin)
    dE = deltaEnergy(n, spins, i, new_spin, Js, xys, zs, h)

    if dE < 0.
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return E_tot + dE
    elseif rand() < exp(-dE * beta)
        # println(E_tot, " + ", dE, " = ", E_tot + dE, " \t ", dE*beta)
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return E_tot + dE
    end

    E_tot
end


# anisotropic Js, returns updated E_tot
function kernel(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        h::Point3{Float64}=Point3(0.)
    )

    # energy difference implemented for anisotropy
    # calculate new scalar products
    @inbounds n = sgraph.nodes[i]
    xys, zs = generate_scalar_products(sgraph, spins, i, new_spin)
    dE = deltaEnergy(n, spins, i, new_spin, Js, xys, zs, h)

    if dE < 0.
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return E_tot + dE
    end

    E_tot
end


# (Js, h, g)
function kernel(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64},
        g::Float64
    )

    @inbounds n = sgraph.nodes[i]
    xys, zs = generate_scalar_products(sgraph, spins, i, new_spin)
    dE = deltaEnergy(sgraph, spins, i, new_spin, Js, xys, zs, h, g)

    if dE < 0.
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return E_tot + dE
    elseif rand() < exp(-dE * beta)
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return E_tot + dE
    end

    E_tot
end


# (Js, h, g) w/o E updates
function kernel(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64},
        g::Float64
    )

    @inbounds n = sgraph.nodes[i]
    xys, zs = generate_scalar_products(sgraph, spins, i, new_spin)
    dE = deltaEnergy(sgraph, spins, i, new_spin, Js, xys, zs, h, g)

    if dE < 0.
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return nothing
    elseif rand() < exp(-dE * beta)
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return nothing
    end

    nothing
end


# (reduced Js, h, g)
function kernel_no_paths(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64},
        g::Float64
    )

    @inbounds n = sgraph.nodes[i]
    xys, zs = generate_scalar_products(sgraph, spins, i, new_spin)
    dE = deltaEnergy_no_paths(sgraph, spins, i, new_spin, Js, xys, zs, h, g)

    if dE < 0.
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return E_tot + dE
    elseif rand() < exp(-dE * beta)
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return E_tot + dE
    end

    E_tot
end


# (reduced Js, h, g) w/o E updates
function kernel_no_paths(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64},
        g::Float64
    )

    @inbounds n = sgraph.nodes[i]
    xys, zs = generate_scalar_products(sgraph, spins, i, new_spin)
    dE = deltaEnergy_no_paths(sgraph, spins, i, new_spin, Js, xys, zs, h, g)

    if dE < 0.
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return nothing
    elseif rand() < exp(-dE * beta)
        @inbounds spins[i] = new_spin
        update_edges!(n, xys, zs)
        return nothing
    end

    nothing
end


################################################################################
#### sweep functions
################################################################################


"""
    sweep(sgraph, spins[, E_tot], Js, beta[, h[, g]])

Attempts as many spin flips as there are sites in the lattice. The current total
energy E_tot can be given optionally. If done so, it will be updated after each
successful spin flip and returned in the end.
"""
function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64},
        g::Float64
    )
    for (i, new_spin) in zip(
            rand(1:sgraph.N_nodes, sgraph.N_nodes),
            rand_spin(sgraph.N_nodes)
        )
        E_tot = kernel(sgraph, spins, i, new_spin, E_tot, Js, beta, h, g)
    end

    E_tot
end


function sweep_no_paths(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64},
        g::Float64
    )
    for (i, new_spin) in zip(
            rand(1:sgraph.N_nodes, sgraph.N_nodes),
            rand_spin(sgraph.N_nodes)
        )
        E_tot = kernel_no_paths(sgraph, spins, i, new_spin, E_tot, Js, beta, h, g)
    end

    E_tot
end


function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64},
        g::Float64
    )
    for (i, new_spin) in zip(
            rand(1:sgraph.N_nodes, sgraph.N_nodes),
            rand_spin(sgraph.N_nodes)
        )
        kernel(sgraph, spins, i, new_spin, Js, beta, h, g)
    end

    nothing
end


function sweep_no_paths(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64},
        g::Float64
    )
    for (i, new_spin) in zip(
            rand(1:sgraph.N_nodes, sgraph.N_nodes),
            rand_spin(sgraph.N_nodes)
        )
        kernel_no_paths(sgraph, spins, i, new_spin, Js, beta, h, g)
    end

    nothing
end


# no beta
function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Vector{Tuple{Float64, Float64}},
        h::Point3{Float64}=Point3(0.)
    )
    for (i, new_spin) in zip(
        rand(1:sgraph.N_nodes, sgraph.N_nodes),
        rand_spin(sgraph.N_nodes)
        )
        kernel(sgraph, spins, i, new_spin, Js, h)
    end

    nothing
end


function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        h::Point3{Float64}=Point3(0.)
    )
    for (i, new_spin) in zip(
            rand(1:sgraph.N_nodes, sgraph.N_nodes),
            rand_spin(sgraph.N_nodes)
        )
        E_tot = kernel(sgraph, spins, i, new_spin, E_tot, Js, h)
    end

    E_tot
end



function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64}
    )
    for (i, new_spin) in zip(
            rand(1:sgraph.N_nodes, sgraph.N_nodes),
            rand_spin(sgraph.N_nodes)
        )
        kernel(sgraph, spins, i, new_spin, Js, beta, h)
    end

    nothing
end

function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64,
        h::Point3{Float64}
    )
    for (i, new_spin) in zip(
            rand(1:sgraph.N_nodes, sgraph.N_nodes),
            rand_spin(sgraph.N_nodes)
        )
        E_tot = kernel(sgraph, spins, i, new_spin, E_tot, Js, beta, h)
    end

    E_tot
end

# # (Js, h, g) w/o E update
# function sweep(
#         sgraph::SGraph,
#         spins::Vector{Point3{Float64}},
#         Js::Vector{Tuple{Float64, Float64}},
#         beta::Float64,
#         h::Point3{Float64},
#         g::Float64
#     )
#
#     for (i, new_spin) in zip(
#             rand(1:sgraph.N_nodes, sgraph.N_nodes),
#             rand_spin(sgraph.N_nodes)
#         )
#         kernel(sgraph, spins, i, new_spin, Js, beta, h, g)
#     end
#
#     nothing
# end
#
#
# # (UJs, h, g) w/ E update
# function sweep(
#         sgraph::SGraph,
#         spins::Vector{Point3{Float64}},
#         E_tot::Float64,
#         Js::Vector{Tuple{Float64, Float64}},
#         beta::Float64,
#         h::Point3{Float64},
#         g::Float64
#     )
#
#     for (i, new_spin) in zip(
#             rand(1:sgraph.N_nodes, sgraph.N_nodes),
#             rand_spin(sgraph.N_nodes)
#         )
#         E_tot = kernel(sgraph, spins, i, new_spin, E_tot, Js, beta, h, g)
#     end
#
#     E_tot
# end
=#
