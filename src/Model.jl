const sqrt3 = sqrt(3.0)

# The purpose of this type is to simplify branching.
# Without this, you may want to call sweep(... J1, g), but there is no efficient
# way to make sure that J1 propagates as J1 and g as g (instead of, say, K).
# This type presents a way to auto-complete parameter information.
# The branching is then done based on whether a field is 0.0 in sweep_picker()
struct Parameters
    J1::Tuple{Float64, Float64}
    J2::Tuple{Float64, Float64}
    J3::Tuple{Float64, Float64}
    K::Float64
    g::Float64
    h::Point3{Float64}
    zeta::Float64
end

"""
    Paramaters(kwargs...)

Creates a Parameters type for a given set of parameters. Supported are J1, J2,
lambda, K, g, h as well as J1s and J2s, while are built from J1/J2 and lambda.
"""
function Parameters(;
        J1::Float64 = 0.0,
        J2::Float64 = 0.0,
        J3::Float64 = 0.0,
        lambda::Float64 = 1.0,
        J1s::Tuple{Float64, Float64} = (J1, lambda*J1),
        J2s::Tuple{Float64, Float64} = (J2, lambda*J2),
        J3s::Tuple{Float64, Float64} = (J3, lambda*J3),
        K::Float64 = 0.0,
        g::Float64 = 0.0,
        h::Point3{Float64} = Point3(0.0),
        zeta::Float64 = 0.0
    )
    Parameters(J1s, J2s, J3s, K, g, h, zeta)
end


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

# To avoid re-allocating these vectors every spin_flip, keep them as globals
const __xys__ = Vector{Float64}(undef, 4)
const __zs__  = Vector{Float64}(undef, 4)
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
    @fastmath @inbounds return new_spin[1] * spins[j][1] +
                               new_spin[2] * spins[j][2],
                               new_spin[3] * spins[j][3]
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

    for j in eachindex(sgraph.nodes[i].first)
        @inbounds __xys__[j], __zs__[j] = scalar_prod(
            sgraph.nodes[i].first[j], i, new_spin, spins
        )
    end

    __xys__, __zs__
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
#### Sweep stack
################################################################################


# This defines all the functions that are created for branching.
# Rules to make live simpler:
# - The first entry is the most general.
# - The first entry defines "normal order".
# - Every entry should follow "normal order"
# Every sweep/spin_flip/deltaEnergy function will be named according to the
# order given here. sweep_picker has to follow this convention.
const param_groups = [
    [:J1, :J2, :J3, :K, :g, :h, :zeta],
    [:J1, :J2, :K, :g, :h, :zeta],
    [:J1, :J2, :K, :g, :h], # Tested
    [:J1, :J2, :K, :g],
    [:J1, :J2, :K, :h],     # Tested
    [:J1, :J2, :K],
    [:J1, :g, :h],
    [:J1, :g],      # Tested
    [:J1],          # Tested
    [:J2],          # Tested
    [:J3],
    [:K],           # Tested
    [:h],           # Tested
    [:g],
    [:zeta],
    # ...
]

"""
    sweep_picker(parameters)

Returns a suitable sweep function for a given set of parameters. Will default to
the most general sweep function if no specialized code is available (with Warning).
"""
function sweep_picker(param::Parameters)
    doJ1 = param.J1 != (0.0, 0.0)
    doJ2 = param.J2 != (0.0, 0.0)
    doJ3 = param.J3 != (0.0, 0.0)
    doK = param.K != 0.0
    dog = param.g != 0.0
    doh = param.h != Point3(0.0)
    dozeta = param.zeta != 0.0

    param_group = Symbol[]
    # Normal order!
    doJ1 && push!(param_group, :J1)
    doJ2 && push!(param_group, :J2)
    doJ3 && push!(param_group, :J3)
    doK && push!(param_group, :K)
    dog && push!(param_group, :g)
    doh && push!(param_group, :h)
    dozeta && push!(param_group, :zeta)

    if param_group in param_groups
        return eval(:($(Symbol(:sweep_, param_group...))))
    else
        @warn(
            "No method generated for (" *
            reduce((a, b) -> a * ", " * b, map(string, param_group)) *
            "). Using the most general method instead. Consider implementing" *
            " a specialized method by adding the parameters to param_groups!"
        )
        param_group = param_groups[1]
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
                spins::Vector{Point3{Float64}},
                E_tot::Float64,
                beta::Float64,
                param::Parameters
            )
            for _ in 1:sgraph.N_nodes
                E_tot = $(Symbol(:spin_flip_, param_group...))(
                    sgraph, spins,
                    trunc(Int64, 1 + sgraph.N_nodes * rand()), rand_spin(),
                    E_tot, beta, param
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
    doJ3 = :J3 in param_group
    doK = :K in param_group
    dog = :g in param_group
    doh = :h in param_group
    dozeta = :zeta in param_group

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

            # J1, K loop
            $((doJ1 || doK) && quote #------------------------------------------
                for j in 1:4
                    @inbounds e = n.first[j]
                    @fastmath @inbounds dxy = xys[j] - e.xy
                    @fastmath @inbounds dz = zs[j] - e.z

                    # J1 only part - accumulate xy, z
                    $(doJ1 && quote
                        @fastmath xy += dxy        # xy - e.xy
                        @fastmath z += dz          # z - e.z
                    end)

                    # K only part
                    # H_K = K * ∑ijkl (Si * Sj) * (Sk * Sl)
                    # H_K = K * ∑ij (Si * Sj) * (∑kl(j) Sk * Sl)
                    #               (dxy, dz)   (temp_xy, temp_z)
                    $(doK && quote
                        @fastmath temp_xy = 0.
                        @fastmath temp_z = 0.
                        for p in n.paths[j]
                            @fastmath temp_xy += p.xy
                            @fastmath temp_z += p.z
                        end
                        @fastmath xyz += dxy * temp_z + dz * temp_xy  # (K-term)
                        # @fastmath K_xy += dxy * temp_xy     # (general 4-spin)
                        # @fastmath K_z += dz * temp_z        # (general 4-spin)
                    end)
                end
            end) #--------------------------------------------------------------

            # Add to dE - J1
            $(doJ1 && quote #---------------------------------------------------
                @fastmath @inbounds dE += param.J1[1] * xy + param.J1[2] * z
            end) #--------------------------------------------------------------

            # Add to dE - K
            $(doK && quote #----------------------------------------------------
                @fastmath @inbounds dE += param.K * xyz
            end) #--------------------------------------------------------------
            # General 4-spin term:
            # 2 * (                               # K here
            #     Js[3][1] * Js[4][1] * xy3 +     # K * 0     (x)
            #     Js[3][2] * Js[4][2] * z3        # 0 * 1     (x)
            # ) + (
            #     Js[3][1] * Js[4][2] +           # K * 1     (keep)
            #     Js[3][2] * Js[4][1]             # 0 * 1     (x)
            # ) * xyz3

            # ΔS for g, h, J2
            $((doJ2 || doJ3 || dog || doh || dozeta) && quote #-----------------
                @fastmath @inbounds delta_s = new_spin .- spins[i]
            end) #--------------------------------------------------------------

            # J2/NNN
            $(doJ2 && quote #---------------------------------------------------
                # xy and z can be overwritten, since dE(J1) has been added to dE
                xy = 0.
                z = 0.
                for j in n.second
                    # @inbounds will fail if the line is split normally ?
                    # (meaning no begin ... end block)
                    @fastmath @inbounds begin
                        xy += delta_s[1] * spins[j][1] +
                              delta_s[2] * spins[j][2]
                    end
                    @fastmath @inbounds z += delta_s[3] * spins[j][3]
                end
                @fastmath @inbounds dE += param.J2[1] * xy + param.J2[2] * z
            end) #--------------------------------------------------------------

            # J3 / third neighbor
            $(doJ3 && quote #---------------------------------------------------
                xy = 0.
                z = 0.
                for j in n.third
                    @fastmath @inbounds begin
                        xy += delta_s[1] * spins[j][1] +
                              delta_s[2] * spins[j][2]
                    end
                    @fastmath @inbounds z += delta_s[3] * spins[j][3]
                end
                @fastmath @inbounds dE += param.J3[1] * xy + param.J3[2] * z
            end) #--------------------------------------------------------------

            $(dozeta && quote #-------------------------------------------------
                id = n.second
                @inbounds dE += param.zeta * (
                    delta_s[1] * (
                        2.0 * ( # xy-plane terms
                            spins[id[1]][1] + spins[id[5]][1] +
                            spins[id[8]][1] + spins[id[12]][1]
                        ) - ( # yz, xz -plane terms
                            spins[id[2]][1] + spins[id[3]][1] +
                            spins[id[4]][1] + spins[id[6]][1] +
                            spins[id[7]][1] + spins[id[9]][1] +
                            spins[id[10]][1] + spins[id[11]][1]
                        ) + sqrt3 * ( # xz-plane
                            - spins[id[2]][2] - spins[id[4]][2] -
                            spins[id[9]][2] - spins[id[11]][2] +
                                     # yz-plane
                            spins[id[3]][2] + spins[id[6]][2] +
                            spins[id[7]][2] + spins[id[10]][2]
                        )
                    ) + delta_s[2] * (
                        -2.0 * ( # xy-plane terms
                            spins[id[1]][2] + spins[id[5]][2] +
                            spins[id[8]][2] + spins[id[12]][2]
                        ) + ( # yz, xz -plane terms
                            spins[id[2]][2] + spins[id[3]][2] +
                            spins[id[4]][2] + spins[id[6]][2] +
                            spins[id[7]][2] + spins[id[9]][2] +
                            spins[id[10]][2] + spins[id[11]][2]
                        ) + sqrt3 * ( # xz-plane
                            - spins[id[2]][1] - spins[id[4]][1] -
                            spins[id[9]][1] - spins[id[11]][1] +
                                      # yz-plane
                            spins[id[3]][1] + spins[id[6]][1] +
                            spins[id[7]][1] + spins[id[10]][1]
                        )
                    )
                )
            end) #--------------------------------------------------------------

            #=
            g-paths
            original: ∑∑∑ xxy + xyx + yxx - yyy
            implies: ∑i ∑j(i) ∑k(j)=/=i xi xj xk + xi yj xk + ...

            affected spins: o-a-b (= b-a-o), a-o-b, where o is flipped
            Notation: dx is flipped, sum over paths (o-a-b, a-o-b) implied
            o-a-b: ∑∑ dx(xy)   + dx(yx)   + dy(xx)   - dy(yy)  +
            a-o-b: ∑∑ (x)dx(y) + (x)dy(x) + (y)dx(x) - y(dy)y
            notation: 1 means o-a-b term, 2 means a-o-b term
            dx * ∑∑ [(xy)1 + (yx)1 + (xy)2 + (yx)2]
            dy * ∑∑ [(xx)1 - (yy)1 + (xx)2 - (yy)2]
            Note symmetry - we can loop over (o-a-b) and (a-o-b) terms here
            further bracketing:
            dx * ∑ [x1 * ∑ (y1 + y2) + y1 * ∑ (x1 + x2)]
            dy * ∑ [x1 * ∑ (x1 + x2) - y1 * ∑ (y1 + y2)]
            Note dublicate terms - we only need to compute two inner brackets
            =#
            # g terms
            $(dog && quote #----------------------------------------------------
                g_x__ = 0.  # ∑ [x1 * g___y + y1 * g___x]
                g_y__ = 0.  # ∑ [x1 * g___x - y1 * g___y]
                for gedge in n.gpaths   # 4 elements
                    g___x = 0.  # ∑ (x1 + x2)
                    g___y = 0.  # ∑ (y1 + y2)
                    for bi in gedge.bs  # 6 5 4 3 elements
                        @inbounds s2 = spins[bi]
                        @fastmath @inbounds g___x += s2[1]
                        @fastmath @inbounds g___y += s2[2]
                    end
                    @inbounds s1 = spins[gedge.a]
                    @fastmath @inbounds g_x__ += s1[1] * g___y + s1[2] * g___x
                    @fastmath @inbounds g_y__ += s1[1] * g___x - s1[2] * g___y
                end
                # commit to dE
                @fastmath @inbounds dE += param.g *
                    (delta_s[1] * g_x__ + delta_s[2] * g_y__)
            end) #--------------------------------------------------------------

            # h
            $(doh && quote #----------------------------------------------------
                @fastmath dE -= dot(param.h, delta_s)
            end) #--------------------------------------------------------------

            return dE
        end
    end
end


################################################################################
#### total Energy functions
################################################################################

# NOTE
# The order of calculations is not optimized here. This is fine because the
# function only gets called once. This also helps verifying the correctness of
# totalEnergy and deltaEnergy.



"""
    totalEnergy(sgraph, spins, parameters)

Calculates the total energy of the current system.
"""
function totalEnergy(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        param::Parameters
    )
    E = 0.
    for e in sgraph.second
        # second enighbor - J2
        E += param.J2[1] * (
            spins[e.n1][1] * spins[e.n2][1] +
            spins[e.n1][2] * spins[e.n2][2]
        ) + param.J2[2] * spins[e.n1][3] * spins[e.n2][3]
        # zeta
        if e.plane == :xy
            E += 2.0 * param.zeta * (
                spins[e.n1][1] * spins[e.n2][1] -
                spins[e.n1][2] * spins[e.n2][2]
            )
        elseif e.plane == :xz
            E += param.zeta * (
                -(
                    spins[e.n1][1] * spins[e.n2][1] -
                    spins[e.n1][2] * spins[e.n2][2]
                ) - sqrt3 * (
                    spins[e.n1][1] * spins[e.n2][2] +
                    spins[e.n1][2] * spins[e.n2][1]
                )
            )
        elseif e.plane == :yz
            E += param.zeta * (
                -(
                    spins[e.n1][1] * spins[e.n2][1] -
                    spins[e.n1][2] * spins[e.n2][2]
                ) + sqrt3 * (
                    spins[e.n1][1] * spins[e.n2][2] +
                    spins[e.n1][2] * spins[e.n2][1]
                )
            )
        else
            error("Second neighbor plane defined incorrectly as $(e.plane).")
        end
    end

    # J3
    for e in sgraph.third
        E += param.J3[1] * (
            spins[e.n1][1] * spins[e.n2][1] +
            spins[e.n1][2] * spins[e.n2][2]
        ) + param.J3[2] * spins[e.n1][3] * spins[e.n2][3]
    end

    # J1, K
    for i in eachindex(sgraph.first)
        e = sgraph.first[i]
        E += param.J1[1] * e.xy + param.J1[2] * e.z
        for p in sgraph.paths[i]
            E += param.K * (e.xy * p.z + e.z * p.xy)
        end
    end

    # g
    # Checked:
    # factor 0.5 (overcounting edges)
    # number of edges (12x a-b-c, 2x (overcounting) 6x c-a-b)
    # edges/paths for first node correct
    for ei in eachindex(sgraph.first)
        e12 = sgraph.first[ei]
        for e23 in sgraph.nodes[e12.n2].first   # LID: 1 -> 2 -> 1
            e12 == e23 && continue
            n3 = e23.n1 != e12.n2 ? e23.n1 : e23.n2
            E += 0.5param.g * (
                spins[e12.n1][1] * spins[e12.n2][1] * spins[n3][2] +
                spins[e12.n1][1] * spins[e12.n2][2] * spins[n3][1] +
                spins[e12.n1][2] * spins[e12.n2][1] * spins[n3][1] -
                spins[e12.n1][2] * spins[e12.n2][2] * spins[n3][2]
            )
        end
        for e23 in sgraph.nodes[e12.n1].first   # LID: 2 -> 1 -> 2
            e12 == e23 && continue              # overcounting cause 2 <- 1 <- 2
            n3 = e23.n1 != e12.n1 ? e23.n1 : e23.n2
            E += 0.5param.g * (
                spins[e12.n2][1] * spins[e12.n1][1] * spins[n3][2] +
                spins[e12.n2][1] * spins[e12.n1][2] * spins[n3][1] +
                spins[e12.n2][2] * spins[e12.n1][1] * spins[n3][1] -
                spins[e12.n2][2] * spins[e12.n1][2] * spins[n3][2]
            )
        end
    end

    # h
    for S in spins
        E -= dot(param.h, S)
    end

    E
end


# NOTE DEPRECATED

function totalEnergy(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Vector{Tuple{Float64, Float64}},
        h::Point3{Float64}=Point3(0.)
    )
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


# DEPRECATED
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
