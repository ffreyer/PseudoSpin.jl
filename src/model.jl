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
        Js::Vector{Float64},
        h::Point3{Float64}=Point3(0.)
    )

    E = 0.

    for e in sgraph.second
        E += Js[2] * dot(spins[e.n1], spins[e.n2])
    end

    for i in eachindex(sgraph.first)
        e = sgraph.first[i]
        E += Js[1] * (e.xy + e.z)
        for p in sgraph.paths[i]
            E += Js[3] * (e.xy + e.z) * Js[4] * (p.xy + p.z)
        end
    end

    if !(Tuple(h) == (0., 0., 0.))
        for S in spins
            E -= dot(h, S)
        end
    end

    E
end


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


# total Energy for (Js, h, g)
function totalEnergy(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Vector{Tuple{Float64, Float64}},
        h::Point3{Float64},
        g::Float64
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
#### single spin flip kernels
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
#### delta Energy Functions
################################################################################
# Note:
# I tried to optimize this stuff as much as I could...


# Energy shift for isotropic Js
# This function is not really used anymore and runs about as fast as the
# anisotropic equivalent below.
"""
    deltaEnergy(node, spins, index, new_spin, Js, xys, zs, h)

Calculates the energy difference for a spin at some index changing to new_spin.
"""
function deltaEnergy(
        n::SNode,                           # the changing node
        spins::Vector{Point3{Float64}},     # current spins
        i::Int64,                           # index of node/changed spin
        new_spin::Point3{Float64},          # new spin
        Js::Vector{Float64},                # bond parameters
        xys::Vector{Float64},               # new scalar products
        zs::Vector{Float64},                # new scalar products
        h::Point3{Float64}=Point3(0.)       # external field
    )

    # Runtime:
    # xys, zs: 4x 3 multiplication, 4 additions -> (12, 4) (used for NN, paths)

    # NN: [12, 4] + (1, 4x 4) = (1, 16) + [12, 4]
    # NNN: 12x (3, 3) + (0, 3) = (36, 39)
    # paths: [12, 4] + 36x (0, 2) + 4x (1, 4) + (2, 0) = (6, 88) + [12, 4]
    #                                       but [12, 4] is reused from NN
    # field: (3, 6)
    # Assuming multiplication is about 4 times as expensive as addition/subtraction
    # scalar products: 52*
    # NN: 20 [+ 52]*
    # NNN: 183
    # paths: 112 [+ 52]*
    # field: 18

    # calculate energy difference
    dE = 0.
    temp = 0.
    for (xy, z, j) in zip(xys, zs, eachindex(n.first))
        @inbounds e = n.first[j]
        temp2 = 0.
        @fastmath dE += xy - e.xy + z - e.z
        for p in n.paths[j]
            @fastmath temp2 += p.xy + p.z # p.xy and p.z don't change
        end
        @fastmath temp += temp2 * (xy - e.xy + z - e.z)
    end
    @fastmath @inbounds dE *= Js[1]
    @fastmath @inbounds dE += Js[3] * Js[4] * temp

    temp = 0.
    @fastmath @inbounds delta_s = new_spin - spins[i]
    for j in n.second
        @fastmath @inbounds temp += dot(delta_s, spins[j])
    end
    @fastmath @inbounds dE += Js[2] * temp

    if !(Tuple(h) == (0., 0., 0.))
        @fastmath @inbounds dE -= dot(h, (new_spin - spins[i]))
    end

    return dE
end


# Energy shift for anisotropic Js
# Checked vs totalEnergy() for anisotropic J3, J4 (aka K) in reordered 144-term
# form. Error ~1e-12% (1e-14 as factor)
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
        2 * (
            Js[3][1] * Js[4][1] * xy3 +
            Js[3][2] * Js[4][2] * z3
        ) + (
            Js[3][1] * Js[4][2] +
            Js[3][2] * Js[4][1]
        ) * xyz3 +
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
# Isotropic version, no update to total Energy
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
        Js::Vector{Float64},
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


# Anisotropic version, no update to total Energy
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


# isotropic Js, returns updated E_tot
function kernel(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        i::Int64,
        new_spin::Point3{Float64},
        E_tot::Float64,
        Js::Vector{Float64},
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


################################################################################
#### sweep functions
################################################################################


"""
    sweep(sgraph, spins[, E_tot], Js, beta, h)

Attempts as many spin flips as there are sites in the lattice. The current total
energy E_tot can be given optionally. If done so, it will be updated after each
successful spin flip and returned in the end.
"""
function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        beta::Float64,
        h::Point3{Float64}=Point3(0.)
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
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
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
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        beta::Float64,
        h::Point3{Float64}=Point3(0.)
    )

    for (i, new_spin) in zip(
            rand(1:sgraph.N_nodes, sgraph.N_nodes),
            rand_spin(sgraph.N_nodes)
        )
        E_tot = kernel(sgraph, spins, i, new_spin, E_tot, Js, beta, h)
    end

    E_tot
end


function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
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


# (Js, h, g) w/o E update
function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
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


# (UJs, h, g) w/ E update
function sweep(
        sgraph::SGraph,
        spins::Vector{Point3{Float64}},
        E_tot::Float64,
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
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
