# To reduce the computational cost of the K-term (paths of 4 spins) the results
# of scalar products between neighbors are safed. This function calculates the
# inital values for those. (see also SGraph.jl)
# This has to happen before any energy calculation
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
# function only gets called once.
function totalEnergy(sgraph::SGraph, spins::Vector{Point3{Float64}},
        Js::Vector{Float64}, h::Point3{Float64}=Point3(0.))

    E = 0.

    for e in sgraph.second
        E += Js[2] * dot(spins[e.n1], spins[e.n2])
    end

    for i in eachindex(sgraph.first)
        e = sgraph.first[i]
        E += Js[1] * (e.xy + e.z)
        for p in sgraph.paths[i]
            E += Js[3] * (e.xy + e.z) * (p.xy + p.z)
        end
    end

    if !(Tuple(h) == (0., 0., 0.))
        for S in spins
            E -= dot(h, S)
        end
    end

    E
end


function totalEnergy(sgraph::SGraph, spins::Vector{Point3{Float64}},
        Js::Vector{Tuple{Float64, Float64}}, h::Point3{Float64}=Point3(0.))

    E = 0.
    for e in sgraph.second
        E += Js[2][1] * (
            spins[e.n1][1] * spins[e.n2][1] +
            spins[e.n1][2] * spins[e.n2][2]
        ) + Js[2][2] * spins[e.n1][3] * spins[e.n2][3]
    end

    for i in eachindex(sgraph.first) # eachindex(sgraph.paths)
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


################################################################################
#### single spin flip kernels
################################################################################


# calculates new values to be allocated as scalar product on Edge e
function scalar_prod(e::Edge1, i::Int64, new_spin::Point3{Float64},
        spins::Vector{Point3{Float64}})

    j = i == e.n1 ? e.n2 : e.n1
    @fastmath @inbounds return new_spin[1] * spins[j][1] + new_spin[2] * spins[j][2], new_spin[3] * spins[j][3]
end


# calculates all new scalar products for a node at index i with spin new_spin
function generate_scalar_products(sgraph::SGraph, spins::Vector{Point3{Float64}},
        i::Int64, new_spin::Point3{Float64})

    xys = Array(Float64, 4)
    zs = Array(Float64, 4)

    for j in eachindex(sgraph.nodes[i].first)
        @inbounds xys[j], zs[j] = scalar_prod(sgraph.nodes[i].first[j], i, new_spin, spins)
    end

    xys, zs
end


# pushes new scalar products to the relevant edges of Node n
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
function deltaEnergy(n::SNode, spins::Vector{Point3{Float64}}, i::Int64,
        new_spin::Point3{Float64}, Js::Vector{Float64},
        xys::Vector{Float64}, zs::Vector{Float64}, h::Point3{Float64}=Point3(0.))

    # println("Don't use this! (isotropic Js)")
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
        # @fastmath @inbounds dE += Js[2] * dot(delta_s, spins[j])
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
function deltaEnergy(n::SNode, spins::Vector{Point3{Float64}}, i::Int64,
        new_spin::Point3{Float64}, Js::Vector{Tuple{Float64, Float64}},
        xys::Vector{Float64}, zs::Vector{Float64}, h::Point3{Float64}=Point3(0.))

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


################################################################################
#### single spin flip kernels
################################################################################


# This attempts a single spin flip for a given new_spin and array position i.
# Isotropic version, no update to total Energy
function kernel(sgraph::SGraph, spins::Vector{Point3{Float64}}, i::Int64,
        new_spin::Point3{Float64},
        Js::Vector{Float64}, beta::Float64, h::Point3{Float64}=Point3(0.)
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
function kernel(sgraph::SGraph, spins::Vector{Point3{Float64}}, i::Int64,
        new_spin::Point3{Float64},
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64, h::Point3{Float64}=Point3(0.))

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
function kernel(sgraph::SGraph, spins::Vector{Point3{Float64}}, i::Int64,
        new_spin::Point3{Float64}, E_tot::Float64,
        Js::Vector{Float64}, beta::Float64, h::Point3{Float64}=Point3(0.))
        #Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}})

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
function kernel(sgraph::SGraph, spins::Vector{Point3{Float64}}, i::Int64,
        new_spin::Point3{Float64}, E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}},
        beta::Float64, h::Point3{Float64}=Point3(0.)) #new_spin::Point3{Float64},

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
function kernel(sgraph::SGraph, spins::Vector{Point3{Float64}}, i::Int64,
        new_spin::Point3{Float64}, E_tot::Float64,
        Js::Vector{Tuple{Float64, Float64}}, h::Point3{Float64}=Point3(0.)) #new_spin::Point3{Float64},

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


################################################################################
#### sweep functions
################################################################################


function sweep(sgraph::SGraph, spins::Vector{Point3{Float64}},
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        beta::Float64, h::Point3{Float64}=Point3(0.))

    for (i, new_spin) in zip(rand(1:sgraph.N_nodes, sgraph.N_nodes), rand_spin(sgraph.N_nodes)) # IndexT
        kernel(sgraph, spins, i, new_spin, Js, beta, h)
    end

    nothing
end


function sweep(sgraph::SGraph, spins::Vector{Point3{Float64}},
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        h::Point3{Float64}=Point3(0.))

    for (i, new_spin) in zip(rand(1:sgraph.N_nodes, sgraph.N_nodes), rand_spin(sgraph.N_nodes)) # IndexT
        kernel(sgraph, spins, i, new_spin, Js, h)
    end

    nothing
end


function sweep(sgraph::SGraph, spins::Vector{Point3{Float64}}, E_tot::Float64,
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        beta::Float64, h::Point3{Float64}=Point3(0.))

    for (i, new_spin) in zip(rand(1:sgraph.N_nodes, sgraph.N_nodes), rand_spin(sgraph.N_nodes)) # IndexT
        E_tot = kernel(sgraph, spins, i, new_spin, E_tot, Js, beta, h)
    end

    E_tot
end


function sweep(sgraph::SGraph, spins::Vector{Point3{Float64}}, E_tot::Float64,
        Js::Union{Vector{Float64}, Vector{Tuple{Float64, Float64}}},
        h::Point3{Float64}=Point3(0.))

    for (i, new_spin) in zip(rand(1:sgraph.N_nodes, sgraph.N_nodes), rand_spin(sgraph.N_nodes)) # IndexT
        E_tot = kernel(sgraph, spins, i, new_spin, E_tot, Js, h)
    end

    E_tot
end
