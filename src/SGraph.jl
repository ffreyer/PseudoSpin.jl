# TODO TODO TODO TODO TODO TODO TODO TODO
# This stuff may fail for small lattices (L = 1 & 2)
# TODO TODO TODO TODO TODO TODO TODO TODO

abstract type SEdge end

# Edge for first order neighbours & paths
# (which are made up of first order edges)
# preallocating scalar products (xy, z) should speed up path calculations
type SEdge1 <: SEdge
    xy::Float64
    z::Float64

    n1::Int64
    n2::Int64
end

# Edge for second (...) order neighbours.
# pre-allocating scalar products would be a net loss here
# (+45 allocations, +21 additions if my math is correct)
immutable SEdge2 <: SEdge
    n1::Int64
    n2::Int64
end


immutable SNode
    first::Vector{SEdge1}
    second::Vector{Int64}
    # includes all additonal edges in paths containing neighbours[1][j]
    paths::Vector{Vector{SEdge1}}
    # gpaths contains NN edges a-b for x-a-b terms and NNN edges a-b for a-x-b
    gpaths::Tuple{Vector{SEdge1}, Vector{SEdge2}}
end


type SGraph
    N_nodes::Int64
    K_edges::Int64
    K_paths::Int64

    nodes::Vector{SNode}

    first::Vector{SEdge1}
    second::Vector{SEdge2}
    paths::Vector{Vector{SEdge1}} # [continues edges[1][1], edges[1][2], ...]
end



################################################################################
### Functions for Simulation Graphs
################################################################################


function ==(e1::SEdge, e2::SEdge)
    if e1.n1 == e2.n1
        return e1.n2 == e2.n2
    elseif e1.n1 == e2.n2
        return e1.n2 == e2.n1
    end
    false
end


function in(edge::SEdge, v::AbstractArray{SEdge})
    for e in v
        if edge == e
            return true
        end
    end
    false
end


################################################################################
### Graph Utility Functions
################################################################################


function define_first!{T <: SEdge}(A::Vector{T}, item::T)
    for i in eachindex(A)
        if !isdefined(A, i)
            A[i] = item
            break
        end
    end
    nothing
end

function define_first!(A::Vector{Int64}, item::Int64)
    for i in eachindex(A)
        if A[i] == 0
            A[i] = item
            break
        end
    end
    nothing
end



# fills SNode.neighbours with relevant nodes if they exist (no lattice-edge recursion)
"""
    connect_nodes!(
        rgraph::RGraph, nodes::Vector{SNode},
        first::Vector{SEdge1}, second::Vector{SEdge2}, paths::Vector{Vector{SEdge1}},
        IDs::Vector{Tuple{Vec3i, Int64}}, flat_index::Function
    )

Fills *first*, *second* and *paths* with the relevant edges. The edge lists
should be empty.
"""
function connect_nodes!(
        rgraph::RGraph,
        nodes::Vector{SNode},
        first::Vector{SEdge1},
        second::Vector{SEdge2},
        paths::Vector{Vector{SEdge1}},
        IDs::Vector{Tuple{Vec3i, Int64}},
        flat_index::Function
    )

    # N: number of created nodes
    # M: whatever find needs
    # E: number of edges per node
    # O(N * M * E)
    # NN
    for i in 1:Int64(length(nodes)) # IndexT
        rnode = rgraph.nodes[IDs[i][2]]
        lvl = 1
        for re in rnode.edges[lvl]
            l = flat_index(IDs[i][1] + re.dir, re.to)

            if l != 0
                if lvl == 1     # <=> NN
                    e = SEdge1(-1., -1., i, l)
                    k = findfirst(first, e)
                    if k != 0
                        push!(nodes[i].first, first[k])
                    else
                        push!(nodes[i].first, e)
                        push!(first, e)
                    end
                # elseif lvl == 2     # <=> NNN
                #     e = SEdge2(i, l)
                #     k = findfirst(second, e)
                #     if k != 0
                #         push!(nodes[i].second, l)
                #     else
                #         push!(nodes[i].second, l)
                #         push!(second, e)
                #     end
                else
                    throw(ErrorException(
                        "More than second neighbours not implemented!"
                    ))
                end
            end
        end
    end

    # NNN
    if length(rgraph.nodes[1].distances) >= 2
        for i in 1:Int64(length(nodes)) # IndexT
            rnode = rgraph.nodes[IDs[i][2]]
            lvl = 2
            for re in rnode.edges[lvl]

                l = flat_index(IDs[i][1] + re.dir, re.to)

                if l != 0
                    # if lvl == 1
                    #     e = SEdge1(-1., -1., i, l)
                    #     k = findfirst(first, e)
                    #     if k != 0
                    #         push!(nodes[i].first, first[k])
                    #     else
                    #         push!(nodes[i].first, e)
                    #         push!(first, e)
                    #     end
                    # else
                    if lvl == 2
                        e = SEdge2(i, l)
                        k = findfirst(second, e)
                        if k != 0
                            push!(nodes[i].second, l)
                        else
                            push!(nodes[i].second, l)
                            push!(second, e)
                        end
                    else
                        throw(ErrorException(
                            "More than second neighbours not implemented!"
                        ))
                    end
                end

            end
        end
    end

    # Array with a spot for each lvl 1 edge
    append!(paths, [SEdge1[] for _ in first])

    # if paths exist
    if length(rgraph.nodes[1].edges) > length(rgraph.nodes[1].distances) # 3

        for ei in eachindex(first)
            # n node, e edge, rp RPath, se SEdge, i index
            # get edge, starting node, rnode
            e = first[ei]
            n1 = nodes[e.n1]
            n2 = nodes[e.n2]
            rnode = e.n1 <= div(length(nodes), 2) ? rgraph.nodes[1] : rgraph.nodes[2] # TODO check this

            # get neighbour index from nodes   ..., path?
            nei1 = findfirst(n1.first, e)
            nei2 = findfirst(n2.first, e)

            # get paths in rnode starting edge e. (This works because node.first
            # is sorted like rnode.edges[1] & rnode.edges[3])
            # this pushes all paths starting with (n1 -> n2)
            # reverse paths, such as (l -> k) (n2 -> n1) won't be pushed
            for rp in rnode.edges[end][(nei1-1) * 9 + 1 : nei1 * 9]

                ni3 = flat_index(IDs[e.n1][1] + rp.dirs[2], rp.tos[2])  # TODO: check if not 0 (for open bounds)
                ni4 = flat_index(IDs[e.n1][1] + rp.dirs[3], rp.tos[3])  # TODO: check if not 0 (for open bounds)
                if (ni3 == 0) || (ni4 == 0); continue; end
                ei2 = findfirst(first, SEdge1(-1., -1, ni3, ni4))

                # save in nodes
                push!(n1.paths[nei1], first[ei2])
                push!(n2.paths[nei2], first[ei2])

                # save in edges
                if ei2 > ei
                    push!(paths[ei], first[ei2])
                end
            end

            # This is required to catch all nodes (including x -> this -> y -> z)
            rnode = e.n2 <= div(length(nodes), 2) ? rgraph.nodes[1] : rgraph.nodes[2] # TODO check this

            # this does all paths starting with (n2 -> n1)
            for rp in rnode.edges[end][(nei2-1) * 9 + 1 : nei2 * 9]

                # get Second Edge in path
                ni3 = flat_index(IDs[e.n2][1] + rp.dirs[2], rp.tos[2])  # TODO: check if not 0 (for open bounds)
                ni4 = flat_index(IDs[e.n2][1] + rp.dirs[3], rp.tos[3])  # TODO: check if not 0 (for open bounds)
                if (ni3 == 0) || (ni4 == 0); continue; end
                ei2 = findfirst(first, SEdge1(-1., -1., ni3, ni4))

                # save in nodes
                push!(n1.paths[nei1], first[ei2])
                push!(n2.paths[nei2], first[ei2])

                # save in edges
                if ei2 > ei
                    push!(paths[ei], first[ei2])
                end
            end
        end
    end

    for (xi, x) in enumerate(nodes)
        # x - a - b paths
        for skipped_edge in x.first
            ai = skipped_edge.n1 == xi ? skipped_edge.n2 : skipped_edge.n1
            for e in nodes[ai].first
                e == skipped_edge && continue # no reversal
                push!(x.gpaths[1], e)
            end
        end

        # a - x - b
        for i in 1:4        # picks a
            axe = x.first[i]
            ai = axe.n1 == xi ? axe.n2 : axe.n1
            for j in i+1:4  # picks b
                xbe = x.first[j]
                bi = xbe.n1 == xi ? xbe.n2 : xbe.n1
                abe = SEdge2(ai, bi)
                abi = findfirst(second, abe)
                if abi == 0
                    push!(x.gpaths[2], abe)
                else
                    push!(x.gpaths[2], second[abi])
                end
            end
        end
    end

    nothing
end


################################################################################
### Graph Constructors
################################################################################

# About "index order"
# With that I want to describe when which index changes. [1, 2, 3] would imply
# that the first index changes the fastes, while the third changes the slowest.

# Create points based on their Bravais lattice. This may result in incomplete
# regions for combined lattices

"""
    Basisfill(rgraph, N[; border = :periodic])

Generates an extensive graph with 2*N^3 nodes and edges connecting neighbors.
The lattice is shaped like the primitive unit cell. Returns the graph, a function
mapping matrix indices to list indices and a function mapping flat to matrix
indices.
"""
function Basisfill(rgraph::RGraph, N::Integer; border::Symbol=:periodic)
    # Indexing
    # index order [(3, 2, 1), 4]
    LSize = N^3
    flattening_vector = Vec3i(N^2, N, 1)

    # Build indexing functions

    # tested completely for 32-edge lattice
    function periodic(uvw, LID) # (uvw::Vec3i, LID::Int64)
        # index order: [(3, 2, 1), 4]
        x = 0
        y = 0
        z = 0

        # This does the border recursion
        uvw[1] < 0  && (x = N)
        uvw[1] >= N && (x = - N)
        uvw[2] < 0  && (y = N)
        uvw[2] >= N && (y = - N)
        uvw[3] < 0  && (z = N)
        uvw[3] >= N && (z = - N)
        # print(x, ", ", y, ", ", z, " \t ")

        1 + LSize * (LID-1) + dot(flattening_vector, uvw + Vec3i(x, y, z)) # IndexT
    end

    # Should always be in bounds; should fit loop
    function open(uvw, LID) # (uvw::Vec3i, LID::Int64)
        # uvw, LID = args
        # index order: [(3, 2, 1), 4]

        (0 > uvw[1] || uvw[1] >= N) && return 0
        (0 > uvw[2] || uvw[2] >= N) && return 0
        (0 > uvw[3] || uvw[3] >= N) && return 0

        1 + LSize * (LID-1) + dot(flattening_vector, uvw) # IndexT
    end

    if border == :periodic
        flat_index = periodic
    elseif border == :open
        flat_index = open
    else
        throw(ErrorException("border not available"))
    end

    # Decomposition should always be the same, right?
    function decompose_index(i::Int64)
        # index order [(3, 2, 1), 4]
        LID = div(i-1, LSize) + 1                       # 1-based
        sublattice_index = (i-1) % LSize                # 0-based

        u = div(sublattice_index, flattening_vector[1]) # 0-based
        r = sublattice_index % flattening_vector[1]     # 0-based
        v = div(r, flattening_vector[2])                # 0-based
        w = r % flattening_vector[2]                    # 0-based

        Vec3i(u, v, w), LID
    end


    # Initialisation
    # The order of this arrays is irrelevant. The nodes will be set according
    # to IDs and the indice functions
    nodes = [
        SNode(
            SEdge1[],
            Int64[],
            [SEdge1[] for _ in 1:4],
            (SEdge1[], SEdge2[])
        ) for __ in 1:N^3, _ in 1:length(rgraph.nodes)
    ][:]

    # The ordering of this HAS to match the indice function!
    # index order: [(3, 2, 1), 4]

    IDs = [
        (Vec3i(u, v, w), i)
        for i in 1:length(rgraph.nodes)
        for u in 0:N-1 for v in 0:N-1 for w in 0:N-1
    ]       # the last loop runs through first

    # populate neighbours
    first = SEdge1[] # [SEdge[], SEdge[]]
    second = SEdge2[]  # [SEdge[], SEdge[]]
    paths = Vector{SEdge1}[]
    connect_nodes!(rgraph, nodes, first, second, paths, IDs, flat_index)


    g = SGraph(
        length(nodes), length(first) + length(second), length(paths),
        nodes, first, second, paths
    )

    # sim, _, lattice_indices
    g, flat_index, decompose_index
end


################################################################################
### Other Constructors
################################################################################


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

    Point3{Float64}(
        st * cos(phi),
        st * sin(phi),
        ct
    )
end

function rand_spin(N::Int64)
    phis = 2 * pi * rand(Float64, N)
    cts = 2 * rand(Float64, N) - 1
    sts = sqrt.(1. - cts .* cts) # sin(acos(cts)) # max(0., )

    [Point3{Float64}(
        sts[i] .* cos.(phis[i]),
        sts[i] .* sin.(phis[i]),
        cts[i]
    ) for i in 1:N]
end


# Not as fast as rand_spin(N) iirc
# function rand_spin!(A::Vector{Point3{Float64}})
#     N = length(A)
#     phis = 2 * pi * rand(Float64, N)
#     cts = 2 * rand(Float64, N) - 1
#     sts = sin(acos(cts))
#
#     @inbounds for i in 1:N
#         A[i] = Point3{Float64}(
#             sts[i] * cos(phis[i]),
#             sts[i] * sin(phis[i]),
#             cts[i]
#         )
#     end
#
#     nothing
# end


# unused?
# function get_positions(rgraph::RGraph, sgraph::SGraph)
#     map(1:sgraph.N_nodes) do i # IndexT
#         uvw, LID = lattice_indices(i)
#         uvw * rgraph.nodes[LID].bravais
#     end
# end
