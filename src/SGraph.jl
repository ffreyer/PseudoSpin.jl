abstract type SEdge end

# Edge for first order neighbours & paths
# (which are made up of first order edges)
# preallocating scalar products (xy, z) should speed up path calculations
mutable struct SEdge1 <: SEdge
    xy::Float64
    z::Float64

    n1::Int64
    n2::Int64
end

# Edge for second (...) order neighbours.
# pre-allocating scalar products would be a net loss here
# (+45 allocations, +21 additions if my math is correct)
struct SEdge2 <: SEdge
    n1::Int64
    n2::Int64
    plane::Symbol
end

struct GEdge <: SEdge
    a::Int64
    bs::Vector{Int64}
end


struct SNode
    first::Vector{SEdge1}
    second::Vector{Int64}
    third::Vector{Int64}
    # includes all additonal edges in paths containing neighbours[1][j]
    paths::Vector{Vector{SEdge1}}
    gpaths::Vector{GEdge}
end


mutable struct SGraph
    N_nodes::Int64
    K_edges::Int64
    K_paths::Int64

    nodes::Vector{SNode}

    first::Vector{SEdge1}
    second::Vector{SEdge2}
    third::Vector{SEdge2}
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


function define_first!(A::Vector{T}, item::T) where {T <: SEdge}
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
        rgraph::RGraph,
        nodes::Vector{SNode},
        first::Vector{SEdge1},
        second::Vector{SEdge2},
        third::Vector{SEdge2},
        paths::Vector{Vector{SEdge1}},
        IDs::Vector{Tuple{SVector{3, Int64}, Int64}},
        flat_index::Function
    )

Fills *first*, *second* and *paths* with the relevant edges. The edge lists
should be empty.
"""
function connect_nodes!(
        rgraph::RGraph,
        nodes::Vector{SNode},
        first::Vector{SEdge1},
        second::Vector{SEdge2},
        third::Vector{SEdge2},
        paths::Vector{Vector{SEdge1}},
        IDs::Vector{Tuple{SVector{3, Int64}, Int64}},
        flat_index::Function
    )

    # N: number of created nodes
    # M: whatever find needs
    # E: number of edges per node
    # O(N * M * E)

    # NNN
    for lvl in eachindex(rgraph.nodes[1].distances)
        for i in 1:Int64(length(nodes)) # IndexT
            rnode = rgraph.nodes[IDs[i][2]]
            # lvl = 2
            for (ri, re) in enumerate(rnode.edges[lvl])

                l = flat_index(IDs[i][1] + re.dir, re.to)

                if l != 0
                    if lvl == 1
                        e = SEdge1(-1., -1., i, l)
                        k = findfirst(isequal(e), first)
                        if k != nothing
                            push!(nodes[i].first, first[k])
                        else
                            push!(nodes[i].first, e)
                            push!(first, e)
                        end
                    elseif lvl == 2
                        plane = [
                            :xy, :xz, :yz, :xz,
                            :xy, :yz, :yz, :xy,
                            :xz, :yz, :xz, :xy
                        ][ri]
                        e = SEdge2(i, l, plane)
                        k = findfirst(isequal(e), second)
                        if k != nothing
                            push!(nodes[i].second, l)
                        else
                            push!(nodes[i].second, l)
                            push!(second, e)
                        end
                    elseif lvl == 3
                        e = SEdge2(i, l, :None)
                        k = findfirst(isequal(e), third)
                        if k != nothing
                            push!(nodes[i].third, l)
                        else
                            push!(nodes[i].third, l)
                            push!(third, e)
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
            # TODO is this ever nothing/0?
            nei1 = something(findfirst(isequal(e), n1.first), 0)
            nei2 = something(findfirst(isequal(e), n2.first), 0)

            # get paths in rnode starting edge e. (This works because node.first
            # is sorted like rnode.edges[1] & rnode.edges[3])
            # this pushes all paths starting with (n1 -> n2)
            # reverse paths, such as (l -> k) (n2 -> n1) won't be pushed
            for rp in rnode.edges[end][(nei1-1) * 9 + 1 : nei1 * 9]

                ni3 = flat_index(IDs[e.n1][1] + rp.dirs[2], rp.tos[2])  # TODO: check if not 0 (for open bounds)
                ni4 = flat_index(IDs[e.n1][1] + rp.dirs[3], rp.tos[3])  # TODO: check if not 0 (for open bounds)
                if (ni3 == 0) || (ni4 == 0); continue; end
                ei2 = findfirst(isequal(SEdge1(-1., -1, ni3, ni4)), first)

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
                ei2 = findfirst(isequal(SEdge1(-1., -1, ni3, ni4)), first)


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
            bs = Int64[]
            for e in nodes[ai].first
                e == skipped_edge && continue # no reversal
                bi = e.n1 == ai ? e.n2 : e.n1
                push!(bs, bi)
            end
            push!(x.gpaths, GEdge(ai, bs))
        end

        # a - x - b
        l = length(x.first)
        for i in 1:l-1        # picks a
            axe = x.first[i]
            ai = axe.n1 == xi ? axe.n2 : axe.n1
            bs = Int64[]
            for j in i+1:l  # picks b
                xbe = x.first[j]
                bi = xbe.n1 == xi ? xbe.n2 : xbe.n1
                push!(bs, bi)
            end
            append!(x.gpaths[i].bs, bs)
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
    @assert N > 2 "Lattice may not generate correctly for L < 3."
    # Indexing
    # index order [(3, 2, 1), 4]
    LSize = N^3
    flattening_vector = SVector{3, Int64}(N^2, N, 1)

    # Build indexing functions

    # tested completely for 32-edge lattice
    function periodic(uvw, LID) # (uvw::SVector{3, Int64}, LID::Int64)
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

        1 + LSize * (LID-1) + dot(flattening_vector, uvw + SVector{3, Int64}(x, y, z)) # IndexT
    end

    # Should always be in bounds; should fit loop
    function open(uvw, LID) # (uvw::SVector{3, Int64}, LID::Int64)
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

        SVector{3, Int64}(u, v, w), LID
    end


    # Initialisation
    # The order of this arrays is irrelevant. The nodes will be set according
    # to IDs and the indice functions
    nodes = [
        SNode(
            SEdge1[],
            Int64[],
            Int64[],
            [SEdge1[] for _ in 1:4],
            GEdge[]
        ) for __ in 1:N^3, _ in 1:length(rgraph.nodes)
    ][:]

    # The ordering of this HAS to match the indice function!
    # index order: [(3, 2, 1), 4]

    IDs = [
        (SVector{3, Int64}(u, v, w), i)
        for i in 1:length(rgraph.nodes)
        for u in 0:N-1 for v in 0:N-1 for w in 0:N-1
    ]       # the last loop runs through first

    # populate neighbours
    first = SEdge1[] # [SEdge[], SEdge[]]
    second = SEdge2[]  # [SEdge[], SEdge[]]
    third = SEdge2[]  # [SEdge[], SEdge[]]
    paths = Vector{SEdge1}[]
    connect_nodes!(rgraph, nodes, first, second, third, paths, IDs, flat_index)


    g = SGraph(
        length(nodes),
        length(first) + length(second) + length(third), length(paths),
        nodes, first, second, third, paths
    )

    # sim, _, lattice_indices
    g, flat_index, decompose_index
end


################################################################################
### Other Constructors
################################################################################
