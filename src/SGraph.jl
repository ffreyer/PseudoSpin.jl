# TODO TODO TODO TODO TODO TODO TODO TODO
# This stuff may fail for small lattices (L = 1 & 2)
# TODO TODO TODO TODO TODO TODO TODO TODO

abstract SEdge
# typealias IndexT Int64 # UInt16 # <-- make this UInt16, no array will be longer than 65 535...
# L = 16 has 65 536 lvl 1 and lvl 2 edges

# Edge for first order neighbours & paths
# (which are made up of first order edges)
# preallocating scalar products (xy, z) should speed up path calculations
type Edge1 <: SEdge
    xy::Float64
    z::Float64

    n1::Int64
    n2::Int64
end

# Edge for second (...) order neighbours.
# pre-allocating scalar products would be a net loss here
# (+45 allocations, +21 additions if my math is right)
immutable Edge2 <: SEdge
    n1::Int64
    n2::Int64
end


immutable SNode
    # neighbours::Vector{Vector{SEdge}}
    first::Vector{Edge1}
    second::Vector{Int64}
    paths::Vector{Vector{Edge1}}        # includes all additonal edges in paths containing neighbours[1][j]
end


type SGraph
    N_nodes::Int64
    K_edges::Int64
    K_paths::Int64

    nodes::Vector{SNode}

    first::Vector{Edge1}
    second::Vector{Edge2}
    paths::Vector{Vector{Edge1}} # [continues edges[1][1], edges[1][2], ...]
end



################################################################################
### Functions for Simulation Graphs
################################################################################

# SNode(rnode::RNode) = SNode(Vector{Int64}[])


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
# something like this


# fills SNode.neighbours with relevant nodes if they exist (no lattice-edge recursion)
function connect_nodes!(rgraph::RGraph, nodes::Vector{SNode},
        first::Vector{Edge1}, second::Vector{Edge2}, paths::Vector{Vector{Edge1}},
        IDs::Vector{Tuple{Vec3i, Int64}}, flat_index::Function)

    # N: number of created nodes
    # K: number of neighbours (10-20 per order)
    # M: whatever find needs
    # E: number of edges
    # O(N * K * M * E)
    for i in 1:Int64(length(nodes)) # IndexT
        rnode = rgraph.nodes[IDs[i][2]]
        lvl = 1
        #for lvl in eachindex(rnode.distances)
            for re in rnode.edges[lvl]

                l = flat_index(IDs[i][1] + re.dir, re.to)

                if l != 0
                    if lvl == 1
                        e = Edge1(-1., -1., i, l)
                        k = findfirst(first, e)
                        if k != 0
                            # define_first!(nodes[i].first, first[k])
                            push!(nodes[i].first, first[k])
                        else
                            # define_first!(nodes[i].first, e)
                            push!(nodes[i].first, e)
                            push!(first, e)
                        end
                    elseif lvl == 2
                        e = Edge2(i, l)
                        k = findfirst(second, e)
                        if k != 0
                            # define_first!(nodes[i].second, l)
                            push!(nodes[i].second, l)
                        else
                            # define_first!(nodes[i].second, l)
                            push!(nodes[i].second, l)
                            push!(second, e)
                        end
                    else
                        throw(ErrorException("More than second neighbours not implemented!"))
                    end
                end

            end
        #end
    end

    if length(rgraph.nodes[1].distances) >= 2
        for i in 1:Int64(length(nodes)) # IndexT
            rnode = rgraph.nodes[IDs[i][2]]
            lvl = 2
        #for lvl in eachindex(rnode.distances)
            for re in rnode.edges[lvl]

                l = flat_index(IDs[i][1] + re.dir, re.to)

                if l != 0
                    if lvl == 1
                        e = Edge1(-1., -1., i, l)
                        k = findfirst(first, e)
                        if k != 0
                            # define_first!(nodes[i].first, first[k])
                            push!(nodes[i].first, first[k])
                        else
                            # define_first!(nodes[i].first, e)
                            push!(nodes[i].first, e)
                            push!(first, e)
                        end
                    elseif lvl == 2
                        e = Edge2(i, l)
                        k = findfirst(second, e)
                        if k != 0
                            # define_first!(nodes[i].second, l)
                            push!(nodes[i].second, l)
                        else
                            # define_first!(nodes[i].second, l)
                            push!(nodes[i].second, l)
                            push!(second, e)
                        end
                    else
                        throw(ErrorException("More than second neighbours not implemented!"))
                    end
                end

            end
        end
    end

    # Array with a spot for each lvl 1 edge
    append!(paths, [Edge1[] for _ in first])

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
                ei2 = findfirst(first, Edge1(-1., -1, ni3, ni4))

                # save in nodes
                # define_first!(n1.paths[nei1], first[ei2])
                # define_first!(n2.paths[nei2], first[ei2])
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
                ei2 = findfirst(first, Edge1(-1., -1., ni3, ni4))

                # save in nodes
                # define_first!(n1.paths[nei1], first[ei2])
                # define_first!(n2.paths[nei2], first[ei2])
                push!(n1.paths[nei1], first[ei2])
                push!(n2.paths[nei2], first[ei2])

                # save in edges
                if ei2 > ei
                    push!(paths[ei], first[ei2])
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
            Edge1[], Int64[],
            [Edge1[] for _ in 1:4]
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
    first = Edge1[] # [SEdge[], SEdge[]]
    second = Edge2[]  # [SEdge[], SEdge[]]
    paths = Vector{Edge1}[]
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


# Following http://mathworld.wolfram.com/SpherePointPicking.html
# function rand_spin()
#     phi = 2*pi * rand(Float32)
#     theta = acos(2*rand(Float32) - 1)
#     Point3f0(
#         sin(theta) * cos(phi),
#         sin(theta) * sin(phi),
#         cos(theta)
#     )
# end

# more allocations, less calculations
function rand_spin()
    phi = 2 * pi * rand(Float64)
    ct = 2 * rand(Float64) - 1
    st = sqrt(max(0., 1. - ct*ct)) #sin(acos(ct))

    Point3{Float64}(
        st * cos(phi),
        st * sin(phi),
        ct
    )
end


# This is indeed much faster (2-3 times)
function rand_spin(N::Int64) # Union{IndexT, Int64}
    phis = 2 * pi * rand(Float64, N)
    cts = 2 * rand(Float64, N) - 1
    sts = sqrt(max(0., 1. - cts .* cts)) # sin(acos(cts))

    [Point3{Float64}(
        sts[i] * cos(phis[i]),
        sts[i] * sin(phis[i]),
        cts[i]
    ) for i in 1:N]
end


# This should also be faster
function rand_spin!(A::Vector{Point3{Float64}})
    N = length(A)
    phis = 2 * pi * rand(Float64, N)
    cts = 2 * rand(Float64, N) - 1
    sts = sin(acos(cts))

    @inbounds for i in 1:N
        A[i] = Point3{Float64}(
            sts[i] * cos(phis[i]),
            sts[i] * sin(phis[i]),
            cts[i]
        )
    end

    nothing
end


# well this looks nice
function get_positions(rgraph::RGraph, sgraph::SGraph)
    map(1:sgraph.N_nodes) do i # IndexT
        uvw, LID = lattice_indices(i)
        uvw * rgraph.nodes[LID].bravais
    end
end



# TODO Square creation (much later, maybe)
# - use old uc creation
# - decompose that to sc lattices (take each atom as origin)
# - do basis fill with that






#=
# About "index order"
# With that I want to describe when which index changes. [1, 2, 3] would imply
# that the first index changes the fastes, while the third changes the slowest.

# Create points based on their Bravais lattice. This may result in incomplete
# regions for combined lattices
function Basisfill(rgraph::RGraph, N::Integer; border::Symbol=:periodic)
    # Indexing
    # index order [(3, 2, 1), 4]
    LSize = N^3
    flattening_vector = Vec3i(N^2, N, 1)

    # index order [1, (2, 3, 4)]
    # flattening_vector = Vec3i(2, 2*N, 2*N^2)


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

    # function periodic(LID, uvw) # (uvw::Vec3i, LID::Int64)
    #     # index order: [1, (2, 3, 4)]
    #     x = 0
    #     y = 0
    #     z = 0
    #
    #     # This does the border recursion
    #     uvw[1] < 0  && (x = N)
    #     uvw[1] >= N && (x = - N)
    #     uvw[2] < 0  && (y = N)
    #     uvw[2] >= N && (y = - N)
    #     uvw[3] < 0  && (z = N)
    #     uvw[3] >= N && (z = - N)
    #     # print(x, ", ", y, ", ", z, " \t ")
    #
    #     1 + (LID-1) + dot(flattening_vector, uvw + Vec3i(x, y, z)) # IndexT
    # end

    # Should always be in bounds; should fit loop
    function open(uvw, LID) # (uvw::Vec3i, LID::Int64)
        # uvw, LID = args
        # index order: [(3, 2, 1), 4]

        (0 > uvw[1] || uvw[1] >= N) && return 0
        (0 > uvw[2] || uvw[2] >= N) && return 0
        (0 > uvw[3] || uvw[3] >= N) && return 0

        1 + LSize * (LID-1) + dot(flattening_vector, uvw) # IndexT
    end


    # function open(uvw, LID) # (uvw::Vec3i, LID::Int64)
    #     # uvw, LID = args
    #     # index order: [1, (2, 3, 4)]
    #
    #     (0 > uvw[1] || uvw[1] >= N) && return 0
    #     (0 > uvw[2] || uvw[2] >= N) && return 0
    #     (0 > uvw[3] || uvw[3] >= N) && return 0
    #
    #     1 + (LID-1) + dot(flattening_vector, uvw) # IndexT
    # end

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

    # function decompose_index(i::Int64)
    #     # index order [1, (2, 3, 4)]
    #     LID = (i-1) % 2 + 1                             # 1-based
    #     sublattice_index = i - 1                        # 0-based
    #
    #     # flattening_vector includes the 2x from LID
    #     w = div(sublattice_index, flattening_vector[1]) # 0-based
    #     r = sublattice_index % flattening_vector[1]     # 0-based
    #     v = div(r, flattening_vector[2])                # 0-based
    #     u = r % flattening_vector[2]                    # 0-based
    #
    #     Vec3i(u, v, w), LID
    # end


    # Initialisation
    # The order of this arrays is irrelevant. The nodes will be set according
    # to IDs and the indice functions
    # nodes = [
    #     SNode(
    #         Edge1[], Int64[],          # zeros(IndexT, 12),                 # this would fail for oepn borders, so we don't do it :(
    #         [Edge1[] for _ in 1:4]      # [Array(Edge1, 18) for _ in 1:4]   # 18 including x -> this -> ...
    #     ) for __ in 1:N^3 for _ in 1:length(rgraph.nodes)
    # ]
    nodes = [
        SNode(
            Edge1[], Int64[],
            [Edge1[] for _ in 1:4]
        ) for __ in 1:N^3, _ in 1:length(rgraph.nodes)
    ][:]

    # The ordering of this HAS to match the indice function!
    # index order: [(3, 2, 1), 4]

    IDs = [
        (Vec3i(u, v, w), i)
        for i in 1:length(rgraph.nodes)
        for u in 0:N-1 for v in 0:N-1 for w in 0:N-1
    ]       # the last loop runs through first
    # IDs = [
    #     (Vec3i(u, v, w), i)
    #     for w in 0:N-1, v in 0:N-1, u in 0:N-1, i in 1:length(rgraph.nodes)
    # ][:] # checked (explicitly for two nodes, first neighbours)

    # Potentially more efficient because Edge1's are closer together
    # index order: [1, (2, 3, 4)]
    # IDs = [
    #     (i, Vec3i(u, v, w))
    #     for i in 1:length(rgraph.nodes), u in 0:N-1, v in 0:N-1, w in 0:N-1
    # ]


    # populate neighbours
    first = Edge1[] # [SEdge[], SEdge[]]
    second = Edge2[]  # [SEdge[], SEdge[]]
    paths = Vector{Edge1}[]
    connect_nodes!(rgraph, nodes, first, second, paths, IDs, flat_index)


    g = SGraph(
        length(nodes), length(first) + length(second), length(paths),
        nodes, first, second, paths
    )

    # sim, _, lattice_indices
    g, flat_index, decompose_index
end
=#
