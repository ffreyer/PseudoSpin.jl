#=
################################################################################
### Types for Recursive Graphs
################################################################################

A recursive Graph is a minimal Graph describing a Crystal cell through nodes
(atoms) and their (first to nth) nearest neighbours.

For all equivalent atoms in a cubic cell, one Node is created. Two atoms are the
same if their neighbours are at the same relative positions and are of the same
type of atom.

Each Edge describes the direct path to a neighbour (of some order) through a
vector. The edge may point to the same node it originates from. The distance of
and edge is given in Node, as multiple edges of same length a grouped there.
=#

abstract type AbstractREdge end

mutable struct SimpleREdge <: AbstractREdge
    dir::Vec3i
    from::Integer # AbstractNode
    to::Integer   # AbstractNode
end

mutable struct PathREdge <: AbstractREdge
    dirs::Vector{Vec3i} # should contain 3 Vec3i's
    from::Integer
    tos::Vector{Integer} # should contain 3 indices
end

mutable struct RNode #<: AbstractNode
    atom::String
    bravais::Bravais
    edges::Vector{Vector{AbstractREdge}} #Abstract
    distances::Vector{AbstractFloat}
end

mutable struct RGraph
    nodes::Vector{RNode}
    edges::Vector{AbstractREdge}
end


################################################################################
### Constructor for Crystal -> RGraph
################################################################################

"""
    RGraph(c::Crystal, N::Integer)

Creates a minimal Graph describing a cubic cell with N levels of neighbours.

Each Node in the resulting Graph represents a Bravais lattice. Each Edge
describes the direction to the next neighbour (up to order N) as well as the
node at that point.
"""
function RGraph(c::Crystal, N::Integer)

    function dist(uvw::Vec3i, B_::Bravais, B::Bravais)
        v = uvw * B_ - B.pos
        dot(v, v)
    end

    nodes = [
        RNode(
            key,
            c[key],
            [AbstractREdge[] for _ in 1:N], Float64[]
        ) for key in keys(c)
    ]
    edges = SimpleREdge[]

    Bs = collect(values(c))
    N_Atoms = length(Bs)

    # For every Bravais lattice:
    #   - center every Bravais lattice next to it
    #   - gather directions and distances
    #   - get order of distances
    #   - create edges up to order N
    for i in 1:N_Atoms
        dirs = Vec3i[]
        sq_distances = Float64[]
        B = Bs[i]

        for j in 1:N_Atoms
            B_ = Bs[j]

            # If a Bravais lattice that already has edges happens to be moved,
            # update the dir
            delta_uvw_j = center_around!(B_, B.pos)
            if j > i
                for el in nodes[j].edges
                    for e in el
                        e.dir += delta_uvw_j
                    end
                end
            end

            _dirs = [Vec3i(x, y, z) for x in -N:N for y in -N:N for z in -N:N]
            append!(dirs, _dirs)
            append!(sq_distances, map(uvw -> dist(uvw, B_, B), _dirs))
        end

        s_ind = sortperm(sq_distances)
        order = 0
        j = 1

        while true
            if order != 0
                e = SimpleREdge(
                    dirs[s_ind[j]],
                    i,
                    div(s_ind[j]-1, (2*N +1)^3) + 1 # sq_distances contains atoms
                    # nodes[i],                       [1, 1, 1, 1, ..., 2, 2, 2, ...]
                    # nodes[div(s_ind[j]-1, (2*N +1)^3) + 1]
                )
                push!(nodes[i].edges[order], e)
                push!(edges, e)
            end

            j += 1

            if abs(sq_distances[s_ind[j-1]] - sq_distances[s_ind[j]]) > 1e-10
                order += 1
                if order <= N
                    push!(nodes[i].distances, sqrt(sq_distances[s_ind[j]]))
                else
                    break
                end
            end
        end
    end

    RGraph(nodes, edges)
end


"""
    generate_paths!(RGraph)

Adds edges for the 4-spin term to an existing RGraph.
"""
function generate_paths!(rgraph::RGraph)
    for i1 in eachindex(rgraph.nodes)       # for each node n1
        paths = PathREdge[]                 # n1 is at Vec3i(0, 0, 0) per definition
        for e1 in rgraph.nodes[i1].edges[1]     # each edge n1 -> n2
            i2 = e1.to                          # n2 is at dir1
            dir1 = e1.dir
            for e2 in rgraph.nodes[i2].edges[1]     # each edge n2 -> n3 with n3 != n1
                i3 = e2.to                          # n3 is at dir1 + dir2 (because the unit vectors are the same)
                dir2 = dir1 + e2.dir                # may not be at Vec3i(0, 0, 0) (no reversal)
                # no reversal
                if dir2 == Vec3i(0, 0, 0); continue end
                for e3 in rgraph.nodes[i3].edges[1]     # each edge n3 -> n4 with n4 != n3
                    i4 = e3.to                          # n4 is at dir2 + dir3
                    dir3 = dir2 + e3.dir                # may not be at dir1
                    # no reversal
                    if dir3 == dir1; continue end

                    # create and push path
                    path = PathREdge(
                        [dir1, dir2, dir3],
                        i1,
                        [i2, i3, i4]
                    )

                    push!(paths, path)
                end
            end
        end
        push!(rgraph.nodes[i1].edges, paths)
        append!(rgraph.edges, paths)
    end

    nothing
end
