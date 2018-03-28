r = RGraph(diamond("A"), 2)


@testset "RGraph creation (Diamond)" begin
    @test length(r.nodes) == 2
    @test length(r.edges) == 32

    edges = r.nodes[1].edges
    for e in edges[1]
        # nearest neighbors belong to different Bravais lattices
        @test (e.from == 1) && (e.to == 2)
        @test e in r.edges
    end
    for e in edges[2]
        # next nearest neighbors belong to the same Bravais lattice
        @test (e.from == 1) && (e.to == 1)
        @test e in r.edges
    end

    @testset "Node 1" begin
        # dir gives the (u, v, w) of the neighboring node
        # the positions and directions have been checked in "Crystals.jl"
        # check (u, v, w) here
        @test edges[1][1].dir == Vec3i(-1, 0, 0)
        @test edges[1][2].dir == Vec3i(0, -1, 0)
        @test edges[1][3].dir == Vec3i(0, 0, -1)
        @test edges[1][4].dir == Vec3i(0, 0, 0)

        @test edges[2][1].dir == Vec3i(-1, 0, 0)
        @test edges[2][2].dir == Vec3i(-1, 0, 1)
        @test edges[2][3].dir == Vec3i(-1, 1, 0)
        @test edges[2][4].dir == Vec3i(0, -1, 0)
        @test edges[2][5].dir == Vec3i(0, -1, 1)
        @test edges[2][6].dir == Vec3i(0, 0, -1)

        @test edges[2][7].dir == Vec3i(0, 0, 1)
        @test edges[2][8].dir == Vec3i(0, 1, -1)
        @test edges[2][9].dir == Vec3i(0, 1, 0)
        @test edges[2][10].dir == Vec3i(1, -1, 0)
        @test edges[2][11].dir == Vec3i(1, 0, -1)
        @test edges[2][12].dir == Vec3i(1, 0, 0)
    end



    edges = r.nodes[2].edges
    for e in edges[1]
        @test (e.from == 2) && (e.to == 1)
        @test e in r.edges
    end
    for e in edges[2]
        @test (e.from == 2) && (e.to == 2)
        @test e in r.edges
    end

    @testset "Node 2" begin
        @test edges[1][1].dir == Vec3i(0, 0, 0)
        @test edges[1][2].dir == Vec3i(0, 0, 1)
        @test edges[1][3].dir == Vec3i(0, 1, 0)
        @test edges[1][4].dir == Vec3i(1, 0, 0)

        @test edges[2][1].dir == Vec3i(-1, 0, 0)
        @test edges[2][2].dir == Vec3i(-1, 0, 1)
        @test edges[2][3].dir == Vec3i(-1, 1, 0)
        @test edges[2][4].dir == Vec3i(0, -1, 0)
        @test edges[2][5].dir == Vec3i(0, -1, 1)
        @test edges[2][6].dir == Vec3i(0, 0, -1)

        @test edges[2][7].dir == Vec3i(0, 0, 1)
        @test edges[2][8].dir == Vec3i(0, 1, -1)
        @test edges[2][9].dir == Vec3i(0, 1, 0)
        @test edges[2][10].dir == Vec3i(1, -1, 0)
        @test edges[2][11].dir == Vec3i(1, 0, -1)
        @test edges[2][12].dir == Vec3i(1, 0, 0)
    end
end

# Might as well test it in SGraph
# generatePaths!(r)
#
# @testset "generatePaths! test" begin
#
# end
