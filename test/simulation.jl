r = RGraph(diamond("A"), 2)
sim, _, _ = Basisfill(r, 3)

@testset "Checking lattice graph parameter" begin
    @test sim.N_nodes == 2*3^3
    @test length(sim.first) == 4*3^3    # 1/2 * 4 * N_nodes
    @test length(sim.second) == 12*3^3
    @test mapreduce(length, +, sim.paths) == 36*3^3

    for n in sim.nodes
        @test length(n.first) == 4
        @test length(n.second) == 12
        @test mapreduce(length, +, n.paths) == 36
    end

end
