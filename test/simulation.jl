r = RGraph(diamond("A"), 2)
generatePaths!(r)
sim, _, _ = Basisfill(r, 3)

@testset "Checking lattice graph parameter" begin
    @test sim.N_nodes == 2*3^3
    @test length(sim.first) == 4*3^3    # 1/2 * 4 * N_nodes
    @test length(sim.second) == 12*3^3
    @test mapreduce(length, +, sim.paths) == 36*3^3

    for n in sim.nodes
        @test length(n.first) == 4
        @test length(n.second) == 12
        @test mapreduce(length, +, n.paths) == 72 # 2*4*3*3
        # 2 - paths like n - x - y - z and x - n - y - z
        # 4, 3 - 4 nearest neighbors; 3 without reversing
    end

end

spins = rand_spin(sim.N_nodes)
ps.init_edges!(sim, spins)
Js = [(-1.0, -0.5), (-1.0, -0.5), (1.0, 0.0), (0.0, -1.0)]
h = Point3(0.)
E_tot = totalEnergy(sim, spins, Js, h)

# high temperature (10), so most spin flip attempts should be successful
for _ in 1:100
    E_tot = sweep(sim, spins, E_tot, Js, 1./10., h)
end
@test E_tot ≈ totalEnergy(sim, spins, Js, h)
