r = RGraph(diamond("A"), 3)
generate_paths!(r)
sim, _, _ = Basisfill(r, 3)

@testset "Checking lattice graph parameter" begin
    @test sim.N_nodes == 2*3^3
    @test length(sim.first) == 4*3^3    # 1/2 * 4 * N_nodes
    @test length(sim.second) == 12*3^3
    @test length(sim.third) == 12*3^3
    @test length(sim.paths) == 4*3^3
    @test mapreduce(length, +, sim.paths) == 36*3^3

    for n in sim.nodes
        @test length(n.first) == 4
        @test length(n.second) == 12
        @test length(n.third) == 12
        @test mapreduce(length, +, n.paths) == 72 # 2*4*3*3
        # 2 - paths like n - x - y - z and x - n - y - z
        # 4, 3 - 4 nearest neighbors; 3 without reversing
    end

end

@testset "Checking Binning Analysis" begin
    N = 1024
    values = rand(N)
    BA = BinnerA(10)
    for v in values
        push!(BA, v)
    end
    @test BA.count[1] == N
    @test mean(BA) ≈ mean(values)
    # @test isapprox(var(BA, 1), var(values) / N, rtol=0.05abs(var(values)))
    @test var(BA) != var(values) / N
end

# Old Test
# @testset "Checking Freezer" begin
#     f1 = Freezer(200, 2.5, N_switch=50)
#     for beta in cool_to(f1, 1.0)
#         @test (1/2.500001 <= beta <= 1.000001)
#     end
# end


spins = rand_spin(sim.N_nodes)
ps.init_edges!(sim, spins)
Js = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
h = Point3(0.)
g = 0.0

function test_param!(param::Parameters)
    @test param.J1 == Js[1]
    @test param.J2 == Js[2]
    @test param.K == Js[3][1] * Js[4][2]
    @test param.g == g
    @test param.h == h
end

# high temperature (10), so most spin flip attempts should be successful
@testset "Testing E_total, sweep stack, parameters" begin
    #-----------------------------------------------------
    ### anisotropic
    # NN
    Js = [(1.0, 0.5), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
    param = Parameters(J1s = Js[1])
    test_param!(param)

    old_E_tot = totalEnergy(sim, spins, Js, h)
    E_tot = totalEnergy(sim, spins, param)
    @test E_tot ≈ old_E_tot
    sweep = sweep_picker(param)
    @test sweep == PseudoSpin.sweep_J1

    E_tot = sweep(sim, spins, E_tot, 1./10., param)
    print("--J1           	"); gc()
    @time for _ in 1:100
        E_tot = sweep(sim, spins, E_tot, 1./10., param)
    end
    @test E_tot ≈ totalEnergy(sim, spins, Js, h)
    @test E_tot ≈ totalEnergy(sim, spins, param)

    #-----------------------------------------------------

    # NNN
    Js = [(0.0, 0.0), (-1.2, 0.8), (0.0, 0.0), (0.0, 0.0)]
    param = Parameters(J2s = Js[2])
    test_param!(param)

    old_E_tot = totalEnergy(sim, spins, Js, h)
    E_tot = totalEnergy(sim, spins, param)
    @test E_tot ≈ old_E_tot
    sweep = sweep_picker(param)
    @test sweep == PseudoSpin.sweep_J2

    E_tot = sweep(sim, spins, E_tot, 1./10., param)
    print("--J2           	"); gc()
    @time for _ in 1:100
        E_tot = sweep(sim, spins, E_tot, 1./10., param)
    end
    @test E_tot ≈ totalEnergy(sim, spins, Js, h)
    @test E_tot ≈ totalEnergy(sim, spins, param)

    #-----------------------------------------------------

    # paths
    Js = [(0.0, 0.0), (0.0, 0.0), (1.5, 0.0), (0.0, 1.0)]
    param = Parameters(K = Js[3][1])
    test_param!(param)

    old_E_tot = totalEnergy(sim, spins, Js, h)
    E_tot = totalEnergy(sim, spins, param)
    @test E_tot ≈ old_E_tot
    sweep = sweep_picker(param)
    @test sweep == PseudoSpin.sweep_K

    E_tot = sweep(sim, spins, E_tot, 1./10., param)
    print("--K            	"); gc()
    @time for _ in 1:100
        E_tot = sweep(sim, spins, E_tot, 1./10., param)
    end
    @test E_tot ≈ totalEnergy(sim, spins, Js, h)
    @test E_tot ≈ totalEnergy(sim, spins, param)

    #-----------------------------------------------------

    # fields
    Js = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
    h = Point3(0.3, 0.6, 0.9)
    param = Parameters(h = h)
    test_param!(param)

    old_E_tot = totalEnergy(sim, spins, Js, h)
    E_tot = totalEnergy(sim, spins, param)
    @test E_tot ≈ old_E_tot
    sweep = sweep_picker(param)
    @test sweep == PseudoSpin.sweep_h

    E_tot = sweep(sim, spins, E_tot, 1./10., param)
    print("--h            	"); gc()
    @time for _ in 1:100
        E_tot = sweep(sim, spins, E_tot, 1./10., param)
    end
    @test E_tot ≈ totalEnergy(sim, spins, Js, h)
    @test E_tot ≈ totalEnergy(sim, spins, param)

    #-----------------------------------------------------

    # all
    Js = [(-0.2, 0.4), (1.0, -0.3), (-0.25, 0.0), (0.0, 0.9)]
    h = Point3(-rand(3))
    param = Parameters(J1s = Js[1], J2s = Js[2], K = Js[3][1]*Js[4][2], h = h)
    test_param!(param)

    old_E_tot = totalEnergy(sim, spins, Js, h)
    E_tot = totalEnergy(sim, spins, param)
    @test E_tot ≈ old_E_tot
    sweep = sweep_picker(param)
    @test sweep == PseudoSpin.sweep_J1J2Kh

    E_tot = sweep(sim, spins, E_tot, 1./10., param)
    print("--J1,J2,K,h    	"); gc()
    @time for _ in 1:100
        E_tot = sweep(sim, spins, E_tot, 1./10., param)
    end
    @test E_tot ≈ totalEnergy(sim, spins, Js, h)
    @test E_tot ≈ totalEnergy(sim, spins, param)

    #-----------------------------------------------------

    ### with g
    g = 0.17
    param = Parameters(J1s=Js[1], J2s=Js[2], K=Js[3][1]*Js[4][2], h=h, g=g)
    test_param!(param)

    old_E_tot = totalEnergy(sim, spins, Js, h, g)
    E_tot = totalEnergy(sim, spins, param)
    @test E_tot ≈ old_E_tot
    sweep = sweep_picker(param)
    @test sweep == PseudoSpin.sweep_J1J2Kgh

    E_tot = sweep(sim, spins, E_tot, 1./10., param)
    print("--J1,J2,K,h,g  	"); gc()
    @time for _ in 1:100
        E_tot = sweep(sim, spins, E_tot, 1./10., param)
    end
    @test E_tot ≈ totalEnergy(sim, spins, Js, h, g)
    @test E_tot ≈ totalEnergy(sim, spins, param)

    #-----------------------------------------------------

    g = 0.85
    Js = [(1.0, 0.5), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
    h = Point3(0.)
    param = Parameters(J1s=Js[1], g=g)
    test_param!(param)

    old_E_tot = totalEnergy(sim, spins, Js, h, g)
    E_tot = totalEnergy(sim, spins, param)
    @test E_tot ≈ old_E_tot
    sweep = sweep_picker(param)
    @test sweep == PseudoSpin.sweep_J1g

    E_tot = sweep(sim, spins, E_tot, 1./10., param)
    print("--J1,g         	"); gc()
    @time for _ in 1:100
        E_tot = sweep(sim, spins, E_tot, 1./10., param)
    end
    @test E_tot ≈ totalEnergy(sim, spins, Js, h, g)
    @test E_tot ≈ totalEnergy(sim, spins, param)

    #-----------------------------------------------------

    param = Parameters(zeta = 1.0)
    E_tot = totalEnergy(sim, spins, param)
    sweep = sweep_picker(param)
    @test sweep == PseudoSpin.sweep_zeta

    E_tot = sweep(sim, spins, E_tot, 1./10., param)
    print("--zeta         	"); gc()
    @time for _ in 1:100
        E_tot = sweep(sim, spins, E_tot, 1./10., param)
    end
    @test E_tot ≈ totalEnergy(sim, spins, param)

    #-----------------------------------------------------

    param = Parameters(J3 = 1.0)
    E_tot = totalEnergy(sim, spins, param)
    sweep = sweep_picker(param)
    @test sweep == PseudoSpin.sweep_J3

    E_tot = sweep(sim, spins, E_tot, 1./10., param)
    print("--J3           	"); gc()
    @time for _ in 1:100
        E_tot = sweep(sim, spins, E_tot, 1./10., param)
    end
    @test E_tot ≈ totalEnergy(sim, spins, param)
end
