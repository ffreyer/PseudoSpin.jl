using LinearAlgebra

@testset "Updates" begin
    # Check if rand_spin(N) creates the correct distribution. Since this is random
    # it may rarely fail.
    @test all(
        abs.(
            reduce(
                +,
                rand_spin(1_000_000)
            ) / 1_000_000
        ) .<= SVector{3, Float64}(0.01, 0.01, 0.01)
    )
    @test all(
        abs.(
            reduce(
                +,
                rand_XY_spin(1_000_000)
            ) / 1_000_000
        ) .<= SVector{3, Float64}(0.01, 0.01, 0.0)
    )

    # Check that M is conserved
    spins = rand_XY_spin(1024)
    U = self_balancing_update(spins)
    for _ in 1:100_000
        idxs, new_spins = PseudoSpin.apply(U, spins)
        spins[idxs] = new_spins
    end
    @test normalize(sum(spins)) ≈ U.eM

    spins = rand_XY_spin(1024)
    U = self_balancing_update2(spins)
    for _ in 1:100_000
        idxs, new_spins = PseudoSpin.apply(U, spins)
        spins[idxs] = new_spins
    end
    @test normalize(sum(spins)) ≈ U.eM

end
