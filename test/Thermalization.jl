@testset "Checking Thermalization Methods" begin
    # ConstantT
    th = ConstantT(200)
    state = initialize(th, 1.0)
    i = 0
    while !done(th, state)
        beta, state = next(th, state)
        i += 1
        @test beta ≈ 1.0
    end
    @test (i == state[1] == 200)
    @test finalize(th, state) ≈ 1.0


    # Freezer
    f = Freezer(200, 2.5, N_switch=50)
    state = initialize(f, 1.0)
    @test state == (0, 1, 1, 2.5 - 1.0, 1.0, 1.0)

    i = 0
    while !done(f, state)
        beta, state = next(f, state)
        i += 1
        if i <= 50
            @test (1/3.00001 <= beta <= 1.000001)
        else
            @test beta ≈ 1.0
        end
    end
    @test (i == state[1] == 200)
    @test finalize(f, state) ≈ 1.0
end
