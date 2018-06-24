@testset "Checking Thermalization Methods" begin
    # ConstantT
    th = ConstantT(N = 200)
    @test length(th) == 200
    beta, state = initialize(th, 1.0)
    i = 0
    while !done(th, state)
        beta, state = next(th, state)
        i += 1
        @test beta ≈ 1.0
    end
    @test (i == current_index(th, state) == 200)
    @test last(th, state) ≈ 1.0


    # Freezer
    f = Freezer(N = 200, T_max = 2.5, N_switch = 50)
    beta, state = initialize(f, 1.0)
    @test length(f) == 200
    @test state == (0, 1, 1, 2.5 - 1.0, 1.0, 1.0) # meh

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
    @test (i == current_index(f, state) == 200)
    @test last(f, state) ≈ 1.0
    @test T_max(f) == 2.5
end
