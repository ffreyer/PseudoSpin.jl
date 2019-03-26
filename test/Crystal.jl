c = diamond("A", SVector{3, Float32}(1.0, 1.0, 1.0))

@testset "diamond constructor test" begin
    @test c["A 1"].pos ≈ SVector{3, Float32}(1.0, 1.0, 1.0)
    @test c["A 1"].x ≈ SVector{3, Float32}(0.5, 0.5, 0.0)
    @test c["A 1"].y ≈ SVector{3, Float32}(0.5, 0.0, 0.5)
    @test c["A 1"].z ≈ SVector{3, Float32}(0.0, 0.5, 0.5)

    @test c["A 2"].pos ≈ SVector{3, Float32}(1.25, 1.25, 1.25)
    @test c["A 2"].x ≈ SVector{3, Float32}(0.5, 0.5, 0.0)
    @test c["A 2"].y ≈ SVector{3, Float32}(0.5, 0.0, 0.5)
    @test c["A 2"].z ≈ SVector{3, Float32}(0.0, 0.5, 0.5)
end

ps.center_around!(c["A 1"])
ps.center_around!(c["A 2"])

@testset "center_around! test" begin
    @test c["A 1"].pos ≈ SVector{3, Float32}(0.0, 0.0, 0.0)
    @test c["A 2"].pos ≈ SVector{3, Float32}(0.25, 0.25, 0.25)
end
