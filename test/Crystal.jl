c = diamond("A", Point3f0(1.0, 1.0, 1.0))

@testset "diamond constructor test" begin
    @test c["A 1"].pos ≈ Poin3f0(1.0, 1.0, 1.0)
    @test c["A 1"].x ≈ Vec3f0(0.5, 0.5, 0.0)
    @test c["A 1"].y ≈ Vec3f0(0.5, 0.0, 0.5)
    @test c["A 1"].z ≈ Vec3f0(0.0, 0.5, 0.5)

    @test c["A 2"].pos ≈ Poin3f0(1.25, 1.25, 1.25)
    @test c["A 2"].x ≈ Vec3f0(0.5, 0.5, 0.0)
    @test c["A 2"].y ≈ Vec3f0(0.5, 0.0, 0.5)
    @test c["A 2"].z ≈ Vec3f0(0.0, 0.5, 0.5)
end

center_around!(c["A 1"])
center_around!(c["A 2"])

@testset "center_around! test" begin
    @test c["A 1"].pos ≈ Poin3f0(0.0, 0.0, 0.0)
    @test c["A 2"].pos ≈ Poin3f0(0.25, 0.25, 0.25)
end
