using PseudoSpin
using Base.Test

const ps = PseudoSpin

const Vec3f0 = ps.Vec3f0
const Vec3i = ps.Vec3i
const Vec3 = ps.Vec3
const Point3f0 = ps.Point3f0
const Point3 = ps.Point3
include("Crystal.jl")

include("RGraph.jl")

@test all(abs.(reduce(+, rand_spin(1_000_000)) / 1_000_000) .<= Point3(0.01, 0.01, 0.01))

include("simulation.jl")


ARGS = ["parameters/test.param"]
println("Attempting simulation with 20k + 80k sweeps. Runtime: (estimate: ~4min)")
@time include("mainf.jl")


# @test begin
#     # ...
#     true
# end
#
#
# @testset "name" begin
#     @test f(1) = 1
#     @test f(10) = 10
# end
