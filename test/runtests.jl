using PseudoSpin
using Base.Test

isfile("output/full_test.part") && rm("output/full_test.part")

const ps = PseudoSpin

const Vec3f0 = ps.Vec3f0
const Vec3i = ps.Vec3i
const Vec3 = ps.Vec3
const Point3f0 = ps.Point3f0
const Point3 = ps.Point3
include("Crystal.jl")

include("RGraph.jl")

# Check if rand_spin(N) creates the correct distribution. Since this is random
# it may rarely fail.
@test all(
    abs.(
        reduce(
            +,
            rand_spin(1_000_000)
        ) / 1_000_000
    ) .<= Point3(0.01, 0.01, 0.01)
)

include("Thermalization.jl")

include("simulation.jl")


ARGS = ["parameters/test.param"]
println("Attempting simulation with 10k + 10k sweeps. Runtime: (estimate: ~40s)")
@time include("mainf.jl")
rm("output/full_test.part")
