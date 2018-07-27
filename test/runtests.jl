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


exit(0)
################################################################################

# Bunch of random garbage


begin
    println("----------------------------")
    code_warntype(
        PseudoSpin.measure!,
        (
            PseudoSpin.SGraph,
            Vector{PseudoSpin.Point3{Float64}},
            Float64,
            Parameters,
            IOStream,
            typeof(PseudoSpin.sweep_J1),
            Int64,
            Bool,
            Int64,
            typeof(PseudoSpin.flip1)
            # PseudoSpin.Point3{Float64},
            # Float64,
            # Int64,
        )
    )
    println("----------------------------")
end


begin
    println("----------------------------")
    code_warntype(
        sweep_picker(Parameters(J1 = 1.0)),
        (
            PseudoSpin.SGraph,
            Vector{PseudoSpin.Point3{Float64}},
            Float64,
            Float64,
            Parameters,
            # IOStream,
            # typeof(PseudoSpin.sweep_J1),
            # Int64,
            # Bool,
            # Int64,
            # typeof(PseudoSpin.flip1)
            # PseudoSpin.Point3{Float64},
            # Float64,
            # Int64,
        )
    )
    println("----------------------------")
end


begin
    println("----------------------------")
    code_warntype(
        PseudoSpin.deltaEnergy_J1J2J3Kghzeta,
        (
            PseudoSpin.SNode,
            Vector{PseudoSpin.Point3{Float64}},
            Int64,
            PseudoSpin.Point3{Float64},
            Vector{Float64},
            Vector{Float64},
            Parameters,
            # IOStream,
            # typeof(PseudoSpin.sweep_J1),
            # Bool,
            # Int64,
            # typeof(PseudoSpin.flip1)
            # Float64,
            # Int64,
        )
    )
    println("----------------------------")
end

f(x, y) = x - y^2
begin
    println("----------------------------")
    code_warntype(
        PseudoSpin.jackknife,
        (
            typeof(f),
            Vector{Float64},
            Vector{Float64}
        )
    )
    println("----------------------------")
end

Profile.clear_malloc_data()
@profile simulate!(
    J1s = (-1.0, -0.4),
    J2s = (0.3, 0.12),
    K = -0.1,
    g = 0.52,
    h = PseudoSpin.Point3{Float64}(0.01, 0., 0.5),
    T = 0.624,
    TH_sweeps = 10_000,
    ME_sweeps = 10_000,
    L = 6,
    path = "output/",
    filename = "full_test"
)
rm("output/full_test.part")
Profile.print()
using ProfileView
ProfileView.view()



--J1              0.003635 seconds (17.94 k allocations: 1.962 MiB)
--J2              0.003721 seconds (34.14 k allocations: 2.457 MiB)
--K               0.004088 seconds (17.94 k allocations: 1.962 MiB)
--h               0.003570 seconds (34.14 k allocations: 2.457 MiB)
--J1,J2,K,h       0.004241 seconds (34.14 k allocations: 2.456 MiB)
--J1,J2,K,h,g     0.004660 seconds (34.14 k allocations: 2.456 MiB)
--J1,g            0.003892 seconds (34.14 k allocations: 2.456 MiB)
--zeta            0.004015 seconds (109.74 k allocations: 3.610 MiB)
--J3              0.004002 seconds (34.14 k allocations: 2.456 MiB)



--J1              0.003160 seconds (17.94 k allocations: 1.962 MiB)
--J2              0.003656 seconds (17.94 k allocations: 1.962 MiB)
--K               0.003797 seconds (17.94 k allocations: 1.962 MiB)
--h               0.003099 seconds (17.94 k allocations: 1.962 MiB)
--J1,J2,K,h       0.003837 seconds (17.94 k allocations: 1.962 MiB)
--J1,J2,K,h,g     0.004277 seconds (17.94 k allocations: 1.962 MiB)
--J1,g            0.003618 seconds (17.94 k allocations: 1.962 MiB)
--zeta            0.004116 seconds (93.54 k allocations: 3.116 MiB)
--J3              0.003409 seconds (17.94 k allocations: 1.962 MiB)
