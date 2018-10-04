using PseudoSpin
if VERSION >= v"0.7.0"
    using Test
    const gc = GC.gc
else
    using Base.Test
end

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
println("Attempting simulation with 10k + 10k sweeps. Runtime: (estimate: ~8s)")
@time include("mainf.jl")
rm("output/full_test.part")

# addprocs(3)
# ARGS = ["parameters/test_mpi.param"]
# println("Attempting simulation with 10k + 10k sweeps. Runtime: (estimate: ~40s)")
# @time include("mainf.jl")
# rm("output/mpi_test1.part")
# rmprocs(workers())

exit(0)
################################################################################

# Bunch of random garbage
addprocs(3)
nprocs()
@everywhere using PseudoSpin
cd("test/")
@profile simulate!(
    path = "output/",
    filename = "mpi_test",
    L = 4,
    J1s = (-0.4, 0.0),
    J2s = (0.0, 0.0),
    K = 0.0,
    T = 0.7,
    # Ts = [0.01, 0.1, 0.7, 0.8],
    TH_sweeps = 10000,
    Freeze_temperature = 2.0,
    ME_sweeps = 10000,
    h = PseudoSpin.Point3{Float64}(0.010000, 0.0, 0.5),
    # thermalizer_method = ParallelTempering,
)
Juno.profiler()
rm("output/mpi_test1.part")
rmprocs(workers())


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




--J1              0.006557 seconds (3.00 k allocations: 140.625 KiB)
--J2              0.007380 seconds (3.00 k allocations: 140.625 KiB)
--K               0.011251 seconds (3.00 k allocations: 140.625 KiB)
--h               0.005898 seconds (3.00 k allocations: 140.625 KiB)
--J1,J2,K,h       0.012720 seconds (3.00 k allocations: 140.625 KiB)
--J1,J2,K,h,g     0.015402 seconds (3.00 k allocations: 140.625 KiB)
--J1,g            0.009772 seconds (3.00 k allocations: 140.625 KiB)
--zeta            0.013739 seconds (759.00 k allocations: 11.673 MiB)
--J3              0.007747 seconds (3.00 k allocations: 140.625 KiB)

Attempting simulation with 10k + 10k sweeps. Runtime: (estimate: ~8s)
  5.918624 seconds (2.02 M allocations: 535.722 MiB, 0.64% gc time)
