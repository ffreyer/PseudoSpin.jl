using PseudoSpin
if VERSION >= v"0.7.0"
    using Test, Distributed
    const gc = GC.gc
else
    using Base.Test
end

isfile("output/full_test.part") && rm("output/full_test.part")

const ps = PseudoSpin
const SVector = ps.SVector

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


include("Thermalization.jl")

include("simulation.jl")


ARGS = ["parameters/test.param"]
println("Attempting simulation with 10k + 10k sweeps. Runtime: (estimate: ~6s)")
@time include("mainf.jl")
rm("output/full_test.part")

println("Attempting simulation with global updates, reduced XY-plane spins, 10k + 10k sweeps. Runtime: (estimate: ~6s)")
spins = rand_red_XY_spin(2*6^3)
@test all(S -> S[2] < 0.2, spins)
@time simulate!(
    path = "output/",
    filename = "full_test",
    L = 6,
    sampler = rand_red_XY_spin,
    J1 = 2rand()-1,
    J2 = 2rand()-1,
    J3 = 2rand()-1,
    lambda = 2rand()-1,
    K = 2rand()-1,
    g = 2rand()-1,
    zeta = 2rand()-1,
    h = SVector((2rand(3) .- 1)...),
    T = rand() + 0.1,
    TH_sweeps = 10_000,
    ME_sweeps = 10_000,
    do_global_updates = true,
    global_update = yaxis_mirror
)
rm("output/full_test.part")
@test all(S -> S[2] < 0.2, spins)


ARGS = ["parameters/test_mpi.param"]
println("Attempting parallel simulation with 10k + 10k sweeps. Runtime: (estimate: ~17s)")
@time include("main_cluster.jl")
rm("output/mpi_test1.part")
rm("output/mpi_test2.part")
rm("output/mpi_test3.part")
rm("output/mpi_test4.part")


exit(0)


using Revise
using PseudoSpin
using PseudoSpinUtils
# using SphereSurfaceHistogram

PseudoSpin.SSHBinner
cell = @datawith J1 == -1 hy == 0 hz == 0
data = read_merged(cell[:path][1])

@time simulate!(
    path = "/home/frederic/.julia/dev/PseudoSpin/test/",
    L = data[:L],
    spins = map(PseudoSpin.Point3{Float64}, data[:spins][50]),
    J1 = data[:J1],
    J2 = data[:J2],
    J3 = data[:J3],
    lambda = data[:lambda],
    K = data[:K],
    g = data[:gamma],
    zeta = data[:zeta],
    T = data[:T][50],
    TH_sweeps = 10_000,
    ME_sweeps = 10_000
)

L = 6
@time simulate!(
    path = "/files/home/part2/ffreyer/Xenial/.julia/dev/PseudoSpin/test/output/",
    J1 = -1.0,
    L = L,
    # spins = PseudoSpin.SVector{3, Float64}[[1, 0, 0] for _ in 1:2L^3],
    T = 0.5,
    TH_sweeps = 10_000,
    ME_sweeps = 100_000
)


function areas(B)
    output = Float64[]
    for i in eachindex(B.phi_divisions)
        t1 = B.thetas[i+1]
        t2 = B.thetas[i]
        for _ in 1:B.phi_divisions[i]
            dphi = 2pi / B.phi_divisions[i]
            push!(output, dphi * (cos(t2) - cos(t1)))
        end
    end
    output
end
dAs = areas(B)
mean(dAs)
mean(dAs) |> eps
std(dAs)
maximum(dAs) - minimum(dAs)


################################################################################

# Bunch of random garbage
addprocs(3)
nprocs()
cd("test/")
@everywhere using PseudoSpin, Profile
cd("/home/frederic/.julia/v0.7/PseudoSpin/test/")
@profile simulate!(
    path = "output/",
    filename = "mpi_test",
    L = 4,
    J1s = (-0.4, 0.1),
    J2s = (1.0, 0.5),
    J3s = (0.3, 1.0),
    g = 0.44,
    zeta = 0.7,
    K = 0.3,
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


@code_warntype simulate!(
    path = "output/",
    filename = "mpi_test",
    L = 4,
    J1s = (-0.4, 0.1),
    J2s = (1.0, 0.5),
    J3s = (0.3, 1.0),
    g = 0.44,
    zeta = 0.7,
    K = 0.3,
    T = 0.7,
    # Ts = [0.01, 0.1, 0.7, 0.8],
    TH_sweeps = 10000,
    Freeze_temperature = 2.0,
    ME_sweeps = 10000,
    h = PseudoSpin.Point3{Float64}(0.010000, 0.0, 0.5),
    # thermalizer_method = ParallelTempering,
)


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

begin
    println("----------------------------")
    code_warntype(
        PseudoSpin._parallel_tempering!,
        (
            SGraph,
            Vector{PseudoSpin.Point3{Float64}},
            Float64,
            Float64
        )
    )
    println("----------------------------")
end




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

r = RGraph(diamond("A"), 3)
generate_paths!(r)
sim, _, __ = Basisfill(r, 8)
Js = [(-0.2, 0.4), (1.0, -0.3), (-0.25, 0.0), (0.0, 0.9)]
h = PseudoSpin.Point3(-rand(3))
param = Parameters(J1s = Js[1], J2s = Js[2], K = Js[3][1]*Js[4][2], h = h, zeta=0.3, g=0.7)

spins = rand_spin(1024)
E_tot = totalEnergy(sim, spins, param)
using BenchmarkTools
sweep = sweep_picker(param)
test(sweep, sim, spins, E_tot, param) = sweep(sim, spins, E_tot, 1.0/10., param)

@benchmark test($sweep, $sim, $spins, $E_tot, $param)
@benchmark sweep(sim, spins, E_tot, 1.0/10., param)
@benchmark append!($B, $spins)


f(x) = a - bx - cx^3
f(1) = 0 => a - b - c = 0
f(-1) = 1 => a + b + c = 1

g(x, p) = @. (p[1] + p[2] - p[1]*x - p[2]*x^3) / (2p[1] + 2p[2])
param = [5., 20.]
ys = eachindex(B.zs) / length(B.zs)
fit = curve_fit(g, B.zs, ys, param)
f(x) = 0.5 - 7/22 * x - 7/132 * x^3
f(x) = 0.5 - 7/22 * x - 1/22 * x^3 - 1/22 * x^5- 1/22 * x^15 - 1/22 * x^33
f(x) = 0.5 - 7/22 * x - 1x^5 / 22 - 3x^11 / 22
f(x) = 0.5 - 7/22 * x - 1x^5 / 22 - 1x^7 / 22 - 1x^21 / 22 - 1x^31 / 22
f(x) = 0.5 - 7/22 * x - 1x^5 / 22 - 1*x^9 / 22 - 1*x^11 / 22 - x^31 / 22
f(x) = 0.5 - 28/88 * x - x^3 / 88 - 3*x^5 / 88 - x^9 / 22 - 1*x^11 / 22 - x^31 / 22
f(x) = 0.5 - 1/88 * (28x + x^3 + 5x^5 + 4x^9 + 0*x^17 + 6x^31)

const c = 1/22
f(x) = 0.5 - c*x*(
    7.0 + (((x*x)*x)*x)*(
        1.0 + (((x*x)*x)*x) * (
            1.0 + (((((((x*x)*x)*x)*x)*x)*x)*x)*(
                1.0 + (((((((((((((x*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)
            )
        )
    )
)
f(x) = 0.5 - 7/22 * x - 1x^5 / 22 - 1*x^9 / 22 - 1*x^17 / 22 - x^31 / 22

g(x, l) = trunc(Int64,
    l * (0.5 - c*x*(
        7.0 + (((x*x)*x)*x)*(
            1.0 + (((x*x)*x)*x) * (
                1.0 + (((((((x*x)*x)*x)*x)*x)*x)*x)*(
                    1.0 + (((((((((((((x*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)
                )
            )
        )
    ))
)
g2(x, l) = trunc(Int64, l*acos(x))
