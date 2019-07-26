module PseudoSpin

using StaticArrays
using SphereSurfaceHistogram


using LinearAlgebra, Printf, Distributed
# I'm implementing new methods for these functions
import Base: *, ==, in, findfirst, push!#, Base.show
# import Base.start, Base.next, Base.length, Base.done, Base.eltype
import Base: length, last
if VERSION >= v"0.7.0"
    import Statistics: mean, var
else
    import Base: mean, var, next, done
end

# Bravais lattice vectors, positions
include("Crystal.jl")
export sc, bcc, fcc, diamond, combine

# Neighbor search and path generation
include("RGraph.jl")
export RGraph, generate_paths!

# Lattice graph
include("SGraph.jl")
export SGraph, Basisfill, get_positions

include("Updates.jl")
export rand_spin, rand_spin!
export rand_XY_spin, rand_XY_spin!
export rand_red_XY_spin, rand_XY_rot_matrix
export rand_3fold_XY_rotation, yaxis_mirror, flipflop_rot

# Essential Metropolis functions (energy, sweep, ...)
include("Model.jl")
export init_edges!
export totalEnergy#, deltaEnergy
# export kernel, sweep
export Parameters, sweep_picker

# Simulated Annealing, parallel tempering
include("ThermalizationMethods.jl")
export parallel_tempering
export Freezer, ConstantT
export NoParallelTempering, ParallelTempering
export ProbabilityEqualizer
export initialize, next, done, last, length, T_max, current_index

# Binning Analysis, Jackknife, histograms
include("DataAnalysis.jl")
export BinnerA, mean, var, tau
export BinnerH, jackknife

# File writing (duh)
include("FileWriter.jl")
export write_header!, write_BA!, write_JK!, write_HB!, write_SC!

# full Measurement (no thermalization)
include("Measure.jl")
export measure!, measure_no_paths!

# full simulation (parameters, thermalization)
include("Simulation.jl")
export thermalize!, thermalize_no_paths!, simulate!

end # module
