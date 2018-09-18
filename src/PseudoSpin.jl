module PseudoSpin

using StaticArrays.FixedSizeArrays

const Vec3 = Vec{3}
const Point3 = Point{3}
const Point3f0 = Point{3, Float32}
const Vec3f0 = Vec{3, Float32}


# I'm implementing new methods for these functions
import Base: *, ==, in, findfirst, push!#, Base.show
# import Base.start, Base.next, Base.length, Base.done, Base.eltype
import Base: length, next, done, last
import Base: mean, var, range


# Bravais lattice vectors, positions
include("Crystal.jl")
export sc, bcc, fcc, diamond, combine

# Neighbor search and path generation
include("RGraph.jl")
export RGraph, generate_paths!

# Lattice graph
include("SGraph.jl")
export SGraph, Basisfill
export rand_spin, rand_spin!, get_positions

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
