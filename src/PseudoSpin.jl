module PseudoSpin

__precompile__(true)


# Bleh
# if VERSION == v"0.4.5"
#     typealias _String ASCIIString
# else
#     typealias _String String
# end


# Small FixedSizeArrays, e.g. 3-component vectors, outperform Julia Arrays.
try
    using FixedSizeArrays
catch e
    # println("FixedSizeArrays not available. Trying StaticArrays.")
    try
        using StaticArrays.FixedSizeArrays
    catch e
        println("StaticArrays not available!")
        throw(e)
    end
end
const Vec3 = Vec{3}
const Point3 = Point{3}
const Point3f0 = Point{3, Float32}
const Vec3f0 = Vec{3, Float32}


# I'm implementing new methods for these functions
import Base.*, Base.==, Base.in, Base.findfirst#, Base.show
import Base.start, Base.next, Base.length, Base.done, Base.eltype
import Base.push!, Base.var


# Files used in simulation
include("Crystal.jl")
export sc, bcc, fcc, diamond, combine


include("RGraph.jl")
export RGraph, generate_paths!


include("SGraph.jl")
export SGraph, Basisfill
export rand_spin, rand_spin!, get_positions


include("model.jl")
export init_edges!
export totalEnergy, deltaEnergy
export kernel, sweep


include("simulation.jl")    # TODO rename this
include("Measure.jl")       # TODO name this simulation.jl... probably
export Freezer, ConstantT, cool_to
export BinnerA, var, tau
export BinnerH, jackknife
export measure!, simulate!


#=
    Cubic.jl
Overly general methods to create cubic unit cells.

    RGraph.jl
Defines a cyclic (or recursive) Graph, which describes the neighbors/paths used
in the simulation. The graph only has one node per Bravais lattice (from Cubic)
and one edge per neighbor-connection (i.e. 4 nearest neighbor edges for each
node on the diamond lattice).

    SGraph.jl
Defines a Simulation-Graph, which describes an extensive lattice used in the
simulation. The generation of this graph requires an RGraph.

    model.jl
This file contains most of the "low-level" Monte-Carlo code, i.e. the functions
required to perform a sweep.

    simulation.jl
Here, higher level functions around the Monte-Carlo code are defined. This
includes file writers, Binning Analysis, Jackknife and functions to run the
simulation.
Note: main.jl and mainf.jl contain more functionality to simplify running the
      simulation.
=#


end # module
