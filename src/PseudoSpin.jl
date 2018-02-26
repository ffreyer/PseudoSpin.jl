module PseudoSpin

__precompile__(true)


# Bleh
if VERSION == v"0.4.5"
    typealias _String ASCIIString
else
    typealias _String String
end


# Small FixedSizeArrays, e.g. 3-component vectors, outperform Julia Arrays.
if Pkg.status("FixedSizeArrays") != nothing
    using FixedSizeArrays
elseif Pkg.status("StaticArrays") != nothing
    using StaticArrays.FixedSizeArrays
end
typealias Vec3 Vec{3}
typealias Point3 Point{3}
typealias Point3f0 Point{3, Float32}
typealias Vec3f0 Vec{3, Float32}


# I'm implementing new methods for these functions
import Base.*, Base.==, Base.in, Base.findfirst#, Base.show
import Base.start, Base.next, Base.length, Base.done, Base.eltype
import Base.push!, Base.var


# Files used in simulation
include("Cubic.jl")
include("RGraph.jl")
include("SGraph.jl")
include("model.jl")
include("simulation.jl")

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
