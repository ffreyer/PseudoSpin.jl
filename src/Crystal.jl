# As in R = R_0 + u e_x + v e_y + w e_z
mutable struct Bravais
    pos::SVector{3, Float32}
    x::SVector{3, Float32}
    y::SVector{3, Float32}
    z::SVector{3, Float32}
end

# typealias Crystal Dict{String, Bravais}

# Crystal are combination of Bravais lattices
const Crystal = Dict{String, Bravais}


################################################################################
#### Constructors
################################################################################


"""
    Crystal(atom, Bravais)

Returns a new crystal with atoms positioned according to the given Bravais
lattice.
"""
Crystal(atom::String, B::Bravais) = Crystal(atom => B)


# Bravais concstructors
"""
    sc(, [pos::SVector{3, Float32}, scaling::Float32])

Generates a simple cubic Bravais lattice. Without arguments, the Bravais lattice
will start at (0, 0, 0) and be scaled by 1. Arguments are evaluated based on
their type.
"""
function sc(pos::SVector{3, Float32})
    Bravais(
        pos,
        SVector{3, Float32}(1, 0, 0),
        SVector{3, Float32}(0, 1, 0),
        SVector{3, Float32}(0, 0, 1)
    )
end
function sc(pos::SVector{3, Float32}, scaling::Float32)
    Bravais(
        pos,
        scaling * SVector{3, Float32}(1, 0, 0),
        scaling * SVector{3, Float32}(0, 1, 0),
        scaling * SVector{3, Float32}(0, 0, 1)
    )
end
sc(scaling::Float32) = sc(SVector{3, Float32}(0, 0, 0), scaling)
sc() = sc(SVector{3, Float32}(0, 0, 0))

# Crystal constructors
"""
    sc(atom::String, [pos::SVector{3, Float32}, scaling::Float32])

Generates a Crystal with atoms positioned according to a simple cubic lattice.
Without arguments, the Bravais lattice will start at (0, 0, 0) and be scaled by
1. Arguments are evaluated based on their type.
"""
sc(atom::String, pos::SVector{3, Float32}, scaling::Float32) = Crystal(atom, sc(pos, scaling))
sc(atom::String, pos::SVector{3, Float32}) = Crystal(atom, sc(pos))
sc(atom::String, scaling::Float32) = Crystal(atom, sc(scaling))
sc(atom::String) = Crystal(atom, sc())


# Bravais constructor
"""
    bcc(, [pos::SVector{3, Float32}, scaling::Float32])

Generates a body centered cubic Bravais lattice. Without arguments, the Bravais
lattice will start at (0, 0, 0) and be scaled by 1. Arguments are evaluated
based on their type.
"""
function bcc(pos::SVector{3, Float32})
    Bravais(
        pos,
        SVector{3, Float32}(0.5, 0.5, -0.5),
        SVector{3, Float32}(0.5, -0.5, 0.5),
        SVector{3, Float32}(-0.5, 0.5, 0.5)
    )
end
function bcc(pos::SVector{3, Float32}, scaling::Float32)
    Bravais(
        pos,
        scaling * SVector{3, Float32}(0.5, 0.5, -0.5),
        scaling * SVector{3, Float32}(0.5, -0.5, 0.5),
        scaling * SVector{3, Float32}(-0.5, 0.5, 0.5)
    )
end
bcc(scaling::Float32) = bcc(SVector{3, Float32}(0, 0, 0), scaling)
bcc() = bcc(SVector{3, Float32}(0, 0, 0))

# Crystal constructors
"""
    bcc(atom::String, [pos::SVector{3, Float32}, scaling::Float32])

Generates a Crystal with atoms positioned according to a body centered cubic
lattice. Without arguments, the Bravais lattice will start at (0, 0, 0) and be
scaled by 1. Arguments are evaluated based on their type.
"""
bcc(atom::String, pos::SVector{3, Float32}, scaling::Float32) = Crystal(atom, bcc(pos, scaling))
bcc(atom::String, pos::SVector{3, Float32}) = Crystal(atom, bcc(pos))
bcc(atom::String, scaling::Float32) = Crystal(atom, bcc(scaling))
bcc(atom::String) = Crystal(atom, bcc())


# Bravais constructor
"""
    fcc(, [pos::SVector{3, Float32}, scaling::Float32])

Generates a face centered cubic Bravais lattice. Without arguments, the Bravais
lattice will start at (0, 0, 0) and be scaled by 1. Arguments are evaluated
based on their type.
"""
function fcc(pos::SVector{3, Float32})
    Bravais(
        pos,
        SVector{3, Float32}(0.5, 0.5, 0.0),
        SVector{3, Float32}(0.5, 0.0, 0.5),
        SVector{3, Float32}(0.0, 0.5, 0.5)
    )
end
function fcc(pos::SVector{3, Float32}, scaling::Float32)
    Bravais(
        pos,
        scaling * SVector{3, Float32}(0.5, 0.5, 0.0),
        scaling * SVector{3, Float32}(0.5, 0.0, 0.5),
        scaling * SVector{3, Float32}(0.0, 0.5, 0.5)
    )
end
fcc(scaling::Float32) = fcc(SVector{3, Float32}(0, 0, 0), scaling)
fcc() = fcc(SVector{3, Float32}(0, 0, 0))

# Crystal constructors
"""
    fcc(atom::String, [pos::SVector{3, Float32}, scaling::Float32])

Generates a Crystal with atoms positioned according to a face centered cubic
lattice. Without arguments, the Bravais lattice will start at (0, 0, 0) and be
scaled by 1. Arguments are evaluated based on their type.
"""
fcc(atom::String, pos::SVector{3, Float32}, scaling::Float32) = Crystal(atom, fcc(pos, scaling))
fcc(atom::String, pos::SVector{3, Float32}) = Crystal(atom, fcc(pos))
fcc(atom::String, scaling::Float32) = Crystal(atom, fcc(scaling))
fcc(atom::String) = Crystal(atom, fcc())


"""
    diamond(atom::String, [pos::SVector{3, Float32}, scaling::Float32])

Generates a Crystal with atoms positioned according to a diamond lattice.
Without arguments, the Bravais lattice will start at (0, 0, 0) and be scaled by
1. Arguments are evaluated based on their type.
"""
function diamond(atom::String, pos::SVector{3, Float32})
    #println(typeof(pos))
    Crystal(
        atom*" 1" => fcc(pos),
        atom*" 2" => fcc(pos + SVector{3, Float32}(0.25, 0.25, 0.25))
    )
end
function diamond(atom::String, pos::SVector{3, Float32}, scaling::Float32)
    #println(typeof(pos))
    Crystal(
        atom*" 1" => fcc(pos, scaling),
        atom*" 2" => fcc(pos + SVector{3, Float32}(0.25, 0.25, 0.25), scaling)
    )
end
diamond(atom::String, scaling::Float32) = diamond(atom, SVector{3, Float32}(0, 0, 0), scaling)
diamond(atom::String) = diamond(atom, SVector{3, Float32}(0, 0, 0))


################################################################################
#### Utilities
################################################################################


# computes R = R_0 + u e_x + v e_y + w e_z for a given Bravais lattice
function *(uvw::SVector{3, Int64}, B::Bravais)
    B.pos + uvw[1] * B.x + uvw[2] * B.y + uvw[3] * B.z
end


# Put the origin of Bravais B.pos as close to Point p as possible.
"""
    center_around!(B::Bravais, p::SVector{3, Float32}=SVector{3, Float32}(0, 0, 0))

Changes the basis of Bravais to be as close to p as possible.
"""
function center_around!(B::Bravais, p::SVector{3, Float32}=SVector{3, Float32}(0, 0, 0))
    return center_rec(B, p, SVector{3, Int64}(0, 0, 0))
end

function center_rec(B::Bravais, p::SVector{3, Float32}, uvw::SVector{3, Int64})
    dist2 = dot(B.pos - p, B.pos - p)

    for x in -1:1, y in -1:1, z in -1:1
        if ((x == 0) && (y == 0) && (z == 0)) continue end
        new_pos = B.pos + x * B.x + y * B.y + z * B.z
        if dot(new_pos - p, new_pos - p) < dist2
            #println("(", x, ", ", y, ", ", z, ") -> ", new_pos)
            uvw += SVector{3, Int64}(x, y, z)
            B.pos = new_pos
            return center_rec(B, p, uvw)
        end
    end

    return uvw
end


################################################################################
#### Combining
################################################################################


# Combines multiple Crystals
function combine(cubics::Crystal...)
    out = Crystal()

    for c in cubics
        for key in keys(c)
            if key in keys(out)
                throw(ErrorException(string(
                    "Failed to combine Crystals for atom ", key, ". ",
                    "Use different names (X1, X2, ...) if the atoms do not form",
                    " a Bravais lattice or initialise one instead."
                )))
            else
                push!(out, Pair(key, c[key]))
            end
        end
    end

    out
end
