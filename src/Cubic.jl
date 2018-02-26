#=
This file is kinda unnecessary.

Here a lattice is defined through a cubic unit cell, given by type Cubic. Each
Cubic consists out of one or more types "Bravais", which represent a cubic
"Bravais" lattice. Each type "Bravais" is given by an anchor position (basis)
"pos" and three vectors x, y, z generating lattice positions.

A few basic cubic lattices are already defined, namely simple cubic "sc()", body
centered cubic "bcc()", face centered cubic "fcc()" and "diamond()".
=#

type Bravais
    pos::Point3f0
    x::Vec3f0
    y::Vec3f0
    z::Vec3f0
end

typealias Cubic Dict{_String, Bravais}


################################################################################
#### Constructors
################################################################################

Cubic(key::_String, B::Bravais) = Cubic(key => B)

# """
#     sc(atom::_String)
#     sc(atom::_String, size::Integer)
#     sc(atom::_String, w::Integer, h::Integer, l::Integer)
#
# Constructor for a simple cubic cell. A cubic cell will always be constructed with
# scaling factor w, h, l, which are 1 by default.
# """

# Bravais concstructors
function sc(pos::Point3f0)
    Bravais(
        pos,
        Vec3f0(1, 0, 0),
        Vec3f0(0, 1, 0),
        Vec3f0(0, 0, 1)
    )
end

function sc(pos::Point3f0, scaling::Float32)
    Bravais(
        pos,
        scaling * Vec3f0(1, 0, 0),
        scaling * Vec3f0(0, 1, 0),
        scaling * Vec3f0(0, 0, 1)
    )
end

sc(scaling::Float32) = sc(Point3f0(0, 0, 0), scaling)
sc() = sc(Point3f0(0, 0, 0))

# Cubic constructors
sc(atom::_String, pos::Point3f0, scaling::Float32) = Cubic(atom, sc(pos, scaling))
sc(atom::_String, pos::Point3f0) = Cubic(atom, sc(pos))
sc(atom::_String, scaling::Float32) = Cubic(atom, sc(scaling))
sc(atom::_String) = Cubic(atom, sc())


# """
#     bcc(atom::_String)
#     bcc(atom::_String, size::Integer)
#     bcc(atom::_String, w::Integer, h::Integer, l::Integer)
#
# Constructor for a body centric cubic cell. A cubic cell will always be
# constructed with scaling factor w, h, l, which are 1 by default.
# """
# Bravais constructor
function bcc(pos::Point3f0)
    Bravais(
        pos,
        Vec3f0(0.5, 0.5, -0.5),
        Vec3f0(0.5, -0.5, 0.5),
        Vec3f0(-0.5, 0.5, 0.5)
    )
end

function bcc(pos::Point3f0, scaling::Float32)
    Bravais(
        pos,
        scaling * Vec3f0(0.5, 0.5, -0.5),
        scaling * Vec3f0(0.5, -0.5, 0.5),
        scaling * Vec3f0(-0.5, 0.5, 0.5)
    )
end

bcc(scaling::Float32) = bcc(Point3f0(0, 0, 0), scaling)
bcc() = bcc(Point3f0(0, 0, 0))

# Cubic constructors
bcc(atom::_String, pos::Point3f0, scaling::Float32) = Cubic(atom, bcc(pos, scaling))
bcc(atom::_String, pos::Point3f0) = Cubic(atom, bcc(pos))
bcc(atom::_String, scaling::Float32) = Cubic(atom, bcc(scaling))
bcc(atom::_String) = Cubic(atom, bcc())


# """
#     fcc(atom::_String)
#     fcc(atom::_String, size::Integer)
#     fcc(atom::_String, w::Integer, h::Integer, l::Integer)
#
# Constructor for a face centric cubic cell. A cubic cell will always be
# constructed with scaling factor w, h, l, which are 1 by default.
# """
# Bravais constructor
function fcc(pos::Point3f0)
    Bravais(
        pos,
        Vec3f0(0.5, 0.5, 0.0),
        Vec3f0(0.5, 0.0, 0.5),
        Vec3f0(0.0, 0.5, 0.5)
    )
end

function fcc(pos::Point3f0, scaling::Float32)
    Bravais(
        pos,
        scaling * Vec3f0(0.5, 0.5, 0.0),
        scaling * Vec3f0(0.5, 0.0, 0.5),
        scaling * Vec3f0(0.0, 0.5, 0.5)
    )
end

fcc(scaling::Float32) = fcc(Point3f0(0, 0, 0), scaling)
fcc() = fcc(Point3f0(0, 0, 0))

# Cubic constructors
fcc(atom::_String, pos::Point3f0, scaling::Float32) = Cubic(atom, fcc(pos, scaling))
fcc(atom::_String, pos::Point3f0) = Cubic(atom, fcc(pos))
fcc(atom::_String, scaling::Float32) = Cubic(atom, fcc(scaling))
fcc(atom::_String) = Cubic(atom, fcc())


# """
#     diamond(atom::_String)
#     diamond(atom::_String, size::Integer)
#     diamond(atom::_String, w::Integer, h::Integer, l::Integer)
#
# Constructor for a diamond cubic cell. A cubic cell will always be constructed
# with scaling factor w, h, l, which are 1 by default.
# """
function diamond(atom::_String, pos::Point3f0)
    #println(typeof(pos))
    Cubic(
        atom*" 1" => fcc(pos),
        atom*" 2" => fcc(pos + Point3f0(0.25, 0.25, 0.25))
    )
end

function diamond(atom::_String, pos::Point3f0, scaling::Float32)
    #println(typeof(pos))
    Cubic(
        atom*" 1" => fcc(pos, scaling),
        atom*" 2" => fcc(pos + Point3f0(0.25, 0.25, 0.25), scaling)
    )
end

diamond(atom::_String, scaling::Float32) = diamond(atom, Point3f0(0, 0, 0), scaling)
diamond(atom::_String) = diamond(atom, Point3f0(0, 0, 0))

################################################################################
#### Utilities
################################################################################

typealias Vec3i Vec{3, Integer}


function *(uvw::Vec3i, B::Bravais)
    B.pos + uvw[1] * B.x + uvw[2] * B.y + uvw[3] * B.z
end


# Put the origin of Bravais B.pos as close to Point p as possible.
"""
    center_around!(B::Bravais, p::Point3f0=Point3f0(0, 0, 0))

Changes the basis of Bravais to be as close to p as possible.
"""
function center_around!(B::Bravais, p::Point3f0=Point3f0(0, 0, 0))
    return center_rec(B, p, Vec3i(0, 0, 0))
end

function center_rec(B::Bravais, p::Point3f0, uvw::Vec3i)
    dist2 = dot(B.pos - p, B.pos - p)

    for x in -1:1, y in -1:1, z in -1:1
        if ((x == 0) && (y == 0) && (z == 0)) continue end
        new_pos = B.pos + x * B.x + y * B.y + z * B.z
        if dot(new_pos - p, new_pos - p) < dist2
            #println("(", x, ", ", y, ", ", z, ") -> ", new_pos)
            uvw += Vec3i(x, y, z)
            B.pos = new_pos
            return center_rec(B, p, uvw)
        end
    end

    return uvw
end

################################################################################
#### Combining
################################################################################

# Combines multiple Cubics
function combine(cubics::Cubic...)
    out = Cubic()

    for c in cubics
        for key in keys(c)
            if key in keys(out)
                throw(ErrorException(string(
                    "Failed to combine Cubics for atom ", key, ". ",
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
