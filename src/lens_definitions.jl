

export Ray, SurfBase, Trace, OptSurface, ModelSurface, ExtendedGeometry
export OpticalSystem, SurfProfileOAConic, SizeLens, RoundAperture, RectAperture
export NoProfile, NoBendIndex, NoAmpParam, updateCoordChange, AbstractSurface 
export AbstractAmplitudeParam, DielectricT,MirrorR,CDiffuser, NoBendIndex
export AmpData


"""
    AbstractRay{N, T<:Real}

Abstract supertype for ray types. `Ray{N,T}` is currently its only
subtype; nothing else in the codebase dispatches on `AbstractRay`
directly -- everything that takes a ray uses the concrete `Ray` type.
"""
abstract type AbstractRay{N, T<:Real} end

"""
    AbstractSurfBase{N, T<:Real}

Abstract supertype for surface-placement/orientation types.
`SurfBase{N,T}` is currently its only subtype; nothing else in the
codebase dispatches on `AbstractSurfBase` directly.
"""
abstract type AbstractSurfBase{N, T<:Real} end

"""
    AbstractSize{T<:Real}

Abstract supertype for surface "size" types -- how large a surface's
aperture/extent is. Subtypes: `SizeLens` (a lens's semi-diameter, not a
clipping aperture) and `TrueAperture` (`RoundAperture`, `RectAperture`,
which do clip). See `isAperture` (`src/aperture.jl`), which dispatches
on this distinction.
"""
abstract type AbstractSize{T<:Real} end

"""
    AbstractSurfProfile{T<:Real}

Abstract supertype for surface-sag profile types (`SurfProfileConic`,
`SurfProfileSphere`, `SurfProfileAsphere`, ...). `sag`, `deltaToSurf`,
and `surfNormal` (`src/tracing.jl`) are implemented per concrete
subtype.
"""
abstract type AbstractSurfProfile{T<:Real} end

"""
    AbstractAmplitudeParam

Abstract supertype for amplitude/coating parameter types (`AmpParam`,
`NoAmpParam`). See also [`APorString`](@ref), the
`Union{AbstractAmplitudeParam, String}` type used to let a coating be
passed either as one of these objects or as a raw string key.
"""
abstract type AbstractAmplitudeParam end

"""
    AbstractBendType{T<:Real}

Abstract supertype for "mod" types describing how a ray's direction
changes at a surface (refract, reflect, diffuse, or pass through
unchanged). Concrete subtypes: `DielectricT`/`MirrorR` (via
`AbstractBendDielectric`/`AbstractBendMirror`), `CDiffuser`,
`NoBendIndex`. `modFunc` (`src/tracing.jl`, `src/surfaces.jl`) is
implemented per subtype.
"""
abstract type AbstractBendType{T<:Real} end #some info on how to do the bending

"""
    AbstractBendDielectric{T}

Abstract supertype for refracting `AbstractBendType`s. `DielectricT` is
its only subtype.
"""
abstract type AbstractBendDielectric{T} <: AbstractBendType{T} end

"""
    AbstractBendMirror{T}

Abstract supertype for reflecting `AbstractBendType`s. `MirrorR` is its
only subtype.
"""
abstract type AbstractBendMirror{T} <: AbstractBendType{T} end

"""
    AbstractOpticalObject{T<:Real}

Abstract type with no subtypes yet. `ExtendedGeometry.geo` is typed
`Array{AbstractSurface}` with a comment noting it "needs to be changed
to `AbstractOpticalObject`" (see `TODO.md`) -- this type appears
reserved for that future use but isn't used anywhere currently.
"""
abstract type AbstractOpticalObject{T<:Real} end

"""
    AbstractTrace{T<:Real}

Abstract type with no subtypes; `Trace` does not subtype it despite the
similar name. Currently unused.
"""
abstract type AbstractTrace{T<:Real} end

"""
    APorString

`Union{AbstractAmplitudeParam, String}`. The type used for `coating`
parameters throughout `src/surfaces.jl`: a coating can be passed either
as an already-built `AbstractAmplitudeParam` object, or as a raw string
name/key to look up (or create) in the `attributesSurfaces` dictionary
-- see `getAmpParams` in `src/tracing.jl`.

Not `const`-declared (unlike the file-level constants in
`src/constants.jl`) -- it depends on `AbstractAmplitudeParam`, defined
above in this same file, so it can't be relocated to `constants.jl`
(which loads first) without also moving the abstract type.
"""
APorString = Union{AbstractAmplitudeParam, String}

"""
    Ray{N, T}

A ray: a `base` point and a `dir` direction vector, both `N`-dimensional.
The fundamental unit propagated through a geometry by
`traceGeometry`/`traceSurf` (`src/tracing.jl`).

Fields:
- `base::Point{N, T}` -- the ray's starting point
- `dir::Vec{N, T}` -- the ray's direction vector
"""
struct Ray{N, T} <: AbstractRay{N, T}
    base::Point{N, T}
    dir::Vec{N, T}

end

"""
    SurfBase{N, T}

The placement and local orientation of a surface: a base point plus the
local z (`dir`) and y (`ydir`) axes of the surface's coordinate frame
(the local x axis is derived from these -- see `findPerpenMap` in
`src/extended_geo.jl`). Used as the `base` field of `OptSurface`/
`ModelSurface`, and updated by `updateCoordChange!`.

Fields:
- `base::Point{N, T}` -- the surface's origin, in global coordinates
- `dir::Vec{N, T}` -- the surface's local z direction (its normal, for
  a plane), in global coordinates
- `ydir::Vec{N, T}` -- the surface's local y direction, in global
  coordinates
"""
mutable struct SurfBase{N, T} <: AbstractSurfBase{N, T}
    base::Point{N, T}
    dir::Vec{N, T}
    ydir::Vec{N, T}
end


"""
    SizeLens{T}

The size of a round optical surface, given only as a semi-diameter --
not a clipping aperture (`isAperture(::SizeLens) == false`, see
`src/aperture.jl`; compare `RoundAperture`/`RectAperture`, which do
clip).

Fields:
- `semiDiameter::T`
"""
mutable struct SizeLens{T} <: AbstractSize{T}  #size for round optics
    semiDiameter::T
end

"""
    TrueAperture{T}

Abstract supertype for size types that represent an actual clipping
aperture/obscuration (as opposed to `SizeLens`, which is just a size).
`isAperture(a::TrueAperture) = true` (`src/aperture.jl`). Subtypes:
`RoundAperture`, `RectAperture`.
"""
abstract type TrueAperture{T} <: AbstractSize{T} end

"""
    RoundAperture{T}

A round aperture with an optional central circular obscuration. Used
by `clipAperture` (`src/aperture.jl`): a point is clipped if its radius
is inside `obscure` or outside `semiDiameter`.

Fields:
- `obscure::T` -- radius of the central obscuration (0 for none)
- `semiDiameter::T` -- outer clear-aperture radius
"""
mutable struct RoundAperture{T} <: TrueAperture{T}
    obscure::T
    semiDiameter::T
end

"""
    RectAperture{T}

A rectangular aperture with an optional central rectangular
obscuration. Used by `clipAperture` (`src/aperture.jl`): a point is
clipped if it falls within the inner `wo`x`lo` obscuration rectangle,
or outside the outer `wclear`x`lclear` clear-aperture rectangle.

Fields:
- `wo::T`, `lo::T` -- half-width/half-length of the central obscuration
  (0 for none)
- `wclear::T`, `lclear::T` -- half-width/half-length of the outer clear
  aperture
"""
mutable struct RectAperture{T} <: TrueAperture{T}
    wo::T
    lo::T
    wclear::T
    lclear::T
end


"""
    SurfProfileSphere{T}

A spherical surface profile.

Fields:
- `curv::T` -- curvature, `1/radius`
"""
mutable struct SurfProfileSphere{T} <: AbstractSurfProfile{T}
    curv :: T
end

"""
    conicToϵ(c)

Convert a Zemax-style conic constant `c` to the `ϵ` conic parameter
used throughout this package's `sag`/`deltaToSurf`/`surfNormal`
formulas (`ϵ = 1 + c`, following Welford's convention). Used e.g. when
importing Zemax surfaces (`src/zemax.jl`).
"""
conicToϵ(c) = 1 + c #convert the conic constant to ϵ used in Welford

"""
    SurfProfileConic{T}

A conic surface profile (sphere, paraboloid, ellipsoid, or hyperboloid,
depending on `ϵ`).

Fields:
- `curv::T` -- curvature, `1/radius`
- `ϵ::T` -- conic parameter (see [`conicToϵ`](@ref); `ϵ=1` is
  spherical)
"""
mutable struct SurfProfileConic{T} <: AbstractSurfProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
end

"""
    SurfProfileOAConic{T}

An off-axis conic surface profile: a `SurfProfileConic` evaluated about
a point offset from the local origin (see the `sag` method for this
type in `src/tracing.jl`, which shifts `x,y` by `offset[1:2]` and adds
`offset[3]` to the result).

Fields:
- `curv::T` -- curvature, `1/radius`
- `ϵ::T` -- conic parameter (see [`conicToϵ`](@ref))
- `offset::Vec{3, T}` -- offset of the conic's own vertex/axis from the
  local origin
"""
mutable struct SurfProfileOAConic{T} <: AbstractSurfProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
    offset::Vec{3, T}
end

"""
    AbstractAsphericProfile{T}

Abstract supertype for aspheric surface-profile types that add
polynomial correction terms to a base conic. Subtypes:
`SurfProfileAsphere`, `SurfProfileEvenAsphere`, `SurfProfileOddAsphere`,
`SurfProfileXYPoly`. Every subtype exposing `curv`/`ϵ` fields (all four
above do) gets `deltaToSurf` (`src/tracing.jl`) and `reverseProfile!`
(`src/surface_manipulation.jl`) for free from generic methods dispatched
on this abstract type; subtypes additionally exposing an `a::Vector{T}`
polynomial-coefficient field also get `surfNormal` for free (via a
generic `ForwardDiff`-based fallback, `src/tracing.jl`) unless they
define a more specific method of their own (`SurfProfileAsphere`/
`SurfProfileEvenAsphere` do; `SurfProfileOddAsphere`/`SurfProfileXYPoly`
don't, relying on the fallback).
"""
abstract type AbstractAsphericProfile{T} <: AbstractSurfProfile{T} end


"""
    SurfProfileAsphere{T}

A general aspheric surface profile: a conic base plus polynomial
correction terms starting from 3rd order (see the `sag` method for this
type in `src/tracing.jl`).

Fields:
- `curv::T` -- curvature, `1/radius`
- `ϵ::T` -- conic parameter (see [`conicToϵ`](@ref))
- `a::Vector{T}` -- polynomial coefficients, starting at 3rd order
"""
mutable struct SurfProfileAsphere{T} <: AbstractAsphericProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
    a::Vector{T} #hope this works...
end


"""
    SurfProfileEvenAsphere{T}

An even-aspheric surface profile: a conic base plus polynomial
correction terms restricted to even orders, starting from 4th order
(the common lens-manufacturing convention; see the `sag` method for
this type in `src/tracing.jl`).

Fields:
- `curv::T` -- curvature, `1/radius`
- `ϵ::T` -- conic parameter (see [`conicToϵ`](@ref))
- `a::Vector{T}` -- even-order polynomial coefficients, starting at 4th
  order
"""
mutable struct SurfProfileEvenAsphere{T} <: AbstractAsphericProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
    a::Vector{T} #even terms starting with 4th order
end

"""
    SurfProfileOddAsphere{T}

An "odd asphere" surface profile (Zemax `TYPE ODDASPHE`): a conic base
plus polynomial correction terms starting from **1st** order (see the
`sag` method for this type in `src/tracing.jl`) -- unlike
`SurfProfileAsphere` (starts at 3rd order) or `SurfProfileEvenAsphere`
(starts at 4th order, even only), this one can represent a pure linear
(`r¹`) term, needed e.g. for an axicon's conical profile. Confirmed
against Zemax's own `ODDASPHE` samples: `Parameter i` multiplies `r^i`
directly for `i = 1..8`, so no reindexing is needed when importing
(compare `SurfProfileEvenAsphere`'s Zemax import, which drops
`Parameter 1`).

Fields:
- `curv::T` -- curvature, `1/radius`
- `ϵ::T` -- conic parameter (see [`conicToϵ`](@ref))
- `a::Vector{T}` -- polynomial coefficients: `a[i]` is the coefficient
  of `r^i`, starting at `i=1`
"""
mutable struct SurfProfileOddAsphere{T} <: AbstractAsphericProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
    a::Vector{T} #a[i] is the coefficient of r^i, starting at i=1
end

"""
    SurfProfileXYPoly{T}

A general 2D (non-rotationally-symmetric) polynomial surface profile
(Zemax `TYPE XPOLYNOM`): a conic base plus a sum of `x^m y^n` terms in
*normalized* coordinates `(x/normRadius, y/normRadius)` (see the `sag`
method for this type in `src/tracing.jl`, and [`xyPolyTermPowers`](@ref)
for how a 1-based term index maps to its `(m,n)` exponent pair). Unlike
every other `AbstractAsphericProfile` subtype, this one is not
rotationally symmetric about the local z axis.

Fields:
- `curv::T` -- curvature of the base conic, `1/radius`
- `ϵ::T` -- conic parameter of the base conic (see [`conicToϵ`](@ref))
- `normRadius::T` -- normalization radius the polynomial's `x`/`y`
  arguments are divided by before being raised to a power
- `a::Vector{T}` -- polynomial term coefficients, in Zemax's own
  bivariate term ordering (`a[1]` -> `x`, `a[2]` -> `y`, `a[3]` -> `x²`,
  `a[4]` -> `xy`, `a[5]` -> `y²`, `a[6]` -> `x³`, ... -- see
  [`xyPolyTermPowers`](@ref))
"""
mutable struct SurfProfileXYPoly{T} <: AbstractAsphericProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
    normRadius::T
    a::Vector{T} #Zemax XY-polynomial term coefficients, see xyPolyTermPowers
end

"""
    SurfProfileCyl{T}

A cylindrical surface profile: a conic cross-section along one axis,
extruded along the other, plus polynomial correction terms (see the
`sag` method for this type in `src/tracing.jl`).

Fields:
- `curv::T` -- curvature, `1/radius`
- `ϵ::T` -- conic parameter (see [`conicToϵ`](@ref))
- `a::Vector{T}` -- polynomial correction coefficients
"""
mutable struct SurfProfileCyl{T} <: AbstractSurfProfile{T}
    curv :: T
    ϵ::T #see Welford for definition of ϵ
#   ϵy::T #see Welford for definition of ϵ
    a::Vector{T} #hope this works...
end

"""
    SurfProfileToroid{T}

A toroidal surface profile (Zemax `TYPE TOROIDAL`): a conic curve in
the local y-z plane (`curvY`/`ϵY`), independently swept by a circular
arc along local x (`curvX`) -- see the `sag` method for this type in
`src/tracing.jl` for the exact formula. `curvX == 0` means no x sweep
at all (a pure extrusion of the y-z curve along x, i.e. a cylinder) --
matches Zemax's own convention where a `TOROIDAL` surface's "Radius of
Rotation" parameter of `0` means the same thing.

Fields:
- `curvY::T` -- curvature of the base curve along the local y
  direction, `1/Ry`
- `ϵY::T` -- conic parameter of the base y-z curve (see
  [`conicToϵ`](@ref))
- `curvX::T` -- curvature of the circular sweep along the local x
  direction, `1/Rx` (`0` disables the x sweep entirely, rather than
  meaning a literal zero-radius sweep)
"""
mutable struct SurfProfileToroid{T} <: AbstractSurfProfile{T}
    curvY :: T
    ϵY :: T
    curvX :: T
end

"""
    NoProfile{T}

Placeholder flat-plane profile (normal `(0,0,1)`, trivial
`Δ = -z0/dir_z` intersection) used for non-refracting reference and
aperture surfaces -- e.g. `referencePlane`, `cDiffuser`
(`src/surfaces.jl`), and the `Washer`/`Disk` aperture shapes
(`src/mesh_primitives.jl`). See the dedicated `sag`/`surfNormal`/
`deltaToSurf` methods for this type in `src/surfaces.jl`/`src/tracing.jl`.

Fields:
- `curv::T` -- present for type-parameter uniformity with other
  profiles, but never actually used (always constructed as `0.`)
"""
mutable struct NoProfile{T} <: AbstractSurfProfile{T}
    curv :: T #but never used
end

"""
    ParaxialProfile{T}

Placeholder flat-plane profile for an ideal thin lens (Zemax `TYPE
PARAXIAL`) -- identical `sag`/`deltaToSurf` to [`NoProfile`](@ref) (a
paraxial lens has no real sag), but a **deliberately different**
`surfNormal` (`src/tracing.jl`): it returns the local intersection
*coordinates* `(x, y, 0)`, not a true unit normal. An ideal thin lens's
ray-bending depends on where a ray hits it, not on any actual surface
normal, so `surfNormal` here repurposes the one channel already routed
from `traceSurf` through to `modFunc` (`s.toGlobalDir` converts
whatever `surfNormal` returns into global coordinates either way) to
carry that position instead -- see `modFunc(ray, offset,
::ParaxialLensT)`, which is the only bend type meant to ever pair with
this profile. This is the one profile type in the package whose
`surfNormal` is *not* a true normal; because `GeometryBasics.normals`
(`src/mesh_primitives.jl`) also calls the generic `surfNormal`, that
file has a more specific override for this profile so mesh shading
still gets the real `(0,0,1)` normal instead of the repurposed value.

Fields:
- `curv::T` -- present for type-parameter uniformity with other
  profiles, but never actually used (always constructed as `0.`)
"""
mutable struct ParaxialProfile{T} <: AbstractSurfProfile{T}
    curv :: T #but never used
end


"""
    AmpParam

An amplitude/coating parameter identified by name. Constructed directly
for a handful of surfaces (e.g. `cDiffuser`, `src/surfaces.jl`); most
surface constructors instead go through `getAmpParams`
(`src/tracing.jl`), which looks a coating name up in (or adds it to)
the global `attributesSurfaces` dictionary.

Fields:
- `type::String` -- the coating name/key. Stored but not otherwise
  read or dispatched on anywhere in the codebase.
"""
mutable struct AmpParam <: AbstractAmplitudeParam
    type::String
end

"""
    AbstractAmpData{T<:Real}

Abstract supertype for amplitude/polarization transfer data attached to
a `Trace`. `AmpData` is its only subtype; used as the
`S<:AbstractAmpData{T}` type parameter of [`Trace`](@ref).
"""
abstract type AbstractAmpData{T <:Real} end

"""
    AmpData{T}

Polarization/amplitude transfer data for one trace step: the `p`
(polarization) and `o` (orthogonal/"other") 3x3 transfer matrices, plus
a scalar transmission coefficient per wavelength. The identity case (a
surface that doesn't alter polarization) is `AMPPERFECTSURFACE`
(`src/tracing.jl`), built from [`identityPol`](@ref).

Fields:
- `p::SMatrix{3,3,T}` -- polarization transfer matrix
- `o::SMatrix{3,3,T}` -- orthogonal/"other" transfer matrix
- `trans::Array{T,1}` -- transmission coefficient(s), for intensity
  systems
"""
mutable struct AmpData{T} <: AbstractAmpData{T}
    p::SMatrix{3,3,T}
    o::SMatrix{3,3,T}
    trans::Array{T,1} #for intensity systems, the "transmission" coefficient of the surface
end


"""
    Trace(ray, nIn, delta, pmatrix)
    Element of a raytrace.
        ray     location of interesection with surface and exiting direction
        nIn     refractive index leaving surface
        delta   length of ray leaving surface
        pmatrix polarization p & o matrices
"""
mutable struct Trace{T <:Real, S <:AbstractAmpData{T}}
    ray::Ray{3, T}
    nIn::T
    delta::T
    pmatrix::S
end

#=

AbstractBendTypes

=#


"""
    DielectricT{T}

The `mod` ("bend") type for a refracting surface: a boundary between
two media with the given refractive indices. `modFunc`
(`src/tracing.jl`) computes the refracted ray direction (and flags
total internal reflection) for this type.

Fields:
- `refIndexIn::T` -- refractive index on the entering side
- `refIndexOut::T` -- refractive index on the exiting side
"""
mutable struct DielectricT{T} <: AbstractBendDielectric{T}
    refIndexIn :: T
    refIndexOut :: T
end

"""
    MirrorR{T}

The `mod` ("bend") type for a reflecting (mirror) surface. `modFunc`
(`src/tracing.jl`) computes the reflected ray direction for this type.
The two refractive-index fields mirror `DielectricT`'s shape so mirror
and dielectric surfaces share the same `OptSurface` field layout, even
though a mirror doesn't change refractive index.

Fields:
- `refIndexIn::T`
- `refIndexOut::T`
"""
mutable struct MirrorR{T} <: AbstractBendMirror{T}
    refIndexIn :: T
    refIndexOut :: T
end

"""
    ParaxialLensT{T}

The `mod` ("bend") type for an ideal thin lens (Zemax `TYPE PARAXIAL`):
bends a ray according to the paraxial thin-lens transfer law, using
`focalLength` alone, rather than Snell's law at an index boundary.
Meant to pair exclusively with [`ParaxialProfile`](@ref) -- `modFunc`
(`src/tracing.jl`) treats its `normal`-slot argument as the ray's
transverse offset from the optical axis (which is what
`ParaxialProfile`'s `surfNormal` actually returns), not a true normal,
so pairing this bend type with any other profile would be a bug.

Fields:
- `focalLength::T` -- the lens's focal length, from Zemax `PARM 1`
- `refIndexIn::T`, `refIndexOut::T` -- refractive indices either side
  (kept for interface consistency with every other `AbstractBendType`;
  the bending itself doesn't depend on them -- an ideal lens's power is
  given directly as a focal length, not derived from curvature+index)
"""
mutable struct ParaxialLensT{T} <: AbstractBendType{T}
    focalLength::T
    refIndexIn::T
    refIndexOut::T
end

"""
    CDiffuser{T}

The `mod` ("bend") type for a conical (random-scatter) diffuser
surface, as built by `cDiffuser` (`src/surfaces.jl`). `modFunc`
(`src/tracing.jl`) samples a random point in the unit disk, scales it
by `tanθ`, and uses it to perturb the transmitted ray direction within
a cone.

Fields:
- `tanθ::T` -- tangent of the diffuser's half-deflection-angle cone
- `refIndexIn::T`, `refIndexOut::T` -- refractive indices either side
"""
mutable struct CDiffuser{T} <: AbstractBendType{T}
    tanθ::T
    refIndexIn::T
    refIndexOut::T
end

"""
    NoBendIndex{T}

The `mod` ("bend") type for a pass-through surface: the ray direction
is unchanged (`modFunc`/`surfAmpFunc` in `src/surfaces.jl` are no-ops
for this type). Used by `referencePlane` (`src/surfaces.jl`) for
reference/model surfaces used in characterization rather than real
optics.

Fields:
- `refIndexIn::T`, `refIndexOut::T` -- typically equal, since no
  refraction occurs; see the `NoBendIndex(n)` convenience constructor
"""
mutable struct NoBendIndex{T} <: AbstractBendType{T}
    refIndexIn :: T
    refIndexOut::T
end

"""
    NoBendIndex(n)

Convenience constructor: `NoBendIndex(n, n)`, i.e. the same refractive
index on both sides.
"""
NoBendIndex(n) = NoBendIndex(n,n)

#=
function NoBendIndex(n::Float64)
    #NoBendIndex(Float64(n), Float64(n))
    NoBendIndex(n, n)
end
=#

"""
    NoAmpParam

The amplitude/coating parameter type used alongside [`NoBendIndex`](@ref)
for pass-through reference surfaces (e.g. `referencePlane`,
`src/surfaces.jl`); paired with a no-op `surfAmpFunc` method
(`src/surfaces.jl`).

Fields:
- `type::String` -- coating name/key, as in [`AmpParam`](@ref); not
  read or dispatched on anywhere in the codebase.
"""
mutable struct NoAmpParam <: AbstractAmplitudeParam
    type::String
end

"""
    AbstractSurface{N,T}

Abstract supertype for all surface types (`OptSurface`, `ModelSurface`),
subtyping `GeometryBasics.GeometryPrimitive{N,T}` so surfaces can be
rendered by Makie via the `GeometryBasics` overloads in
`src/mesh_primitives.jl`.
"""
abstract type AbstractSurface{N,T} <: GeometryBasics.GeometryPrimitive{N, T} end
#abstract type AbstractSurface <: GeometryBasics.AbstractGeometry{3, Float64} end

"""
    OptSurface{N,T,S,U,V,W}

A full optical surface: placement, aperture, profile, refractive
behavior ("mod"), coating, precomputed local<->global coordinate
transforms, and display color. The primary surface type traced by
`traceSurf`/`traceSurf!` (`src/tracing.jl`), which can refract, reflect,
or diffuse (status 0/1/2 -- normal/missed/TIR; this method does not
check `aperture`, see `TODO.md`).

Fields:
- `surfname::String` -- surface name, used for lookups (e.g.
  `tracenumFromName`) and printing
- `base::SurfBase{N,T}` -- placement/orientation
- `aperture::U` -- size/clipping aperture (`U<:AbstractSize{T}`)
- `profile::V` -- sag profile (`V<:AbstractSurfProfile{T}`)
- `mod::W` -- refract/reflect/diffuse behavior (`W<:AbstractBendType{T}`)
- `coating::S` -- amplitude/coating parameter (`S<:AbstractAmplitudeParam`)
- `toGlobalCoord::AffineMap`, `toLocalCoord::AffineMap` -- point
  transforms between local and global coordinates
- `toGlobalDir::LinearMap`, `toLocalDir::LinearMap` -- direction
  transforms between local and global coordinates
- `color` -- display color, used by `plotSurface3D!` (`src/plotting.jl`)
"""
mutable struct OptSurface{N,T, S <:AbstractAmplitudeParam, U<:AbstractSize{T}, V<:AbstractSurfProfile{T}, W<:AbstractBendType{T}} <: AbstractSurface{N,T} #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase{N, T}
    aperture::U
    profile::V
    mod::W
    coating::S
    toGlobalCoord::AffineMap
    toLocalCoord::AffineMap
    toGlobalDir::LinearMap
    toLocalDir::LinearMap
    color
end

"""
    ModelSurface{N,T,U,V}

A non-refracting reference/model surface: placement, aperture, profile,
a fixed refractive index (for OPD bookkeeping), and coordinate
transforms, but no `mod`/coating -- the ray passes straight through
unchanged. Traced by `traceSurf`/`traceSurf!` (`src/tracing.jl`), which
returns status 0 (normal), 1 (missed), or 3 (aperture-clipped); used
for e.g. reference planes and characterization surfaces.

Fields:
- `surfname::String` -- surface name
- `base::SurfBase{N,T}` -- placement/orientation
- `aperture::U` -- size/clipping aperture (`U<:AbstractSize{T}`)
- `profile::V` -- sag profile (`V<:AbstractSurfProfile{T}`)
- `refIndex::T` -- refractive index at this surface, used for OPD
  calculations
- `toGlobalCoord::AffineMap`, `toLocalCoord::AffineMap` -- point
  transforms between local and global coordinates
- `toGlobalDir::LinearMap`, `toLocalDir::LinearMap` -- direction
  transforms between local and global coordinates
- `color` -- display color
"""
mutable struct ModelSurface{N,T, U<:AbstractSize{T}, V<:AbstractSurfProfile{T}, } <: AbstractSurface{N,T}  #use the data to overload GemoetryBasics
    surfname::String
    base::SurfBase{N, T}
    aperture::U
    profile::V
    refIndex::T #needed for OPD calculations
    toGlobalCoord::AffineMap
    toLocalCoord::AffineMap
    toGlobalDir::LinearMap
    toLocalDir::LinearMap
    color
end



"""
    ExtendedGeometry(geo, funcGeo, funcSetup, surfaceObject, wavelenth, parameters)

    Struct of the information to full characterize an optical system
    including pointers to functions that create the geometry and setup the values
    It includes a dictionary for definitions that can be used by the functions
        geo         array of AbstractSurfaces, the geometry. typically used for analysis
        funcGeo     the function that generates geo. enables adjustable parameters
        funcSetup   function called to set up conditions for funcGeo
        surfaceObject   define an Object to analyze
        wavelength  array of wavelengths in microns
        parameters  dictionary for parameters used by funcSetup and funcGeo

"""
mutable struct ExtendedGeometry
    geo::Array{AbstractSurface}  #needs to be changed to AbstractOpticalObject
    funcGeo
    funcSetup
    surfaceObject::AbstractSurface
    wavelength::Array{Float64}
    parameters::Dict
end

"""
    OpticalSystem{T<:Real}

A materialized optical system: the traceable geometry plus the
system-level metadata that isn't itself a surface. Unlike
[`ExtendedGeometry`](@ref) (a *parametric*, regenerable system --
`funcGeo`/`funcSetup` + a parameter dict, used by `plotOPD!`/
`plotOPD3D!`/`characterization.jl` to rebuild geometry on demand),
`OpticalSystem` represents an already-materialized result with no
function pointers or regeneration story -- currently produced by
[`readZemaxSystem`](@ref) (`src/zemax.jl`), the canonical entry point
for the Zemax-import pipeline.

Fields:
    geo::Vector{AbstractSurface}    - real, traceable optical surfaces; never includes the object surface
    objectSurface::AbstractSurface  - the object surface's own geometry (curvature/aperture, if any),
                                       kept separate from `geo` since it isn't part of the traceable
                                       refracting chain -- see `zemaxObjectToModelSurface`
    name::String                    - system name
    units::String                   - length units
    wavelengths::Vector{T}          - system wavelengths
    primaryWavelengthIndex::Int     - index into `wavelengths` marking the primary/reference wavelength
    objectDistance::T               - distance from the object to the first real surface; can be `Inf`
    objectAtInfinity::Bool          - true iff `objectDistance` is infinite -- stored explicitly
                                       (not just derived via `isinf`) so downstream ray-source code can
                                       branch on it directly
    apertureType::String            - which aperture specification is active: `"ENPD"`, `"OBNA"`,
                                       `"FNUM"`, or `"FLOA"` (float-by-stop, no numeric value)
    apertureValue::T                - the value for whichever `apertureType` is active (`NaN` for `"FLOA"`)
    fieldType::Int                  - field-type code (angle / object height / image height / ...)
    fields::Vector{Point2{T}}       - design field points, one `(x,y)` per field
    fieldWeight::Vector{T}          - per-field weight, index-aligned with `fields`
    glassCatalogs::Vector{String}   - glass catalog names materials should resolve against
    mode::String                    - `"SEQ"` or `"NSC"`
    notes::String                   - freeform system notes/description
"""
mutable struct OpticalSystem{T<:Real}
    geo::Vector{AbstractSurface}
    objectSurface::AbstractSurface
    name::String
    units::String
    wavelengths::Vector{T}
    primaryWavelengthIndex::Int
    objectDistance::T
    objectAtInfinity::Bool
    apertureType::String
    apertureValue::T
    fieldType::Int
    fields::Vector{Point2{T}}
    fieldWeight::Vector{T}
    glassCatalogs::Vector{String}
    mode::String
    notes::String
end

