

export Ray, SurfBase, Trace, OptSurface, ModelSurface, ExtendedGeometry, SurfProfileOAConic, SizeLens, RoundAperture, RectAperture
export NoProfile, NoBendIndex, NoAmpParam, updateCoordChange, AbstractSurface, AbstractAmplitudeParam, DielectricT,MirrorR,CDiffuser, NoBendIndex
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
`SurfProfileAsphere`, `SurfProfileEvenAsphere`.
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

A toroidal surface profile, with independent curvatures along the local
x and y directions. **Currently only partially implemented**: a `sag`
method exists (`src/tracing.jl`) but is explicitly code-commented as
"likely incorrect", and no `deltaToSurf`/`surfNormal` method exists at
all -- toroidal surfaces cannot actually be raytraced yet (see
`TODO.md`).

Fields:
- `curvY::T` -- curvature along the local y direction
- `curvX::T` -- curvature along the local x direction
"""
mutable struct SurfProfileToroid{T} <: AbstractSurfProfile{T}
    curvY :: T
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

