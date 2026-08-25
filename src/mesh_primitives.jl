##



#=
The functions needed to make OptSurf appear as an object compatiable with
GeometryBasics and thus suitable for plotting by Makie.


=#

#=

The functions needed to render optical surfaces as meshes. With luck fancier
solid models can be created using other packages.

=#

"""
    samplePoints - use the information in the aperture struct to define vertices
                    in the local x,y plane that will sample the Surface

    This needs to be overloaded for every decendant of AbstractSize
"""
function samplePoints(aperture::SizeLens{T},nvertices) where T<:Real
    r = LinRange(0., aperture.semiDiameter, nvertices);
    φ = LinRange(0., 2pi, nvertices)
    inner(r, φ) = [r*cos(φ), r*sin(φ)]
    return ivec((inner(r, φ) for r in r, φ in φ))
end

"""
gbWidths - provide widths for GeometryBasics
    returns SVector{3, float64}  of approx size
"""
function gbWidths(a::SizeLens{T}, p::SurfProfileConic{T}) where T<:Real
    diam = 2a.semiDiameter
    SVector(diam, diam, p.curv * a.semiDiameter^2)
end


"""
gbRadius - provide a radius for GeometryBasics
    OptSurface does not have enough info so the types of aperture & @profile
    are used to find the right formula
    returns the radius
"""
function gbRadius(aperture::SizeLens{T}, profile::SurfProfileConic{T}) where T<:Real
    aperture.semiDiameter
end


# Overloaded Functions

"""
    GeometryBasics.origin(c::OptSurface)

The origin of `c`, in global coordinates: `c.base.base`.
"""
GeometryBasics.origin(c::OptSurface) = c.base.base #the ORIGIN of the surface in global coordinates

"""
    GeometryBasics.radius(c::OptSurface)

Approximate bounding radius of `c`, in global coordinates:
`gbRadius(c.aperture, c.profile)`. See `gbRadius`'s docstring above for
the aperture/profile combinations it actually supports -- an
unsupported combination throws a `MethodError` here (see `TODO.md`).
"""
GeometryBasics.radius(c::OptSurface) = gbRadius(c.aperture, c.profile) #fix

"""
    GeometryBasics.widths(c::OptSurface)

Approximate bounding box widths of `c`, in global coordinates:
`gbWidths(c.aperture, c.profile)`. See `gbWidths`'s docstring above for
the aperture/profile combinations it actually supports -- an
unsupported combination throws a `MethodError` here (see `TODO.md`).
"""
GeometryBasics.widths(c::OptSurface) = gbWidths(c.aperture, c.profile)

"""
    coordinates - use sag & samplePoints to compute vertex locations.
    with any luck it can remain for all surfaces

"""
function GeometryBasics.coordinates(s::OptSurface, nvertices=60)
    a= samplePoints(s.aperture, nvertices)
    inner(t) = s.toGlobalCoord(Point3(t[1], t[2] ,sag(t[1], t[2], s.profile)))
    (inner(a) for a in a)
end

"""
    texturecoordinates - this is stolen directly from GeometryBasics for Sphere

"""
function GeometryBasics.texturecoordinates(s::AbstractSurface, nvertices=60)
    ux = LinRange(0, 1, nvertices)
    return ivec(((φ, θ) for θ in reverse(ux), φ in ux))
end

"""
    faces - this is stolen directly from GeometryBasics for Sphere

"""
function GeometryBasics.faces(s::AbstractSurface, nvertices=60)
    return GeometryBasics.faces(GeometryBasics.Rect(0, 0, 1, 1), (nvertices, nvertices))
end

"""
    normals - we need normals to do raytraces. Use them to do fancy rendering

    Works for any `AbstractSurface` subtype: `inOrOut(s)` has a method
    for every concrete surface type the package defines (`OptSurface`,
    `ModelSurface`).
"""
function GeometryBasics.normals(s::AbstractSurface, nvertices=60)
    a= samplePoints(s.aperture,nvertices)
    dir = inOrOut(s)
    #=
    if dir == -1
        println("Normals reverse on : $(s.surfname)")
    end
    =#
    inner(t) = dir .* s.toGlobalDir(surfNormal(Point3(t[1], t[2] ,sag(t[1], t[2], s.profile)),s.profile))
    (inner(a) for a in a)
end


"""
    inOrOut - tell shader which way the normal is pointing. Uses referactive
        index info to make that decision.
        returns 1 for out, -1 for in
"""
function inOrOut(s::OptSurface)

    mod = s.mod
    if mod.refIndexIn > mod.refIndexOut
        r = 1 #outputSpot
    else
        r = -1
    end

    r
end

"""
    inOrOut(s::ModelSurface)

`ModelSurface` has no refractive-index-in/out distinction to make (no
`mod` field, unlike `OptSurface`) -- always returns `1`, so its normals
point outward by convention. See `inOrOut(s::OptSurface)` above for the
refractive-index-based convention used for actual optical surfaces.
"""
function inOrOut(s::ModelSurface)
    1
end


#=
GeometryBasics.widths(c::OptSurf) = SVector(radius(c)*2,radius(c)*2,c.curv*radius(c)^2) #fix
GeometryBasics.radius(c::OptSurf) = c.semiDiam #fix


Base.minimum(c::OptSurf) = Vec{3, Float32}(origin(c)) - Vec{3, Float32}(radius(c), radius(c), 0.)
Base.maximum(c::OptSurf) = Vec{3, Float32}(origin(c)) + Vec{3, Float32}(radius(c), radius(c), c.curv*radius(c)^2)
=#
# overloaded methods for apertures

"""
    Washer{N,T} <: AbstractSurface{N,T}

A flat annular (washer-shaped) surface, used to render an obscuration
or clipping ring around a `ModelSurface` (see `plotModelSurf!` in
`src/plotting.jl`, which builds one from a `RoundAperture`'s `obscure`
size). Not a raytraceable optical surface -- it exists only to be
rendered by Makie via the `GeometryBasics` overloads below.

Fields:
- `base::Point{N}` -- origin, in global coordinates
- `dir::Vec{N}` -- local z direction, in global coordinates
- `semiDiameter::T` -- inner radius of the washer ring
- `toGlobalCoord::AffineMap`, `toGlobalDir::LinearMap` -- local->global
  transforms
- `profile::AbstractSurfProfile` -- always `NoProfile(0.)` in practice

A `Washer(base, dir, semiDiameter, toGlobalCoord, toGlobalDir)`
convenience constructor (below) builds one with `profile = NoProfile(0.)`
filled in automatically. (That constructor isn't separately documented
here: in Julia, a same-named outer constructor's docstring and the
struct's own docstring both bind to the same "no signature" doc slot,
so whichever is written second silently replaces the first -- adding
one here would have discarded this docstring rather than adding a
second, distinct entry.)
"""
struct Washer{N,T} <: AbstractSurface{N,T}
    base::Point{N}
    dir::Vec{N}
    semiDiameter::T
    toGlobalCoord::AffineMap
    toGlobalDir::LinearMap
    profile::AbstractSurfProfile  #NoProfile
end

"""
    samplePoints(aperture::Washer, nvertices, dsize=0.3)

Sample points on `aperture`'s annular ring, from its inner
`semiDiameter` out to `semiDiameter*(1+dsize)`. See
`samplePoints(aperture::SizeLens, nvertices)`'s docstring above for the
general contract shared by every `samplePoints` method.
"""
function samplePoints(aperture::Washer,nvertices, dsize = 0.3)
    r = LinRange(aperture.semiDiameter, aperture.semiDiameter*(1+dsize), nvertices)
    φ = LinRange(0., 2pi, nvertices)
    inner(r, φ) = [r*cos(φ), r*sin(φ)]
    return ivec((inner(r, φ) for r in r, φ in φ))
end

"""
    gbWidths(a::Washer, p::NoProfile)

Bounding-box widths of a `Washer`'s outer extent, as
`SVector(2*semiDiameter, 2*semiDiameter, 0.)`. See
`gbWidths(a::SizeLens, p::SurfProfileConic)`'s docstring above for the
general contract shared by every `gbWidths` method.
"""
function gbWidths(a::Washer, p::NoProfile)
    SVector(2a.semiDiameter, 2a.semiDiameter, 0.)
end

"""
    gbRadius(aperture::Washer, profile::NoProfile)

Bounding radius of a `Washer`: its outer `semiDiameter`. See
`gbRadius(aperture::SizeLens, profile::SurfProfileConic)`'s docstring
above for the general contract shared by every `gbRadius` method.
"""
function gbRadius(aperture::Washer, profile::NoProfile)
    aperture.semiDiameter
end


const Washer(base::Point{3}, dir::Vec{3},
    semiDiameter::Float64,toGlobalCoord::AffineMap,toGlobalDir::LinearMap) =
        Washer{3, Float64}(base, dir, semiDiameter,toGlobalCoord, toGlobalDir,NoProfile(0.))

"""
    GeometryBasics.origin(c::Washer)

The origin of `c`, in global coordinates: `c.base`.
"""
GeometryBasics.origin(c::Washer) = c.base #the ORIGIN of the surface in global coordinates

"""
    GeometryBasics.radius(c::Washer)

`gbRadius(c, c.profile)` -- see
`gbRadius(aperture::Washer, profile::NoProfile)`'s docstring above.
"""
GeometryBasics.radius(c::Washer) = gbRadius(c, c.profile)

"""
    GeometryBasics.widths(c::Washer)

`gbWidths(c, c.profile)` -- see
`gbWidths(a::Washer, p::NoProfile)`'s docstring above.
"""
GeometryBasics.widths(c::Washer) = gbWidths(c, c.profile)

"""
    GeometryBasics.coordinates(s::Washer, nvertices=60)

Vertex locations for rendering `s`'s flat annular ring: samples `s`'s
ring (via `samplePoints`) and maps each local `(x,y,0)` point to global
coordinates. See `GeometryBasics.coordinates(s::OptSurface,
nvertices=60)`'s docstring above for the general contract shared by
every `coordinates` method.
"""
function GeometryBasics.coordinates(s::Washer, nvertices=60)
    a= samplePoints(s, nvertices)
    inner(t) = s.toGlobalCoord(Point3(t[1], t[2] , 0.))
    (inner(a) for a in a)
end


"""
    GeometryBasics.normals(s::Washer, nvertices=60)

Surface normals for rendering `s`: since a `Washer` is flat, this is
the same constant global-space `ZAXIS` direction repeated once per
sample point (the samples themselves, from `samplePoints`, are computed
but only used to determine how many normals to yield).
"""
function GeometryBasics.normals(s::Washer, nvertices=60)
    a= samplePoints(s,nvertices)
    inner(t) = s.toGlobalDir(ZAXIS)
    (inner(a) for a in a)
end


"""
    gbWidths(a::RectAperture, p::NoProfile)

Bounding-box widths of a `RectAperture`: `2*wclear` x `2*lclear` (or
`2*wo`/`2*lo` on whichever axis has no finite clear aperture, i.e.
`wclear`/`lclear == ∞`). See `gbWidths(a::SizeLens,
p::SurfProfileConic)`'s docstring above for the general contract shared
by every `gbWidths` method.

No live code currently constructs an `OptSurface`/`ModelSurface` with a
`RectAperture` aperture and a `NoProfile` profile, so this method has
no confirmed call site (see `TODO.md`'s Tests needed section).
"""
function gbWidths(a::RectAperture, p::NoProfile)
    if a.wclear == ∞
        w = 2a.wo
    else
        w = 2a.wclear
    end
    if a.lclear == ∞
        l = 2a.lo
    else
        l = 2a.lclear
    end
    SVector(w, l, 0.)
end

"""
    gbRadius(a::RectAperture, profile::NoProfile)

Bounding radius of a `RectAperture`: `sqrt(wclear²+lclear²)` if both
have a finite clear aperture, else `sqrt(wo²+lo²)` (the obscuration
half-widths, used as a fallback size measure). See
`gbRadius(aperture::SizeLens, profile::SurfProfileConic)`'s docstring
above for the general contract shared by every `gbRadius` method.
"""
function gbRadius(a::RectAperture, profile::NoProfile)
    if a.wclear == ∞ && a.lclear == ∞
        sqrt(a.wo^2+a.lo^2) #if only a rectangular obscuration, then size that
    else
        sqrt(a.wclear^2+a.lclear^2) #can be infinity if one is finite
    end
end

#overloaded methods for apertures

"""
    Disk{N,T} <: AbstractSurface{N,T}

A flat circular disk surface, used to render a `ModelSurface`'s round
aperture/clipping region (see `plotModelSurf!` in `src/plotting.jl`,
which builds one from a `RoundAperture`). Not a raytraceable optical
surface -- exists only to be rendered by Makie via the `GeometryBasics`
overloads below. Field layout mirrors `Washer`'s (see that struct's
docstring above), except `samplePoints` covers the disk's full interior
(`0` to `semiDiameter`) rather than an annular ring.

A `Disk(base, dir, semiDiameter, toGlobalCoord, toGlobalDir)`
convenience constructor (below) builds one with `profile = NoProfile(0.)`
filled in automatically -- not separately documented, for the same
same-name-doc-slot reason given on `Washer`'s docstring above.
"""
struct Disk{N,T} <: AbstractSurface{N,T}
   base::Point{N}
   dir::Vec{N}
   semiDiameter::T
   toGlobalCoord::AffineMap
   toGlobalDir::LinearMap
   profile::AbstractSurfProfile  #NoProfile
end

"""
    samplePoints(aperture::Disk, nvertices)

Sample points across `aperture`'s full disk, from its center out to
`semiDiameter`. See `samplePoints(aperture::SizeLens, nvertices)`'s
docstring above for the general contract shared by every
`samplePoints` method.
"""
function samplePoints(aperture::Disk,nvertices)
   r = LinRange(0., aperture.semiDiameter, nvertices)
   φ = LinRange(0., 2pi, nvertices)
   inner(r, φ) = [r*cos(φ), r*sin(φ)]
   return ivec((inner(r, φ) for r in r, φ in φ))
end

"""
    gbWidths(a::Disk, p::NoProfile)

Bounding-box widths of a `Disk`, as `SVector(2*semiDiameter,
2*semiDiameter, 0.)`. See `gbWidths(a::SizeLens, p::SurfProfileConic)`'s
docstring above for the general contract shared by every `gbWidths`
method.
"""
function gbWidths(a::Disk, p::NoProfile)
   SVector(2a.semiDiameter, 2a.semiDiameter, 0.)
end

"""
    gbRadius(aperture::Disk, profile::NoProfile)

Bounding radius of a `Disk`: its `semiDiameter`. See
`gbRadius(aperture::SizeLens, profile::SurfProfileConic)`'s docstring
above for the general contract shared by every `gbRadius` method.
"""
function gbRadius(aperture::Disk, profile::NoProfile)
   aperture.semiDiameter
end




const Disk(base::Point3, dir::Vec3, semiDiameter::Float64,toGlobalCoord::AffineMap,toGlobalDir::LinearMap) = Disk(base, dir, semiDiameter, toGlobalCoord, toGlobalDir,NoProfile(0.))

"""
    GeometryBasics.origin(c::Disk)

The origin of `c`, in global coordinates: `c.base`.
"""
GeometryBasics.origin(c::Disk) = c.base #the ORIGIN of the surface in global coordinates

"""
    GeometryBasics.radius(c::Disk)

`gbRadius(c, c.profile)` -- see
`gbRadius(aperture::Disk, profile::NoProfile)`'s docstring above.
"""
GeometryBasics.radius(c::Disk) = gbRadius(c, c.profile)

"""
    GeometryBasics.widths(c::Disk)

`gbWidths(c, c.profile)` -- see
`gbWidths(a::Disk, p::NoProfile)`'s docstring above.
"""
GeometryBasics.widths(c::Disk) = gbWidths(c, c.profile)

"""
    GeometryBasics.coordinates(s::Disk, nvertices=60)

Vertex locations for rendering `s`'s flat disk: samples `s`'s interior
(via `samplePoints`) and maps each local `(x,y,0)` point to global
coordinates. See `GeometryBasics.coordinates(s::OptSurface,
nvertices=60)`'s docstring above for the general contract shared by
every `coordinates` method.
"""
function GeometryBasics.coordinates(s::Disk, nvertices=60)
   a= samplePoints(s, nvertices)
   inner(t) = s.toGlobalCoord(Point3(t[1], t[2] , 0.))
   (inner(a) for a in a)
end


"""
    GeometryBasics.normals(s::Disk, nvertices=60)

Surface normals for rendering `s`: since a `Disk` is flat, this is the
same constant global-space `ZAXIS` direction repeated once per sample
point. The author's own trailing comment on this method ("isn't this
the dumbest thing you've ever seen. Not taking time to fix it.") flags
computing and discarding `samplePoints`' actual coordinates just to get
a repeat count as wasteful -- see `GeometryBasics.normals(s::Washer,
nvertices=60)`'s docstring above for the identical pattern there.
"""
function GeometryBasics.normals(s::Disk, nvertices=60)
    a= samplePoints(s,nvertices)
    inner(t) = s.toGlobalDir(ZAXIS)
    (inner(a) for a in a) #isn't this the dumbest thing you've ever seen. Not taking time to fix it.
end
