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

    Despite the `s::AbstractSurface` signature, this only actually works
    for `OptSurface` in practice: it calls `inOrOut(s)`, which only has
    a method for `OptSurface`. Calling this on a `ModelSurface` (or any
    other `AbstractSurface` subtype) throws a `MethodError` (see
    `TODO.md`). `OptSurface` also has its own identical, redundant
    method below.
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
    GeometryBasics.normals(s::OptSurface, nvertices=60)

Identical in body to `GeometryBasics.normals(s::AbstractSurface,
nvertices=60)` above -- this more-specific `OptSurface` method shadows
the generic one for `OptSurface` but computes exactly the same thing,
so it's redundant (see `TODO.md`).
"""
function GeometryBasics.normals(s::OptSurface, nvertices=60)
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

Intended to return the bounding-box widths of a `Washer`'s outer
extent, as `SVector(2*semiDiameter, 2*semiDiameter, 0.)`. See
`gbWidths(a::SizeLens, p::SurfProfileConic)`'s docstring above for the
general contract shared by every `gbWidths` method.

**This method is broken**: it references `a.SemiDiameter`/
`a.SemiDiamater`, neither of which is a real field of `Washer` (the
actual field is `semiDiameter`, and `SemiDiamater` is also misspelled)
-- calling it throws a field-access error. See `TODO.md`.
"""
function gbWidths(a::Washer, p::NoProfile)
    SVector(a.SemiDiameter, a.SemiDiamater, 0.)
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

Intended to be `gbRadius(c, c.profile)` -- see
`gbRadius(aperture::Washer, profile::NoProfile)`'s docstring above.

**This method is broken**: it calls `gbRadius(c.semiDiameter,
c.profile)`, passing the bare `semiDiameter` number where `gbRadius`
expects the whole `Washer`. No `gbRadius` method matches a number, so
this throws a `MethodError`. Reachable from live plotting
(`plotModelSurf!`/`Makie.mesh!`) whenever a `ModelSurface`'s aperture is
a `RoundAperture` with a nonzero `obscure`, if Makie's mesh conversion
calls this method. See `TODO.md`.
"""
GeometryBasics.radius(c::Washer) = gbRadius(c.semiDiameter, c.profile) #fix

"""
    GeometryBasics.widths(c::Washer)

Intended to be `gbWidths(c, c.profile)` -- see
`gbWidths(a::Washer, p::NoProfile)`'s docstring above.

**This method is broken** the same way as `GeometryBasics.radius(c::Washer)`
above: it calls `gbWidths(c.semiDiameter, c.profile)`, passing a number
instead of the `Washer` itself, which throws a `MethodError` (on top of
the separate bug already in `gbWidths(a::Washer, p::NoProfile)`
itself). See `TODO.md`.
"""
GeometryBasics.widths(c::Washer) = gbWidths(c.semiDiameter, c.profile)

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

**This method is broken** in the *finite*-aperture branch (i.e.
whenever `wclear`/`lclear` are not both infinite -- the ordinary,
everyday case for a `RectAperture`): it references `a.clear`, which is
not a field of `RectAperture` (the actual fields are `wclear`/`lclear`)
-- calling it with a finite clear aperture throws a field-access error.
The `wclear == ∞ && lclear == ∞` branch, by contrast, works correctly.
Like `gbWidths(a::RectAperture, ...)` above, this currently has no live
call site, so the bug hasn't surfaced -- but it would affect the common
case, not just an edge case, if that changed. See `TODO.md`.
"""
function gbRadius(a::RectAperture, profile::NoProfile)
    if a.wclear == ∞ && a.lclear == ∞
        sqrt(a.wo^2+a.lo^2) #if only a rectangular obscuration, then size that
    else
        sqrt(a.wclear^2+a.clear^2) #can be infinity if one is finite
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

Intended to return the bounding-box widths of a `Disk`, as
`SVector(2*semiDiameter, 2*semiDiameter, 0.)`. See
`gbWidths(a::SizeLens, p::SurfProfileConic)`'s docstring above for the
general contract shared by every `gbWidths` method.

**This method has the same bug as `gbWidths(a::Washer, p::NoProfile)`**
above: it references the nonexistent fields `a.SemiDiameter`/
`a.SemiDiamater` instead of `a.semiDiameter`. See `TODO.md`.
"""
function gbWidths(a::Disk, p::NoProfile)
   SVector(a.SemiDiameter, a.SemiDiamater, 0.)
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

Intended to be `gbRadius(c, c.profile)` -- see
`gbRadius(aperture::Disk, profile::NoProfile)`'s docstring above.

**This method has the same bug as `GeometryBasics.radius(c::Washer)`**
above: it calls `gbRadius(c.semiDiameter, c.profile)`, passing a number
where `gbRadius` expects the whole `Disk`, which throws a
`MethodError`. See `TODO.md`.
"""
GeometryBasics.radius(c::Disk) = gbRadius(c.semiDiameter, c.profile) #fix

"""
    GeometryBasics.widths(c::Disk)

Intended to be `gbWidths(c, c.profile)` -- see
`gbWidths(a::Disk, p::NoProfile)`'s docstring above.

**This method has the same bug as `GeometryBasics.widths(c::Washer)`**
above: it calls `gbWidths(c.semiDiameter, c.profile)`, passing a number
instead of the `Disk` itself, on top of the separate bug already in
`gbWidths(a::Disk, p::NoProfile)`. See `TODO.md`.
"""
GeometryBasics.widths(c::Disk) = gbWidths(c.semiDiameter, c.profile)

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
