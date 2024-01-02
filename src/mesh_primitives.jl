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
function samplePoints(aperture::SizeLens,nvertices)
    r = LinRange(0., aperture.semiDiameter, nvertices);
    φ = LinRange(0., 2pi, nvertices)
    inner(r, φ) = [r*cos(φ), r*sin(φ)]
    return ivec((inner(r, φ) for r in r, φ in φ))
end

"""
gbWidths - provide widths for GeometryBasics
    returns SVector{3, float64}  of approx size
"""
function gbWidths(a::SizeLens, p::SurfProfileConic)
    diam = 2a.semiDiameter
    SVector(diam, diam, p.curv * a.semiDiameter^2)
end


"""
gbRadius - provide a radius for GeometryBasics
    OptSurface does not have enough info so the types of aperture & @profile
    are used to find the right formula
    returns the radius
"""
function gbRadius(aperture::SizeLens, profile::SurfProfileConic)
    aperture.semiDiameter
end


# Overloaded Functions

GeometryBasics.origin(c::OptSurface) = c.base.base #the ORIGIN of the surface in global coordinates
GeometryBasics.radius(c::OptSurface) = gbRadius(c.aperture, c.profile) #fix
GeometryBasics.widths(c::OptSurface) = gbWidths(c.aperture, c.profile)

"""
    coordinates - use sag & samplePoints to compute vertex locations.
    with any luck it can remain for all surfaces

"""
function GeometryBasics.coordinates(s::OptSurface, nvertices=60)
    ap= samplePoints(s.aperture, nvertices)
    inner(t) = s.toGlobalCoord(SVector{3}(t[1], t[2] ,sag(t[1], t[2], s.profile)))
    (inner(a) for a in ap)
end

"""
    coordinates - this is stolen directly from GeometryBasics for Sphere

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

"""
function GeometryBasics.normals(s::AbstractSurface, nvertices=60)
    ap= samplePoints(s.aperture,nvertices)
    dir = inOrOut(s)
    #=
    if dir == -1
        println("Normals reverse on : $(s.surfname)")
    end
    =#
    inner(t) = dir .* s.toGlobalDir(surfNormal(SVector{3}(t[1], t[2] ,sag(t[1], t[2], s.profile)),s.profile))
    (inner(a) for a in ap)
end

function GeometryBasics.normals(s::OptSurface, nvertices=60)
    ap= samplePoints(s.aperture,nvertices)
    dir = inOrOut(s)
    #=
    if dir == -1
        println("Normals reverse on : $(s.surfname)")
    end
    =#
    inner(t) = dir .* s.toGlobalDir(surfNormal(SVector{3}(t[1], t[2] ,sag(t[1], t[2], s.profile)),s.profile))
    (inner(a) for a in ap)
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
struct Washer{N,T} <: AbstractSurface{N,T}
    base::SVector{N, T}
    dir::SVector{N,T}
    semiDiameter::T
    toGlobalCoord::AffineMap
    toGlobalDir::LinearMap
    profile::AbstractSurfProfile  #NoProfile
end

function samplePoints(aperture::Washer,nvertices, dsize = 0.3)
    r = LinRange(aperture.semiDiameter, aperture.semiDiameter*(1+dsize), nvertices)
    φ = LinRange(0., 2pi, nvertices)
    inner(r, φ) = [r*cos(φ), r*sin(φ)]
    return ivec((inner(r, φ) for r in r, φ in φ))
end

function gbWidths(a::Washer, p::NoProfile)
    SVector(a.SemiDiameter, a.SemiDiamater, 0.)
end

function gbRadius(aperture::Washer, profile::NoProfile)
    aperture.semiDiameter
end


const Washer(base::SVector{3, Float64}, dir::SVector{3, Float64},
    semiDiameter::Float64,toGlobalCoord::AffineMap,toGlobalDir::LinearMap) =
        Washer{3, Float64}(base, dir, semiDiameter,toGlobalCoord, toGlobalDir,NoProfile(0.))

GeometryBasics.origin(c::Washer) = c.base #the ORIGIN of the surface in global coordinates
GeometryBasics.radius(c::Washer) = gbRadius(c.semiDiameter, c.profile) #fix
GeometryBasics.widths(c::Washer) = gbWidths(c.semiDiameter, c.profile)

function GeometryBasics.coordinates(s::Washer, nvertices=60)
    a= samplePoints(s, nvertices)
    inner(t) = s.toGlobalCoord(SVector{3}(t[1], t[2] , 0.))
    (inner(a) for a in a)
end


function GeometryBasics.normals(s::Washer, nvertices=60)
    a= samplePoints(s,nvertices)
    inner(t) = s.toGlobalDir(ZAXIS)
    (inner(a) for a in a)
end


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

function gbRadius(a::RectAperture, profile::NoProfile)
    if a.wclear == ∞ && a.lclear == ∞
        sqrt(a.wo^2+a.lo^2) #if only a rectangular obscuration, then size that
    else
        sqrt(a.wclear^2+a.clear^2) #can be infinity if one is finite
    end
end

#overloaded methods for apertures
struct Disk{N,T} <: AbstractSurface{N,T}
   base::SVector{N,T}
   dir::SVector{N,T}
   semiDiameter::T
   toGlobalCoord::AffineMap
   toGlobalDir::LinearMap
   profile::AbstractSurfProfile  #NoProfile
end

function samplePoints(aperture::Disk,nvertices)
   r = LinRange(0., aperture.semiDiameter, nvertices)
   φ = LinRange(0., 2pi, nvertices)
   inner(r, φ) = [r*cos(φ), r*sin(φ)]
   return ivec((inner(r, φ) for r in r, φ in φ))
end

function gbWidths(a::Disk, p::NoProfile)
   SVector(a.SemiDiameter, a.SemiDiamater, 0.)
end

function gbRadius(aperture::Disk, profile::NoProfile)
   aperture.semiDiameter
end




const Disk(base::SVector{3, Float64}, dir::SVector{3, Float64}, semiDiameter::Float64,toGlobalCoord::AffineMap,toGlobalDir::LinearMap) = Disk(base, dir, semiDiameter, toGlobalCoord, toGlobalDir,NoProfile(0.))

GeometryBasics.origin(c::Disk) = c.base #the ORIGIN of the surface in global coordinates
GeometryBasics.radius(c::Disk) = gbRadius(c.semiDiameter, c.profile) #fix
GeometryBasics.widths(c::Disk) = gbWidths(c.semiDiameter, c.profile)

function GeometryBasics.coordinates(s::Disk, nvertices=60)
   a= samplePoints(s, nvertices)
   inner(t) = s.toGlobalCoord(SVector{3}(t[1], t[2] , 0.))
   (inner(a) for a in a)
end


function GeometryBasics.normals(s::Disk, nvertices=60)
    a= samplePoints(s,nvertices)
    inner(t) = s.toGlobalDir(ZAXIS)
    (inner(a) for a in a) #isn't this the dumbest thing you've ever seen. Not taking time to fix it.
end
