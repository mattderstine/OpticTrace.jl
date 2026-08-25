#extended geometry functions
export findPerpenMap, updateCoordChange!
"""
    findPerpenMap(planenormal)
Find a perpendicular to plane normal - "local x direction" by using y dir as the
guess. 
"""
function findPerpenMap(planenormal::Vec3)
    tryy = YAXIS
    newx = cross(tryy, planenormal) # if planenormal == ZAXIS this give x axis
    if norm(newx)==0. # planenormal is along y axis
        tryy = XAXIS
        newx = cross(tryy, planenormal)
    end

    newx=normalize(newx)
    newy = Vec3(cross(planenormal, newx)...)
    newy, LinearMap(SMatrix{3,3}([newx; newy; planenormal]))
end

"""
    findPerpenMap(planenormal, ydir)
Find a perpendicular in the plane to planenormal and ydir
"""
function findPerpenMap(planenormal::Vec3, ydir::Union{Nothing, Vec3})
    if isnothing(ydir)
        tryy = YAXIS
        newx = cross(tryy, planenormal) # if planenormal == ZAXIS this give x axis
        if norm(newx) < 0.1 # planenormal is very close to Y axis
            tryy = XAXIS
            newx = cross(tryy, planenormal)
        end
        newx = normalize(newx)
        newy = Vec3(cross(planenormal, newx)...)
    else
        newx = cross(ydir, planenormal) # if planenormal == ZAXIS this give x axis
        newx = normalize(newx)
        newy = ydir
    end
    newy, LinearMap(SMatrix{3,3}([newx; newy; planenormal]))
end



"""
    updateCoordChange(pointInPlane::Point3, planenormal::Vec3)

Compute the local/global coordinate-transform functions for a plane at
`pointInPlane` with normal `planenormal`, choosing the local y direction
via `findPerpenMap(planenormal)` (no explicit `ydir` guess). Does not
mutate anything -- for updating a surface's stored coordinate change
in place, see [`updateCoordChange!`](@ref).

Returns `(ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir)`.
"""
function updateCoordChange(pointInPlane::Point3,
        planenormal::Vec3)
    pip::SVector{3,Float64} = pointInPlane
    ydir, toGlobalDir = findPerpenMap(planenormal)
    toGlobalCoord = compose(Translation(pip), toGlobalDir)
    toLocalCoord = inv(toGlobalCoord)
    toLocalDir = inv(toGlobalDir)
    return ydir,toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir
end

"""
    updateCoordChange(pointInPlane::Point3, planenormal::Vec3, ydir::Union{Vec3,Nothing})

Like `updateCoordChange(pointInPlane, planenormal)`, but the local y
direction is resolved via `findPerpenMap(planenormal, ydir)` -- pass an
explicit `ydir`, or `nothing` to fall back to the same default guess as
the 2-argument method. Does not mutate anything.

Returns `(ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir)`.
"""
function updateCoordChange(pointInPlane::Point3,
        planenormal::Vec3, ydir::Union{Vec3,Nothing})
    pip::SVector{3,Float64} = pointInPlane
    myydir, toGlobalDir = findPerpenMap(planenormal, ydir)
    toGlobalCoord = compose(Translation(pip), toGlobalDir)
    toLocalCoord = inv(toGlobalCoord)
    toLocalDir = inv(toGlobalDir)
    return myydir,toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir
end

"""
    updateCoordChange!(surf::T) where T<:OpticTrace.AbstractSurface

Recompute `surf`'s coordinate-transform fields (`toGlobalCoord`,
`toLocalCoord`, `toGlobalDir`, `toLocalDir`, and `base.ydir`) in place
from its current `base.base`/`base.dir`/`base.ydir`, via
[`updateCoordChange`](@ref). Use this after mutating a surface's base
point or direction directly.
"""
function updateCoordChange!(surf::T) where T<:OpticTrace.AbstractSurface
    surf.base.ydir, surf.toGlobalCoord, surf.toLocalCoord, surf.toGlobalDir, surf.toLocalDir = updateCoordChange(surf.base.base, surf.base.dir, surf.base.ydir)
end

"""
    updateCoordChange!(surf::T, newbase::SurfBase, newydir::Union{Vec3,Nothing}=nothing) where T<:OpticTrace.AbstractSurface

Replace `surf.base` with `newbase`, then recompute `surf`'s
coordinate-transform fields in place (via [`updateCoordChange`](@ref))
using `newbase`'s point/direction and the given `newydir` (or `nothing`
to use the default y-direction guess).
"""
function updateCoordChange!(surf::T, newbase::SurfBase, newydir::Union{Vec3,Nothing}=nothing) where T<:OpticTrace.AbstractSurface
    surf.base = newbase
    surf.base.ydir, surf.toGlobalCoord, surf.toLocalCoord, surf.toGlobalDir, surf.toLocalDir = updateCoordChange(surf.base.base, surf.base.dir, newydir)
end


"""
    EGeo(func, baseObject, size, wavelength; dir = ZAXIS, setup = defaultSetupGeo, parameters = Dict{Any,Any}() )
    simple way to set up ExtendedGeometry

    func - function defining geometry
    baseObject - location of object
    size - size of object
    wavelength - array of wavelengths in microns

    dir - direction of object
    setup - function to set parameters for func
    parameters - dictionary for parameters used by funcSetup and funcGeo

"""
function EGeo(func, baseObject, size, wavelength; dir = ZAXIS, setup = defaultSetupGeo, parameters = Dict{Any,Any}() )
    object = referencePlane("object", baseObject, dir, refIndexDefault, size, "none")
    newgeo = setup(func, object, wavelength, parameters)
    a = ExtendedGeometry(newgeo, func, setup, object, wavelength, parameters)
end


"""
    defaultSetupGeo(func, object, wavelength, parameters; numWL = 1)
    predefined function for simplest setup: nothing but running func

"""
function defaultSetupGeo(func, object, wavelength, parameters; numWL = 1)
    # perform setup functions for the GeometryBasics
    geo = func(parameters, wavelength[numWL])
end


"""
    updateEGeo!(egeo; numWL=1) calls the function necessary to create a static geometry for tracing
        egeo    extended geometry to update
        numWL   which wavelength to use

**This method is broken**: despite its name/docstring implying it
mutates `egeo`, it only calls `defaultSetupGeo(...)` and returns the
result -- it never assigns back into `egeo.geo`. Likely fix:
`egeo.geo = defaultSetupGeo(...)`. See `TODO.md`.
"""
function updateEGeo!(egeo::ExtendedGeometry; numWL = 1)
    defaultSetupGeo(egeo.funcGeo, egeo.surfaceObject, egeo.wavelength, egeo.parameters; numWL = numWL)
end


#=
    ampltude

=#
