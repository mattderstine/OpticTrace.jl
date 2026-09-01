
export refractConic, refractSphere, reflectConic, cDiffuser
export reflectOAConic, reflectOAP, refractAsphere, reflectAsphere
export referencePlane, planeMirror
export lensSinglet, lensASinglet

"""
    refractConic(
        surfname,
        base, dir,
        refractiveindex in, refractive index out,
        curvature,
        ϵ,
        semiDiameter,
        coating)
"""
function refractConic(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    rinOut::T,
    c::T,
    ϵ::T,
    semiDiam::T,
    coating::APorString
    ;color = :yellow, 
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal,ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileConic( c, ϵ),
        DielectricT(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
    refractSphere(
        surfname,
        base, dir,
        refractiveindex in, refractive index out,
        curvature,
        semiDiameter,
        coating)
"""
function refractSphere(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    rinOut::T,
    c::T,
    semiDiam::T,
    coating::APorString
    ;color = :aquamarine2, 
    attributesSurfaces = attributesSurfaces, ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    refractConic(surfname,
        pointInPlane,
        planenormal,
        rinIn,
        rinOut,
        c,
        1., # ϵ
        semiDiam,
        coating
        ;color, attributesSurfaces, ydir)
end



"""
    reflectConic(
        surfname,
        base, dir,
        refractiveindex in, refractive index out,
        curvature,
        ϵ,
        semiDiameter,
        coating)
"""
function reflectConic(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    rinOut::T,
    c::T,
    ϵ::T,
    semiDiam::T,
    coating::APorString
    ;color = :aquamarine2,
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal,ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileConic( c, ϵ),
        MirrorR(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
    cDiffuser(
        surfname,
        base, dir,
        refractiveindex in, refractive index out,
        half deflection angle,
        semiDiameter,
        coating)
"""
function cDiffuser(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    rinOut::T,
    θ::T,
    semiDiam::T,
    coating::APorString
    ; color = :brown, 
    ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        NoProfile(0.),
        CDiffuser(tan(θ), rinIn, rinOut),
        AmpParam(coating),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end



"""
    reflectOAConic(surfname::String,
        pointInPlane::Point3,
        planenormal::Vec3,
        ydir::Vec3,
        offset::Vec3,
        rinIn::Float64,
        rinOut::Float64,
        c::Float64,
        ϵ::Float64,
        semiDiam::Float64,
        coating::APorString
        ;color = :aquamarine2, attributesSurfaces = attributesSurfaces)

Build a reflecting off-axis conic `OptSurface` (`SurfProfileOAConic`,
`MirrorR`). Unlike its sibling refract/reflect constructors above,
`ydir` is a required positional argument here rather than an optional
keyword defaulting to `nothing` -- off-axis conics need an explicit y
direction since `offset` is defined relative to it. See `reflectOAP`
below for a convenience wrapper that computes `offset` from `c` for the
common off-axis-parabola case.
"""
function reflectOAConic(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    ydir::Vec3{T},
    offset::Vec3{T},
    rinIn::T,
    rinOut::T,
    c::T,
    ϵ::T,
    semiDiam::T,
    coating::APorString
    ;color = :aquamarine2, attributesSurfaces = attributesSurfaces) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileOAConic( c, ϵ, offset),
        MirrorR(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
    convention is that rotation around z axis, ray going -z direction reflects to +y dir

    reflectOAP(surfname::String,
        pointInPlane::Point3,   #base point in plane of surface
        planenormal::Vec3,     #normal to plane of surface
        ydir::Vec3, 
        rinIn::Float64,         #refractive index in
        rinOut::Float64,        #refractive index out
        c::Float64,            #curvature
        semiDiam::Float64,     #semi diameter of surface
        coating::APorString #coating or string key to attributesSurfaces dictionary
        ;color = :blue, attributesSurfaces = attributesSurfaces
        )
"""
reflectOAP(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    ydir::Vec3{T},
    rinIn::T,
    rinOut::T,
    c::T,
    semiDiam::T,
    coating::APorString
    ;color = :blue, attributesSurfaces = attributesSurfaces
    ) where T<:Real = reflectOAConic(surfname,
                        pointInPlane,
                        planenormal,
                        ydir,
                        Vec3(0, 1/c, -0.5/c),
                        rinIn,
                        rinOut,
                        c,
                        0.0,
                        semiDiam,
                        coating; color, attributesSurfaces)
#=

    refracting asphere
    reflecting asphere
    OA conic, refract & reflect - not tested  yet
    prism
    hologram/grating


=#
#=

Overloads for referencePlane

=#

"""
    refractAsphere(surfname::String,
        pointInPlane::Point3,
        planenormal::Vec3,
        rinIn::Float64,
        rinOut::Float64,
        c::Float64,
        ϵ::Float64,
        asphere::AbstractVector{T},
        semiDiam::Float64,
        coating::APorString
        ;color = :aquamarine,
        attributesSurfaces = attributesSurfaces,
        ydir::Union{Vec3, Nothing}=nothing)


"""
function refractAsphere(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    rinOut::T,
    c::T,
    ϵ::T,
    asphere::AbstractVector{T},
    semiDiam::T,
    coating::APorString
    ;color = :aquamarine,  
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileAsphere( c, ϵ, asphere),
        DielectricT(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
    refractEvenAsphere(surfname::String,
        pointInPlane::Point3,
        planenormal::Vec3,
        rinIn::Float64,
        rinOut::Float64,
        c::Float64,
        ϵ::Float64,
        asphere::AbstractVector{T},
        semiDiam::Float64,
        coating::APorString
        ;color = :aquamarine,
        attributesSurfaces = attributesSurfaces,
        ydir::Union{Vec3, Nothing}=nothing)


"""
function refractEvenAsphere(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    rinOut::T,
    c::T,
    ϵ::T,
    asphere::AbstractVector{T},
    semiDiam::T,
    coating::APorString
    ;color = :aquamarine,  
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileEvenAsphere( c, ϵ, asphere),
        DielectricT(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
    reflectAsphere
    like refractAsphere but reflects

"""
function reflectAsphere(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    rinOut::T,
    c::T,
    ϵ::T,
    asphere::AbstractVector{T},
    semiDiam::T,
    coating::APorString
    ;color = :aquamarine, 
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileAsphere( c, ϵ, asphere),
        MirrorR(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir, 
        color
        )
end

"""
    reflectEvenAsphere
    like refractEvenAsphere but reflects

"""
function reflectEvenAsphere(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    rinOut::T,
    c::T,
    ϵ::T,
    asphere::AbstractVector{T},
    semiDiam::T,
    coating::APorString
    ;color = :aquamarine, 
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        SurfProfileEvenAsphere( c, ϵ, asphere),
        MirrorR(rinIn, rinOut),
        #AmpParam(coating),
        getAmpParams(coating; attributesSurfaces),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir, 
        color
        )
end

"""
    sag(x, y, s::NoProfile)

Sag of a flat `NoProfile` plane: always `0.`, regardless of `x`/`y`.
See `sag(x, y, s::SurfProfileConic)`'s docstring (`src/tracing.jl`) for
the general x/y/s/return contract shared by every `sag` method.
"""
function sag(x::T, y::T, s::NoProfile) where T<:Real
    0.
end

"""
    sag(x, y, s::ParaxialProfile)

Sag of a flat `ParaxialProfile` plane: always `0.`, regardless of
`x`/`y` -- identical to `sag(x, y, ::NoProfile)`, since an ideal thin
lens has no real sag (all its optical effect is in `modFunc`). See
`sag(x, y, s::SurfProfileConic)`'s docstring (`src/tracing.jl`) for the
general x/y/s/return contract shared by every `sag` method.
"""
function sag(x::T, y::T, s::ParaxialProfile) where T<:Real
    0.
end

"""
    gbRadius(aperture::SizeLens{T}, profile::NoProfile) where T<:Real

Bounding radius for a flat (`NoProfile`) surface sized by a `SizeLens`:
its `semiDiameter`. Needed alongside `gbWidths(a::SizeLens,
p::NoProfile)` below because `referencePlane` (below) pairs
`SizeLens`/`NoProfile`, a combination `mesh_primitives.jl`'s
`(SizeLens, SurfProfileConic)` method doesn't cover. See
`gbRadius(aperture::SizeLens, profile::SurfProfileConic)`'s docstring
(`src/mesh_primitives.jl`) for the general contract shared by every
`gbRadius` method.
"""
function gbRadius(aperture::SizeLens{T}, profile::NoProfile) where T<:Real
    aperture.semiDiameter
end

"""
    gbWidths(a::SizeLens{T}, p::NoProfile) where T<:Real

Bounding-box widths for a flat (`NoProfile`) surface sized by a
`SizeLens`: `(2*semiDiameter, 2*semiDiameter, 0.)`. See
`gbRadius(aperture::SizeLens, profile::NoProfile)`'s docstring above
for why this method exists alongside the `mesh_primitives.jl`
`(SizeLens, SurfProfileConic)` method.
"""
function gbWidths(a::SizeLens{T}, p::NoProfile) where T<:Real
    diam = 2a.semiDiameter
    SVector(diam, diam, 0.)
end

#=
function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newRayBase::Point3, ndex::NoBendIndex,amp::NoAmpParam)
    identityAmpMats(), dirOut
end
=#

"""
    referencePlane(surfname::String,
        pointInPlane::Point3,
        planenormal::Vec3,
        rinIn::Float64,
        semiDiam::Float64,
        coating::APorString
        ;color = :khaki3, ydir::Union{Vec3, Nothing}=nothing)

Build a non-refracting, non-reflecting `OptSurface` (flat `NoProfile`,
`NoBendIndex(rinIn)` -- same index on both sides, `NoAmpParam`): a
reference/model plane used for object/image planes and other
characterization surfaces that shouldn't alter a ray, only mark a
position. See `planeMirror` below for the analogous reflecting (mirror)
reference surface.
"""
function referencePlane(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    semiDiam::T,
    coating::APorString
    ;color = :khaki3, 
    ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    OptSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        SizeLens(semiDiam),
        NoProfile( 0.),
        NoBendIndex(rinIn),
        NoAmpParam(coating),
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end


"""
    planeMirror is just a mirror with zero curvature and conic
    could use NoProfile to speed things up but then would need to
    potentially overload other functions
"""
function planeMirror(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    rinOut::T,
    semiDiam::T,
    coating::APorString
    ;color = :blue, 
    attributesSurfaces = attributesSurfaces,  #attributesSurfaces is a global dictionary
    ydir::Union{Vec3{T}, Nothing}=nothing) where T<:Real

    reflectConic(surfname,
        pointInPlane,
        planenormal,
        rinIn,
        rinOut,
        0.,
        0.,
        semiDiam,
        coating; color, attributesSurfaces, ydir)
end

"""
helper function for singlet lens

Builds a two-surface spherical singlet lens (`refractSphere` x2), in
either `order = "forward"` or `order = "reverse"`.
"""
function lensSinglet(base, dir, curv1, curv2, thick, lambda, riFunc, semiDiam; order = "forward", lensname = "Singlet", coating = "default")
    ri = riFunc(lambda)
    base1 = base + thick .* dir
    if (order != "reverse" )
        lens = [
        refractSphere("$(lensname)_1", base, dir, refIndexDefault, ri,
            curv1, semiDiam, coating),
        refractSphere("$(lensname)_2", base1, dir, ri,refIndexDefault,
            curv2, semiDiam, coating),

        ]
    else
        

        lens = [
        refractSphere("$(lensname)_2", base, dir, refIndexDefault, ri, 
            -curv2, semiDiam, coating),
        refractSphere("$(lensname)_1", base1, dir, ri, refIndexDefault,
            -curv1, semiDiam, coating)
        ]
    end
    lens
end



"""
lensASinglet(base, dir, curv1, ϵ1, aphere1, curv2, ϵ2, aphere2,
      thick, lambda, riFunc, semiDiam; order = "forward", lensname = "ASinglet", coating = "default")
"""
function lensASinglet(base, dir, curv1, ϵ1, aphere1, curv2, ϵ2, aphere2, thick, lambda, riFunc, semiDiam; order = "forward", lensname = "ASinglet", coating="default")
    ri = riFunc(lambda)
    # should check if the input is really an asphere. if not make the surface spherical #ToDo
    base1 = base + thick .* dir
    if (order !="reverse" && order != "forward")
        error("order must be 'forward' or 'reverse', got $order")
    end
    if (order != "reverse" )
        lens = [
        refractAsphere("$(lensname)_1", base, dir, refIndexDefault, ri,
            curv1, ϵ1, aphere1, semiDiam, coating),
        refractAsphere("$(lensname)_2", base1, dir, ri,refIndexDefault,
            curv2, ϵ2, aphere2, semiDiam, coating)
        ]
    else

        lens = [
        refractAsphere("$(lensname)_2", base, dir, refIndexDefault, ri, 
            -curv2, ϵ2, -aphere2, semiDiam, coating),
        refractAsphere("$(lensname)_1", base1, dir, ri, refIndexDefault, 
            -curv1, ϵ1, -aphere1, semiDiam, coating)
        ]
    end
    lens
end

"""
lensEASinglet(base, dir, curv1, ϵ1, aphere1, curv2, ϵ2, aphere2,
      thick, lambda, riFunc, semiDiam; order = "forward", lensname = "ASinglet", coating = "default")
"""
function lensEASinglet(base, dir, curv1, ϵ1, aphere1, curv2, ϵ2, aphere2, thick, lambda, riFunc, semiDiam; order = "forward", lensname = "ASinglet", coating="default")
    ri = riFunc(lambda)
    # should check if the input is really an asphere. if not make the surface spherical #ToDo
    base1 = base + thick .* dir
    if (order !="reverse" && order != "forward")
        error("order must be 'forward' or 'reverse', got $order")
    end
    if (order != "reverse" )
        lens = [
        refractEvenAsphere("$(lensname)_1", base, dir, refIndexDefault, ri,
            curv1, ϵ1, aphere1, semiDiam, coating),
        refractEvenAsphere("$(lensname)_2", base1, dir, ri,refIndexDefault,
            curv2, ϵ2, aphere2, semiDiam, coating)
        ]
    else

        lens = [
        refractEvenAsphere("$(lensname)_2", base, dir, refIndexDefault, ri, 
            -curv2, ϵ2, -aphere2, semiDiam, coating),
        refractEvenAsphere("$(lensname)_1", base1, dir, ri, refIndexDefault, 
            -curv1, ϵ1, -aphere1, semiDiam, coating)
        ]
    end
    lens
end
