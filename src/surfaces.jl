
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
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :yellow, 
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

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
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine2, 
    attributesSurfaces = attributesSurfaces, ydir::Union{Vec3, Nothing}=nothing)

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
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine2,
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

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
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    θ::Float64,
    semiDiam::Float64,
    coating::APorString
    ; color = :brown, 
    ydir::Union{Vec3, Nothing}=nothing)

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



function reflectOAConic(surfname::String,
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
    ;color = :aquamarine2, attributesSurfaces = attributeSurfaces)

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
    pointInPlane::Point3,
    planenormal::Vec3,
    ydir::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :blue, attributesSurfaces = attributesSurfaces
    ) = reflectOAConic(surfname,
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
        asphere::AbstractVector{Float64},
        semiDiam::Float64,
        coating::APorString
        ;color = :aquamarine,  
        attributesSurfaces = attributesSurfaces, 
        ydir::Union{Vec3, Nothing}=nothing)
    

"""
function refractAsphere(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    asphere::AbstractVector{Float64},
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine,  
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

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
        asphere::AbstractVector{Float64},
        semiDiam::Float64,
        coating::APorString
        ;color = :aquamarine,  
        attributesSurfaces = attributesSurfaces, 
        ydir::Union{Vec3, Nothing}=nothing)
    

"""
function refractEvenAsphere(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    asphere::AbstractVector{Float64},
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine,  
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

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
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    asphere::AbstractVector{Float64},
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine, 
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

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
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    c::Float64,
    ϵ::Float64,
    asphere::AbstractVector{Float64},
    semiDiam::Float64,
    coating::APorString
    ;color = :aquamarine, 
    attributesSurfaces = attributesSurfaces, 
    ydir::Union{Vec3, Nothing}=nothing)

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

function sag(x::Float64, y::Float64, s::NoProfile)
    0.
end

function gbRadius(aperture::SizeLens, profile::NoProfile)
    aperture.semiDiameter
end

function gbWidths(a::SizeLens, p::NoProfile)
    diam = 2a.semiDiameter
    SVector(diam, diam, 0.)
end

function surfNormal(r::Point3, s::NoProfile)
    Vec3(0., 0., 1.)
end

function deltaToSurf(r::Ray, p::NoProfile)
    x0, y0, z0 = r.base
    L,M,N = r.dir

    if N ≈ 0.
        Δ =  NaN #ray parallel to flat surface
    else
        Δ = -z0/N
    end
    Δ
end

function modFunc(ray::Ray, normal::Vec3, d::NoBendIndex)
    true, ray.dir, d.refIndexIn
end

#=
function surfAmpFunc(dirIn::Vec3, dirOut::Vec3, normal::Vec3, newRayBase::Point3, ndex::NoBendIndex,amp::NoAmpParam)
    identityAmpMats(), dirOut
end
=#

function referencePlane(surfname::String,
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :khaki3, 
    ydir::Union{Vec3, Nothing}=nothing)

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
    pointInPlane::Point3,
    planenormal::Vec3,
    rinIn::Float64,
    rinOut::Float64,
    semiDiam::Float64,
    coating::APorString
    ;color = :blue, 
    attributesSurfaces = attributesSurfaces,  #attributesSurfaces is a global dictionary
    ydir::Union{Vec3, Nothing}=nothing)

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
            -curv1, semiDiam, Base.compute_assumed_setting)
        ]
    end
    lens
end



"""
lensASinglet(base, dir, curv1, ϵ1, aphere1, curv2, ϵ2, aphere2, 
      thick, lambda, riFunc, semiDiam; order = "forward", lensname = "ASinglet")
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
      thick, lambda, riFunc, semiDiam; order = "forward", lensname = "ASinglet")
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
