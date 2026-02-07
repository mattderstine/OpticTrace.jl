
export  roundAperture, rectAperture, isAperture

function isAperture(a::AbstractSize{T}) where T<:Real
    false
end

function isAperture(a::TrueAperture{T}) where T<:Real
    true
end
"""
clipAperture(interceptPoint, aperture)
    interceptPoint - ray interecept in local coordinates
    aperture - array structured, with apeture and/or obscuration

    returns Boolean

    A single function has both obscuration and aperture. 

"""
function clipAperture(localBase::Point3{T}, aperture::RoundAperture{T}) where T<:Real
    #println("\nclipAperture\nlocalBase = $localBase  aperture=$aperture")
    r2 = localBase[1]^2 + localBase[2]^2 #z is supposed to be zero
    t = (r2 < aperture.obscure^2 )|| (r2 > aperture.semiDiameter^2 )
    #println("t = $t")
    t
end

function clipAperture(localBase::Point3{T}, aperture::RectAperture{T}) where T<:Real
    x = localBase[1] #width
    y = localBase[2] #length
    #println("x = $x  y = $y  $(aperture.wo)  $(aperture.wclear)  $(aperture.lo)  $(aperture.lclear)")
    t = (abs(x)<aperture.wo && abs(y)<aperture.lo) || abs(x)>aperture.wclear || abs(y)>aperture.lclear
    #println("test = $t")
    t
end

"""
define a round aperture

    roundAperture(surfname, pointInPlane, planenormal,rinIn,obscure, semiDiam
        ;color, ydir)

        surfname::String       - name to identify the surface
        pointInPlane::Point3   - point in the plane of the surface
        planenormal::Vec3      - normal vector of the surface
        rinIn::Float64         - refractive index of the surface
        obscure::Float64       - radius of an round obscuration
        semiDiam::Float64      - radius of the round aperture
        ;color = :blue,        - color of the surface
        ydir                   - direction of the y-axis (or nothing if not specified)

"""
function roundAperture(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    rinIn::T,
    obscure::T,
    semiDiam::T
    ;color = :blue, ydir::Union{Vec3{T}, Nothing} = nothing) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    ModelSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        RoundAperture(obscure, semiDiam),
        NoProfile(0.),
        rinIn,
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

"""
define a rectangular aperture
"""
function rectAperture(surfname::String,
    pointInPlane::Point3{T},
    planenormal::Vec3{T},
    ydir::Vec3{T},
    rinIn::T,
    obscurex::T,
    obscurey::T,
    semiDiamx::T,
    semiDiamy::T
    ;color = :blue) where T<:Real

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
                updateCoordChange(pointInPlane, planenormal, ydir)

    ModelSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        RectAperture(obscurex, obscurey, semiDiamx, semiDiamy),
        NoProfile(0.),
        rinIn,
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,
        color
        )
end

#=
"""
function radiusAperture(aperture::AbstractSize)
    return semiDiameter for RoundAperture
    return norm([semiDiamx, semiDimy]) for RectAperture
"""
function radiusAperture(aperture::AbstractSize)
    return aperture.semiDiameter
end

function radiusAperture(aperture::RectAperture)
    return norm([aperture.wclear, aperture.lclear])
end

=#