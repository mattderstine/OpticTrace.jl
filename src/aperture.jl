
export  roundAperture, rectAperture, isAperture

function isAperture(a::AbstractSize)
    false
end

function isAperture(a::TrueAperture)
    true
end
"""
clipAperture(interceptPoint, aperture)
    interceptPoint - ray interecept in local coordinates
    aperture - array structured, with apeture and/or obscuration

    returns Boolean

    A single function has both obscuration and aperture. 

"""
function clipAperture(localBase::SVector{3,Float64}, aperture::RoundAperture)
    #println("\nclipAperture\nlocalBase = $localBase  aperture=$aperture")
    r2 = localBase[1]^2 + localBase[2]^2 #z is supposed to be zero
    t = (r2 < aperture.obscure^2 )|| (r2 > aperture.semiDiameter^2 )
    #println("t = $t")
    t
end

function clipAperture(localBase::SVector{3,Float64}, aperture::RectAperture)
    x = localBase[1] #width
    y = localBase[2] #length
    #println("x = $x  y = $y  $(aperture.wo)  $(aperture.wclear)  $(aperture.lo)  $(aperture.lclear)")
    t = (abs(x)<aperture.wo && abs(y)<aperture.lo) || abs(x)>aperture.wclear || abs(y)>aperture.lclear
    #println("test = $t")
    t
end

"""
define a round aperture

"""
function roundAperture(surfname::String,
    pointInPlane::SVector{3,Float64},
    planenormal::SVector{3,Float64},
    rinIn::Float64,
    obscure::Float64,
    semiDiam::Float64
    ;color = :blue, ydir::Union{SVector{3,Float64}, Nothing} = nothing)

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir=#,toGCMat, toLCMat =
                updateCoordChange(pointInPlane, planenormal, ydir)

    ModelSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        RoundAperture(obscure, semiDiam),
        NoProfile(0.),
        rinIn,
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir, #  toGCMat, toLCMat,
        color
        )
end

"""
define a rectangular aperture
"""
function rectAperture(surfname::String,
    pointInPlane::SVector{3,Float64},
    planenormal::SVector{3,Float64},
    ydir::SVector{3,Float64},
    rinIn::Float64,
    obscurex::Float64,
    obscurey::Float64,
    semiDiamx::Float64,
    semiDiamy::Float64
    ;color = :blue)

    ydir, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir = #,toGCMat, toLCMat =
                updateCoordChange(pointInPlane, planenormal, ydir)

    ModelSurface(surfname,
        SurfBase(pointInPlane, planenormal, ydir),
        RectAperture(obscurex, obscurey, semiDiamx, semiDiamy),
        NoProfile(0.),
        rinIn,
        toGlobalCoord,toLocalCoord,toGlobalDir,toLocalDir,#toGCMat, toLCMat,
        color
        )
end
