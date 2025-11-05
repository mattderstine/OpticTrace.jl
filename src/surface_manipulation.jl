
#reverse a geo

export reverseGeo

function reverseGeo(geo::Vector{AbstractSurface})
    g = deepcopy(geo)
    if geo[1].base.dir != geo[end].base.dir
        error("Can only reverse geometry if first and last surface directions are the same")
    end
    beginpoint = geo[1].base.base
    endpoint = geo[end].base.base
    lengthgeo = norm(endpoint - beginpoint)

    for s in g
        reverseSurface!(s, beginpoint, lengthgeo)
    end
    g
end

function reverseSurface!(surf::OpticTrace.OptSurface, bpoint, epoint)
    reverseBase!(surf.base, bpoint, epoint)
    # Reverse the surface profile if needed
    reverseProfile!(surf.profile)
    reverseMod!(surf.mod)
    -, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
        updateCoordChange(surf.base.base, surf.base.dir,surf.base.ydir)
    surf.toGlobalCoord = toGlobalCoord
    surf.toLocalCoord = toLocalCoord
    surf.toGlobalDir = toGlobalDir
    surf.toLocalDir = toLocalDir

    return surf
end

function reverseBase!(base::SurfBase, bpoint, epoint)
    base.dir = -base.dir
    base.ydir = -base.ydir
    base.base = bpoint + epoint - base.base
    return base
end

function reverseProfile!(profile::SurfProfileSphere)
    # Implement profile-specific reversal logic if needed
    profile.curve = -profile.curve
    return profile
end

function reverseProfile!(profile::SurfProfileConic)
    # Implement profile-specific reversal logic if needed
    profile.curv = -profile.curv
    return profile
end

function reverseProfile!(profile::T) where T<:AbstractSurfProfile
    # Implement profile-specific reversal logic if needed
    profile.curv = -profile.curv
    profile.a = -profile.a
    return profile
end

function reverseProfile!(profile::SurfProfileToroid)
    # Implement profile-specific reversal logic if needed
    profile.curveX = -profile.curveX
    profile.curveY = -profile.curveY
    return profile
end

function reverseProfile!(profile::SurfProfileCyl)
    # Implement profile-specific reversal logic if needed
    profile.curv = -profile.curv
    profile.a = -profile.a
    return profile
end

function reverseMod!(mod::AbstractBendType)
    mod.refIndexIn, mod.refIndexOut = mod.refIndexOut, mod.refIndexIn
    return mod
end