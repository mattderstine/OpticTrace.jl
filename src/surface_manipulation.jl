
export reverseGeo, thickGeo

#reverse a geo



"""
    reverseGeo(geo::Vector{T}) where T<:AbstractSurface

Build a reversed copy of `geo`: reverses the surface order (via
`reverse!(deepcopy(geo))`) and, for each surface, flips its position
about the geometry's own start/end points (`reverseSurface!`) so the
geometry can be traced "backwards" (e.g. object<->image swapped).
Requires `geo[1].base.dir == geo[end].base.dir` (errors otherwise).

Returns the reversed `Vector{T}` (a deep copy; `geo` itself is
untouched other than being deep-copied, not shared).
"""
function reverseGeo(geo::Vector{T}) where T<:AbstractSurface

    if geo[1].base.dir != geo[end].base.dir
        error("Can only reverse geometry if first and last surface directions are the same")
    end
    g = reverse!(deepcopy(geo))
    beginpoint = geo[1].base.base
    endpoint = geo[end].base.base
    lengthgeo = norm(endpoint - beginpoint)

    for s in g
        reverseSurface!(s, beginpoint, endpoint)
    end
    g
end

"""
    reverseSurface!(surf::OptSurface, bpoint, epoint)

Mutate `surf` in place for `reverseGeo`: reflects its base position
about the midpoint of `bpoint`/`epoint` (`reverseBase!`), negates its
profile's curvature (`reverseProfile!`), swaps its entering/exiting
refractive indices (`reverseMod!`), and recomputes its local<->global
coordinate transforms accordingly. Returns `surf`.

Note: `surf.base.dir`/`surf.base.ydir` (the surface's local
orientation) are left unchanged -- only its position and profile/index
data are flipped. Whether that's correct for reversing propagation
direction through the surface, or a gap, isn't obvious from the code
alone; worth confirming before relying on this.
"""
function reverseSurface!(surf::OpticTrace.OptSurface, bpoint, epoint)
    #println("reverseSurface!:  bpoint: $bpoint epoint: $epoint")
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

"""
    reverseSurface!(surf::ModelSurface, bpoint, epoint)

Like `reverseSurface!(surf::OptSurface, bpoint, epoint)` above, but for
a non-refracting `ModelSurface`: reflects its base position
(`reverseBase!`) and negates its profile's curvature
(`reverseProfile!`), then recomputes its coordinate transforms.
`ModelSurface` has no `mod`/coating (unlike `OptSurface`) -- just a
fixed `refIndex` for OPD bookkeeping, not an in/out pair -- so there's
no `reverseMod!` call here, and `refIndex` is left unchanged. Returns
`surf`.

Note: like the `OptSurface` method above, `surf.base.dir`/
`surf.base.ydir` are left unchanged -- see that method's docstring for
the same open question about whether that's correct.
"""
function reverseSurface!(surf::OpticTrace.ModelSurface, bpoint, epoint)
    reverseBase!(surf.base, bpoint, epoint)
    reverseProfile!(surf.profile)
    -, toGlobalCoord, toLocalCoord, toGlobalDir, toLocalDir =
        updateCoordChange(surf.base.base, surf.base.dir, surf.base.ydir)
    surf.toGlobalCoord = toGlobalCoord
    surf.toLocalCoord = toLocalCoord
    surf.toGlobalDir = toGlobalDir
    surf.toLocalDir = toLocalDir

    return surf
end

"""
    reverseBase!(base::SurfBase, bpoint, epoint)

Reflect `base.base` about the midpoint of `bpoint` and `epoint`:
`base.base = bpoint + epoint - base.base`. Used by `reverseSurface!` to
reposition a surface when reversing a geometry end-for-end. Returns
`base`. Does not touch `base.dir`/`base.ydir`.
"""
function reverseBase!(base::SurfBase, bpoint, epoint)
    #println("epoint: $epoint bpoint: $bpoint base before: $(base.base)")
    base.base = bpoint + epoint - base.base
    #println("base after: $(base.base)")
    return base
end

"""
    reverseProfile!(profile::SurfProfileSphere)

Negate `profile`'s curvature in place (matching its sibling methods
below), and return `profile`.
"""
function reverseProfile!(profile::SurfProfileSphere)
    profile.curv = -profile.curv
    return profile
end

"""
    reverseProfile!(profile::SurfProfileConic)

Negate `profile.curv` in place (leaving `profile.ϵ` unchanged) and
return `profile`. See `reverseProfile!(profile::T) where
T<:AbstractSurfProfile` below for the fallback used by profile types
without their own specific method.
"""
function reverseProfile!(profile::SurfProfileConic)
    # Implement profile-specific reversal logic if needed
    profile.curv = -profile.curv
    return profile
end

"""
    reverseProfile!(profile::SurfProfileOAConic)

Negate `profile.curv` in place (leaving `profile.ϵ` unchanged) and
return `profile`. See `reverseProfile!(profile::T) where
T<:AbstractSurfProfile` below for the fallback used by profile types
without their own specific method.

**Not finished**: does not touch `profile.offset`, so the reversed
profile's off-axis offset is left describing the pre-reversal geometry.
See `TODO.md`.
"""
function reverseProfile!(profile::SurfProfileOAConic)
    profile.curv = -profile.curv
    println("reverseProfile!(profile::SurfProfileOAConic) is not finished because it does not handle the offset field properly.")
    return profile
end

"""
    reverseProfile!(profile::T) where T<:AbstractSurfProfile

Fallback for any `AbstractSurfProfile` subtype without its own more
specific `reverseProfile!` method above/below: negates `profile.curv`
and `profile.a` in place, and returns `profile`. Assumes every such
type has an `a` field -- true for `SurfProfileAsphere`/
`SurfProfileEvenAsphere`, the only concrete profile types that actually
reach this fallback (every other type, including `SurfProfileOAConic`
and `NoProfile`, has its own more specific method above/below that
takes precedence).
"""
function reverseProfile!(profile::T) where T<:AbstractSurfProfile
    # Implement profile-specific reversal logic if needed
    profile.curv = -profile.curv
    profile.a = -profile.a
    return profile
end

"""
    reverseProfile!(profile::NoProfile)

No-op: `NoProfile`'s only field (`curv`) is documented as never
actually used, so there's nothing meaningful to reverse. Returns
`profile` unchanged, for consistency with the other `reverseProfile!`
methods' return-the-profile contract.
"""
function reverseProfile!(profile::NoProfile)
    profile
end

"""
    reverseProfile!(profile::SurfProfileToroid)

Negate `profile`'s two curvatures in place (matching its sibling
methods above/below), and return `profile`.
"""
function reverseProfile!(profile::SurfProfileToroid)
    profile.curvX = -profile.curvX
    profile.curvY = -profile.curvY
    return profile
end

"""
    reverseProfile!(profile::SurfProfileCyl)

Negate `profile.curv` and `profile.a` in place (leaving `profile.ϵ`
unchanged) and return `profile`. See `reverseProfile!(profile::T)
where T<:AbstractSurfProfile` above for the shared `curv`/`a` negation
pattern used by profile types that don't need their own specific
method.
"""
function reverseProfile!(profile::SurfProfileCyl)
    # Implement profile-specific reversal logic if needed
    profile.curv = -profile.curv
    profile.a = -profile.a
    return profile
end

"""
    reverseMod!(mod::T) where T<:AbstractBendType

Swap `mod.refIndexIn`/`mod.refIndexOut` in place (every concrete
`AbstractBendType` -- `DielectricT`, `MirrorR`, `CDiffuser`,
`NoBendIndex` -- has these two fields), and return `mod`. Used by
`reverseSurface!` so a surface's entering/exiting index convention
matches the reversed propagation direction.
"""
function reverseMod!(mod::T) where T<:AbstractBendType
    #println("Reversing bend type: $(typeof(mod))")
    #println("Before reversal: rinIn=$(mod.refIndexIn) rinOut=$(mod.refIndexOut)")
    mod.refIndexIn, mod.refIndexOut = mod.refIndexOut, mod.refIndexIn
    #println("After reversal: rinIn=$(mod.refIndexIn) rinOut=$(mod.refIndexOut)")
    return mod
end

"""
    thickGeo(geo)

Compute the overall length of `geo`, projected onto its own initial
direction: `(geo[end].base.base - geo[begin].base.base) ⋅
geo[begin].base.dir` (dotting with the direction, rather than a plain
norm, so the result is still meaningful if the optics have been shifted
off-axis). Requires `geo[begin].base.dir == geo[end].base.dir` (errors
otherwise).
"""
function thickGeo(geo)
    if geo[begin].base.dir != geo[end].base.dir
        error("not in same direction")
    end
    v = geo[end].base.base - geo[begin].base.base
    return v⋅geo[begin].base.dir #project onto direction in case of shifted optics
end

