
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

Despite the `T<:AbstractSurface` signature, only works for geometries
made entirely of `OptSurface`s: `reverseSurface!` (below) has no method
for `ModelSurface` (e.g. surfaces built by `roundAperture`/
`rectAperture`), so a `geo` containing one throws a `MethodError`. See
`TODO.md`.
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

Intended to negate `profile`'s curvature in place (matching its
sibling methods below), and return `profile`.

**This method is broken**: it assigns `profile.curve` (note the extra
`e`), which is not a real field of `SurfProfileSphere` (the actual
field is `curv`) -- calling it throws a field-access error. See
`TODO.md`.
"""
function reverseProfile!(profile::SurfProfileSphere)
    # Implement profile-specific reversal logic if needed
    profile.curve = -profile.curve
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
    reverseProfile!(profile::T) where T<:AbstractSurfProfile

Fallback for any `AbstractSurfProfile` subtype without its own more
specific `reverseProfile!` method above/below: negates `profile.curv`
and `profile.a` in place, and returns `profile`.

**This method is broken for some profile types it applies to**: it
assumes every such type has an `a` field, which isn't true for
`SurfProfileOAConic` (fields: `curv`, `ϵ`, `offset`) or `NoProfile`
(field: `curv` only) -- calling this method on either throws a
field-access error. In practice this fallback is only actually reached
for `SurfProfileAsphere`/`SurfProfileEvenAsphere` (which do have `a`,
so those work), `SurfProfileOAConic`, and `NoProfile` -- every other
concrete profile type has its own more specific method above/below
that takes precedence. See `TODO.md`.
"""
function reverseProfile!(profile::T) where T<:AbstractSurfProfile
    # Implement profile-specific reversal logic if needed
    profile.curv = -profile.curv
    profile.a = -profile.a
    return profile
end

"""
    reverseProfile!(profile::SurfProfileToroid)

Intended to negate `profile`'s two curvatures in place (matching its
sibling methods above/below), and return `profile`.

**This method is broken**: it assigns `profile.curveX`/`profile.curveY`,
neither of which is a real field of `SurfProfileToroid` (the actual
fields are `curvX`/`curvY`, without the extra `e`) -- calling it throws
a field-access error. See `TODO.md`.
"""
function reverseProfile!(profile::SurfProfileToroid)
    # Implement profile-specific reversal logic if needed
    profile.curveX = -profile.curveX
    profile.curveY = -profile.curveY
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

