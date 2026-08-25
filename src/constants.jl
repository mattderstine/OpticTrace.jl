
export ∞, ORIGIN, ZAXIS, YAXIS, XAXIS, refIndexDefault, rInDef, setRInDef

"""
    ORIGIN

The global-coordinate-system origin, `Point3(0., 0., 0.)`.
"""
const ORIGIN = Point3(0., 0., 0.)

"""
    ZAXIS, YAXIS, XAXIS

The global-coordinate-system basis vectors, as `Vec3`s. `ZAXIS` is the
default propagation direction used throughout the package (e.g. as the
default `dir` for a new geometry).
"""
const ZAXIS = Vec3(0., 0., 1.)
const YAXIS = Vec3(0.,1.,0.)
const XAXIS = Vec3(1.,0.,0.)

"""
    ∞

Alias for `Inf`, used as an "infinite" aperture size or radius of
curvature (e.g. a flat surface, or an unbounded aperture).
"""
const ∞ = Inf

"""
    EPSILON_ANGLE

Small angular tolerance (radians). Used to avoid a duplicate endpoint
when building a closed angular range, e.g.
`0.0:deltaangle:(2pi - EPSILON_ANGLE)` in `src/characterization.jl`.
"""
const EPSILON_ANGLE = 1e-7

"""
    refIndexDefault

Mutable global holding the default/ambient refractive index (e.g. air,
1.0) used throughout the package as the default index on either side of
a surface. Read via [`rInDef`](@ref); change it via [`setRInDef`](@ref)
rather than assigning to it directly.
"""
global refIndexDefault::Float64 = 1.0

"""
    rInDef()

Return the current default refractive index, `refIndexDefault`.
"""
rInDef()=refIndexDefault

"""
    rInDef(wave)

Same as `rInDef()`, ignoring `wave`. Exists so `rInDef` has the same
`wavelength -> index` signature as other glass functions, letting it be
used directly as a [`defaultGlassCatalog`](@ref) entry.
"""
rInDef(wave) = refIndexDefault #used as a glass function in the glass dictionary

"""
    setRInDef(in)

Set the default refractive index, `refIndexDefault`, to `in`.
"""
function setRInDef(in)
    global refIndexDefault = in
end

"""
    defaultGlassCatalog

Global `Dict{AbstractString, Any}` mapping material-name strings to
`wavelength -> index` functions. Seeded here with a `"DEFAULT"` entry
([`rInDef`](@ref)); populated further by `loadRICatalog!`
(`src/lens_refractive_index.jl`). Consumed by, e.g.,
`zemaxsurfsToGeo`/`viewZemaxFile` (`src/zemax.jl`) as the `glassCatalog`
default.
"""
global defaultGlassCatalog = Dict{AbstractString, Any}()
defaultGlassCatalog["DEFAULT"]= rInDef #default refractive index function

"""
    identityPol

The 3x3 identity matrix, as an `SMatrix{3,3,Float64}`. Used as the
no-op `p`/`o` polarization matrices in `AmpData` (`src/lens_definitions.jl`)
for a surface that doesn't alter polarization -- see `AMPPERFECTSURFACE`
and `identityAmpMats()` in `src/tracing.jl`.
"""
const identityPol = SMatrix{3, 3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
