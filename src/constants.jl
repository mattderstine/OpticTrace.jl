
export ∞, ORIGIN, ZAXIS, YAXIS, XAXIS, refIndexDefault, rInDef, setRInDef
export LENGTH_UNIT, WAVELENGTH_UNIT, LENGTH_TO_WAVELENGTH

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
    LENGTH_UNIT

The canonical physical-length unit this package's geometry is always
expressed in: `"mm"`. Ray positions, surface curvature/radius/thickness,
aperture sizes, aspheric/polynomial coefficients, and the built-in
Edmund/Thorlabs catalog lenses (`src/lens_edmund.jl`/
`src/lens_thorlabs.jl`) all assume this. `src/zemax.jl`'s Zemax import
pipeline (`readZemax`, via `convertZemaxUnitsToMM!`) converts every
length-dimensioned field from the source `.zmx` file's own `UNIT` to
this unit at parse time, so nothing downstream (tracing, OPD,
characterization, plotting) needs to consult a per-system unit tag --
this constant documents that invariant; nothing in this package branches
on its value today.
"""
const LENGTH_UNIT = "mm"

"""
    WAVELENGTH_UNIT

The canonical wavelength unit: `"μm"` (micrometers) -- Zemax's own
native `WAVM` wavelength unit (independent of a `.zmx` file's `UNIT`
keyword, which governs physical lens dimensions only, never
wavelengths), and this package's convention throughout
(`ExtendedGeometry.wavelength`, `OpticalSystem.wavelengths`, the glass
refractive-index formulas in `src/lens_refractive_index.jl`). Like
[`LENGTH_UNIT`](@ref), this documents an invariant rather than a
switchable setting.
"""
const WAVELENGTH_UNIT = "μm"

"""
    LENGTH_TO_WAVELENGTH

Conversion factor between [`LENGTH_UNIT`](@ref) (mm) and
[`WAVELENGTH_UNIT`](@ref) (μm): `0.001`, i.e. 1 `WAVELENGTH_UNIT` = `0.001`
`LENGTH_UNIT` (1 μm = 0.001 mm). The single source of truth for any
calculation that combines a physical length with a wavelength -- e.g.
converting an optical path difference into a number of waves -- used by
`opdRel`'s callers in `src/plotting.jl` instead of a hardcoded
`1000.0`/`1e-3` literal. Both directions of use:
- multiply a `WAVELENGTH_UNIT`-valued number by this to convert it to
  `LENGTH_UNIT` (μm → mm);
- divide a `LENGTH_UNIT`-valued number by this to convert it to
  `WAVELENGTH_UNIT` (mm → μm).

Name and value are deliberately decoupled from "mm"/"μm" specifically --
if [`LENGTH_UNIT`](@ref)/[`WAVELENGTH_UNIT`](@ref) ever changed, only
this factor's value would need updating, not every call site
referencing it.
"""
const LENGTH_TO_WAVELENGTH = 0.001

"""
    identityPol

The 3x3 identity matrix, as an `SMatrix{3,3,Float64}`. Used as the
no-op `p`/`o` polarization matrices in `AmpData` (`src/lens_definitions.jl`)
for a surface that doesn't alter polarization -- see `AMPPERFECTSURFACE`
and `identityAmpMats()` in `src/tracing.jl`.
"""
const identityPol = SMatrix{3, 3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
