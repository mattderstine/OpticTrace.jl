#=
Wavefront error analysis -- see plans/wavefront_analysis.md.

Step 1 (this file): field/aimer/apodization types, per-ray optical path
length with a per-surface phase hook, chief ray, real-ray exit pupil,
analytic reference sphere, and `wavefrontMap`, which returns a
`WavefrontResult` (OPD on a Cartesian pupil grid, in waves).

Statistics (RMS, PV, best focus, Strehl), the PSF, and plots are later
steps and take a `WavefrontResult`.
=#

export AbstractFieldPoint, ObjectPoint, CollimatedField
export AbstractRayAimer, StopCenterAim, aimRay
export AbstractPupilAmplitude, UniformPupil, GaussianPupil, FunctionPupil, AnglePupil, pupilAmplitude
export surfacePhase, bendPhase, pathLength
export WavefrontResult, wavefrontMap
export wavefrontRMS, wavefrontPV, wavefrontStats, strehlMarechal, strehlRatio
export refocusOPD, bestFocus, rmsWavefrontLegacy
export PSFResult, wavefrontPSF, spotMetrics, geometricSpot
export wavefrontRemove

"""
    WAVEFRONT_OPD_SIGN

Sign applied to `OPL_ray - OPL_chief` to give the reported OPD. `+1`
means a ray with a longer optical path than the chief ray has a positive
OPD (the convention of `opdRel`). Whether this matches Zemax's sign
convention has not been verified (`plans/wavefront_analysis.md`, Step 0
item 5) -- flip this one constant if it turns out to be opposite.
"""
const WAVEFRONT_OPD_SIGN = 1.0

#=
Object / field points
=#

"""
    AbstractFieldPoint{T<:Real}

Abstract supertype for the source of the rays in a wavefront analysis:
[`ObjectPoint`](@ref) (a finite point source) or [`CollimatedField`](@ref)
(an object at infinity, a collimated beam from one direction).
"""
abstract type AbstractFieldPoint{T<:Real} end

"""
    ObjectPoint(point::Point3)

A finite point source at `point`, in global coordinates. All rays start
in phase at `point` (a spherical source wavefront), so the path from the
point to the first surface is part of every ray's optical path.
"""
struct ObjectPoint{T<:Real} <: AbstractFieldPoint{T}
    point::Point3{T}
end

"""
    CollimatedField(dir::Vec3)

An object at infinity: a collimated beam travelling along `dir` (global
coordinates; normalized on construction). All rays start in phase on the
plane through the center of the pupil surface normal to `dir` (a plane
source wavefront). Shifting that plane along `dir` shifts every ray's
optical path equally, so it does not change the OPD.
"""
struct CollimatedField{T<:Real} <: AbstractFieldPoint{T}
    dir::Vec3{T}
    CollimatedField(dir::Vec3{T}) where {T<:Real} = new{T}(Vec3{T}(normalize(dir)))
end

#=
Ray aiming
=#

"""
    AbstractRayAimer

Abstract supertype for the strategy that turns a normalized pupil
coordinate `(u, v)` into a launch ray. This is the hook for chief-ray /
real-ray aiming: a future `IterativeAim <: AbstractRayAimer` would solve
for the launch direction that lands the ray on the stop at `(u, v)`.
The same aimer launches the chief ray and every grid ray.

Interface: `aimRay(aimer, geo, field, u, v) -> (status, Ray)`, plus
`pupilPoint(aimer, geo, u, v)` and `pupilSurfaceIndex(aimer, geo)`.
"""
abstract type AbstractRayAimer end

"""
    StopCenterAim(pupil = 1)

Paraxial-style aiming (no iteration): the launch ray for pupil coordinate
`(u, v)` is aimed straight at the point `(u*R, v*R)` in the local `z = 0`
plane of the pupil surface, where `R` is that surface's `sizeOptic`. The
pupil surface is a flat sampling plane, so surface sag is ignored when
aiming, and refraction by surfaces in front of it is not accounted for
(the pupil surface should be the first surface in `geo`; this is
required by `wavefrontMap` in this version).

`pupil` is the surface's index in `geo` or its `surfname`.
"""
struct StopCenterAim <: AbstractRayAimer
    pupil::Union{Int,String}
end
StopCenterAim() = StopCenterAim(1)

"""
    pupilSurfaceIndex(aimer::StopCenterAim, geo)

Index in `geo` of `aimer`'s pupil surface. Throws an `ArgumentError` if a
named surface is not found.
"""
function pupilSurfaceIndex(aimer::StopCenterAim, geo)
    if aimer.pupil isa Int
        return aimer.pupil
    end
    idx = findfirst(s -> s.surfname == aimer.pupil, geo)
    isnothing(idx) && throw(ArgumentError("pupil surface \"$(aimer.pupil)\" not found in geo"))
    idx
end

"""
    pupilPoint(aimer::StopCenterAim, geo, u, v)

Global coordinates of the point at normalized pupil coordinate `(u, v)`
(`(0, 0)` is the pupil center, the unit circle is the pupil edge) in the
local `z = 0` plane of `aimer`'s pupil surface.
"""
function pupilPoint(aimer::StopCenterAim, geo, u, v)
    s = geo[pupilSurfaceIndex(aimer, geo)]
    R = sizeOptic(s)
    s.toGlobalCoord(Point3(u * R, v * R, zero(R)))
end

"""
    aimRay(aimer::StopCenterAim, geo, field::ObjectPoint, u, v)

Launch ray from the point source `field.point` toward the pupil point
`(u, v)`. Returns `(0, ray)`; the status slot is reserved for aimers that
can fail (an iterative aimer failing to converge).
"""
function aimRay(aimer::StopCenterAim, geo, field::ObjectPoint{T}, u, v) where T<:Real
    P = pupilPoint(aimer, geo, u, v)
    dir = Vec3(normalize(P .- field.point))
    return 0, Ray(field.point, Vec3{T}(dir))
end

"""
    aimRay(aimer::StopCenterAim, geo, field::CollimatedField, u, v)

Launch ray for a collimated beam: parallel to `field.dir`, starting on the
plane through the pupil center normal to `field.dir`, at the point from
which it reaches the pupil point `(u, v)`. Returns `(0, ray)`.
"""
function aimRay(aimer::StopCenterAim, geo, field::CollimatedField{T}, u, v) where T<:Real
    P0 = pupilPoint(aimer, geo, 0, 0)
    P = pupilPoint(aimer, geo, u, v)
    d = field.dir
    t = sum((P .- P0) .* d)
    base = Point3{T}(P .- t .* d)
    return 0, Ray(base, d)
end

#=
Pupil amplitude (apodization / input beam profile)
=#

"""
    AbstractPupilAmplitude

Abstract supertype for the amplitude (not intensity) weighting of rays
across the pupil -- an apodized pupil or an input beam with a non-uniform
profile. Interface: `pupilAmplitude(a, u, v, ray) -> amplitude`, where
`(u, v)` is the normalized pupil coordinate and `ray` is the launch ray.
"""
abstract type AbstractPupilAmplitude end

"""
    UniformPupil()

Uniform illumination: amplitude 1 everywhere.
"""
struct UniformPupil <: AbstractPupilAmplitude end

"""
    GaussianPupil(w)

Gaussian amplitude `exp(-(u^2 + v^2)/w^2)` in normalized pupil
coordinates; `w` is the amplitude `1/e` radius (the intensity `1/e^2`
radius) as a fraction of the pupil radius.
"""
struct GaussianPupil{T<:Real} <: AbstractPupilAmplitude
    w::T
end

"""
    FunctionPupil(f)

Amplitude `f(u, v)` as a function of the normalized pupil coordinate.
"""
struct FunctionPupil{F} <: AbstractPupilAmplitude
    f::F
end

"""
    AnglePupil(f)

Amplitude `f(dir)` as a function of the launch ray's direction (a
`Vec3`), for sources with an angular emission profile.
"""
struct AnglePupil{F} <: AbstractPupilAmplitude
    f::F
end

"""
    pupilAmplitude(a::UniformPupil, u, v, ray)

Returns `1.0`. See [`AbstractPupilAmplitude`](@ref) for the interface.
"""
pupilAmplitude(a::UniformPupil, u, v, ray) = 1.0

"""
    pupilAmplitude(a::GaussianPupil, u, v, ray)

Returns `exp(-(u^2 + v^2)/a.w^2)`. See [`AbstractPupilAmplitude`](@ref).
"""
pupilAmplitude(a::GaussianPupil, u, v, ray) = exp(-(u^2 + v^2) / a.w^2)

"""
    pupilAmplitude(a::FunctionPupil, u, v, ray)

Returns `a.f(u, v)`. See [`AbstractPupilAmplitude`](@ref).
"""
pupilAmplitude(a::FunctionPupil, u, v, ray) = a.f(u, v)

"""
    pupilAmplitude(a::AnglePupil, u, v, ray)

Returns `a.f(ray.dir)`. See [`AbstractPupilAmplitude`](@ref).
"""
pupilAmplitude(a::AnglePupil, u, v, ray) = a.f(ray.dir)

#=
Optical path length with a phase hook
=#

"""
    bendPhase(mod::AbstractBendType, s, rayIn, trcOut, λ)

Phase added by a surface's `mod` ("bend") type, as an optical path length
in `LENGTH_UNIT` (mm), on top of the geometric `nIn * delta` term. The
default is `0`. Add a method for a new `AbstractBendType` to model a
grating (`m*λ*x/d`), a diffractive optic, or an ideal lens
(`TODO.md` Code issues #41). `s` is the surface, `rayIn` the ray arriving
at it, `trcOut` the [`Trace`](@ref) element it produced, and `λ` the
evaluation wavelength in `WAVELENGTH_UNIT` (μm) -- convert with
`LENGTH_TO_WAVELENGTH` when the result is a length.
"""
bendPhase(mod::AbstractBendType, s, rayIn, trcOut, λ) = zero(trcOut.delta)

"""
    surfacePhase(s::AbstractSurface, rayIn, trcOut, λ)

Extra optical path length (mm) a surface adds beyond `nIn * delta`;
`0` for `ModelSurface`s. See [`bendPhase`](@ref).
"""
surfacePhase(s::AbstractSurface, rayIn, trcOut, λ) = zero(trcOut.delta)

"""
    surfacePhase(s::OptSurface, rayIn, trcOut, λ)

Forwards to [`bendPhase`](@ref) for `s.mod`.
"""
surfacePhase(s::OptSurface, rayIn, trcOut, λ) = bendPhase(s.mod, s, rayIn, trcOut, λ)

"""
    pathLength(trc, geo, λ; upto = length(trc))

Optical path length in mm of a trace: `Σ nIn*delta` over `trc[2:upto]`
(`delta` is the distance travelled to reach that surface and `nIn` the
index of the medium it travelled in), plus each surface's
[`surfacePhase`](@ref). `trc[1]` is the start ray and has no path.
`geo[i-1]` is the surface that produced `trc[i]`. `λ` is in μm and only
reaches the phase hooks.
"""
function pathLength(trc, geo, λ; upto = length(trc))
    opl = zero(trc[1].delta)
    for i in 2:upto
        opl += trc[i].nIn * trc[i].delta + surfacePhase(geo[i-1], trc[i-1].ray, trc[i], λ)
    end
    opl
end

#=
Reference sphere
=#

"""
    oplToReference(b, d, opl0, n, C, R, σ)

Optical path length (mm) from the start of the trace to the reference
sphere, up to a constant that is the same for every ray (so that
differences between rays are exact). The ray's last real segment is
`b + t*d` (global; `b` is the point on the last real surface, `d` the unit
direction), `opl0` the path length to `b`, `n` the image-space index. The
reference sphere is centered at `C` with radius `R` (`Inf` for a plane
reference, i.e. an exit pupil at infinity); `σ = +1` if the sphere's
vertex (exit pupil) is before `C` along the propagation, `-1` if after.

Uses `t = τ - σ*sqrt(R^2 - h^2)` with `τ = d·(C - b)` and `h` the ray's
distance of closest approach to `C`, written as
`τ - σR + σq`, `q = h²/(R + sqrt(R² - h²))`, so the large constant `σR`
is dropped without cancellation error. Returns `NaN` if the ray misses
the sphere.
"""
function oplToReference(b, d, opl0, n, C, R, σ)
    w = Vec3(C .- b)
    τ = sum(d .* w)
    h2 = sum(abs2, w .- τ .* d)
    if isinf(R)
        q = zero(τ)
    else
        s2 = R^2 - h2
        s2 < 0 && return oftype(τ, NaN)
        q = h2 / (R + sqrt(s2))
    end
    opl0 + n * (τ + σ * q)
end

"""
    referencePoint(b, d, C, R, σ, bc, dc)

Global point where the ray `b + t*d` crosses the reference sphere
(centered `C`, radius `R`, vertex side `σ`; see [`oplToReference`](@ref)),
or, for `R == Inf`, the plane through the chief ray's last-surface point
`bc` normal to the chief direction `dc`. Used to record the exit-pupil
coordinates of each ray.
"""
function referencePoint(b, d, C, R, σ, bc, dc)
    if isinf(R)
        t = sum(dc .* (bc .- b)) / sum(dc .* d)
    else
        w = Vec3(C .- b)
        τ = sum(d .* w)
        h2 = sum(abs2, w .- τ .* d)
        s2 = R^2 - h2
        s2 < 0 && return oftype(b, Point3(NaN, NaN, NaN))
        t = τ - σ * sqrt(s2)
    end
    Point3(b .+ t .* d)
end

"""
    findExitPupil(geoPre, aimer, launchDir, bc, dc, probeAngle)

Locate the exit pupil on the chief ray, in 3D, in global coordinates.
Traces four probe rays from the pupil center, tilted by `probeAngle` rad
about the chief ray's launch direction `launchDir` in two perpendicular planes, and finds the
point on the chief ray's image-space line that is nearest (least
squares) to the four image-space probe lines -- the image of the stop
center. (Probe rays launched from the object point toward other pupil
points would converge at the image point, not the pupil.)

`bc`, `dc` are the chief ray's base and direction after the last real
surface (image space); `launchDir` is its direction before the system. Returns `(status, s, spread)`: `status` is `0` or `1` if a probe
ray failed to trace; `s` the signed distance of the exit pupil from `bc`
along `dc`; `spread` the RMS sine of the angle between the probe rays
and the chief ray in image space (small when the pupil is at infinity).
"""
function findExitPupil(geoPre, aimer, launchDir, bc, dc, probeAngle)
    P0 = pupilPoint(aimer, geoPre, 0, 0)
    ref = abs(launchDir[2]) < 0.9 ? Vec3(0., 1., 0.) : Vec3(1., 0., 0.)
    e1 = Vec3(normalize(cross(ref, launchDir)))
    e2 = Vec3(cross(launchDir, e1))
    num = 0.0
    den = 0.0
    for dirk in (launchDir + probeAngle * e1, launchDir - probeAngle * e1, launchDir + probeAngle * e2, launchDir - probeAngle * e2)
        status, trc = traceGeometry(Ray(P0, Vec3(normalize(dirk))), geoPre)
        status != 0 && return 1, NaN, NaN
        rk = trc[end].ray
        pdc = dc .- sum(dc .* rk.dir) .* rk.dir
        num += sum(pdc .* (rk.base .- bc))
        den += sum(abs2, pdc)
    end
    return 0, num / den, sqrt(den / 4)
end

#=
Result
=#

"""
    WavefrontResult{T}

Result of [`wavefrontMap`](@ref): the wavefront error on a square
`N x N` Cartesian grid of normalized pupil coordinates, with everything
later analysis steps need (so nothing is re-traced).

Grid and values (`[i, j]` is the point `(u[i], v[j])`; the pupil is the
unit disk):
- `u`, `v` -- grid axes, `range(-1, 1, N)` (`N` odd, so a center row/column exists)
- `opd` -- OPD in waves of `λ`, relative to the chief ray (zero at the
  pupil center); `NaN` where `mask` is `false`
- `mask` -- `true` where the point is in the pupil and the ray traced
  through to the reference sphere
- `weight` -- pupil amplitude (apodization), `0` where masked
- `pupilPoints` -- global point where each ray crosses the reference
  sphere (plane for an exit pupil at infinity): the exit-pupil coordinates

Per-ray image-space data (global coordinates; `NaN` where masked), enough
to re-evaluate the reference without re-tracing (best focus):
- `lastBase`, `lastDir` -- point on the last real surface and direction after it
- `oplToLast` -- optical path length (mm) to `lastBase`

Reference:
- `chiefBase`, `chiefDir`, `chiefOPL` -- chief ray's last real segment and path length to it
- `imagePoint` -- where the chief ray meets the image surface (reference sphere center)
- `exitPupil` -- exit pupil location on the chief ray (`NaN`s if at infinity)
- `refRadius` -- reference sphere radius, `Inf` for a plane reference
- `refSide` -- `+1` if the exit pupil is before the image point along the chief ray, else `-1`
- `nImage` -- image-space refractive index (as seen by the chief ray)

Provenance: `λ` (μm), `geoWavelength` (μm, or `nothing` if unknown), `field`.
"""
struct WavefrontResult{T<:Real}
    u::Vector{T}
    v::Vector{T}
    opd::Matrix{T}
    mask::Matrix{Bool}
    weight::Matrix{T}
    pupilPoints::Matrix{Point3{T}}
    lastBase::Matrix{Point3{T}}
    lastDir::Matrix{Vec3{T}}
    oplToLast::Matrix{T}
    chiefBase::Point3{T}
    chiefDir::Vec3{T}
    chiefOPL::T
    imagePoint::Point3{T}
    exitPupil::Point3{T}
    refRadius::T
    refSide::Int
    nImage::T
    λ::T
    geoWavelength::Union{Nothing,T}
    field::AbstractFieldPoint{T}
end

#=
Entry point
=#

"""
    wavefrontMap(geo, field, λ; N = 65, aimer = StopCenterAim(), apodization = UniformPupil(),
                 imageSurface = "end", geoWavelength = nothing,
                 probeAngle = 1e-3, pupilInfinityRatio = 1e-3)

Wavefront error of the system `geo` for the source `field`, at the
evaluation wavelength `λ` (μm), relative to the chief ray from `field`.

For each point of an `N x N` grid over the unit pupil disk
(`aimer` supplies the launch rays), the ray is traced through the
surfaces before `imageSurface`, and its optical path (`pathLength`,
including any [`surfacePhase`](@ref)) to a reference sphere is compared
with the chief ray's. The reference sphere is centered on the point where
the chief ray meets `imageSurface`, passes through the real-ray exit
pupil ([`findExitPupil`](@ref)), and lies in the actual image-space
refractive index; for an exit pupil at infinity it becomes a plane
(`refRadius == Inf`). Returns a [`WavefrontResult`](@ref), or `nothing`
(after a warning) if the analysis cannot be done.

Wavelength: indices in `geo` are fixed numbers set when `geo` was built,
so `λ` is only used to convert OPD to waves and in the phase hooks.
`geoWavelength` is the wavelength `geo` was built for; a warning is issued
if it differs from `λ` or is not given.

Not supported (warns and returns `nothing`): a `CDiffuser` anywhere in
`geo` (random scatter, so OPD is meaningless). Warns but continues: a
`ParaxialLensT` (adds no phase, `TODO.md` Code issues #41, so OPD
through it is wrong). `OptSurface`s do not clip rays at their `aperture`
(only `ModelSurface`s do), so the mask only reflects misses, TIR, and
`ModelSurface` clipping.

`N` must be odd (so the pupil center is a grid point). `probeAngle` (rad)
is the exit-pupil probe tilt; the pupil counts as at infinity when the RMS
sine of the probe rays' angle to the chief ray in image space is below
`pupilInfinityRatio * probeAngle`. The pupil surface (`aimer`) must be
`geo[1]` in this version.
"""
function wavefrontMap(geo, field::AbstractFieldPoint{T}, λ::Real;
        N::Integer = 65,
        aimer::AbstractRayAimer = StopCenterAim(),
        apodization::AbstractPupilAmplitude = UniformPupil(),
        imageSurface = "end",
        geoWavelength = nothing,
        probeAngle::Real = 1e-3,
        pupilInfinityRatio::Real = 1e-3) where T<:Real

    (N >= 3 && isodd(N)) || throw(ArgumentError("N must be odd and >= 3, got $N"))
    pupilSurfaceIndex(aimer, geo) == 1 ||
        throw(ArgumentError("the pupil surface must be geo[1] in this version"))

    for s in geo
        if s isa OptSurface && s.mod isa CDiffuser
            @warn "wavefrontMap: surface \"$(s.surfname)\" is a CDiffuser (random scatter), which is incompatible with OPD calculations. Stopping."
            return nothing
        end
    end
    if any(s -> s isa OptSurface && s.mod isa ParaxialLensT, geo)
        @warn "wavefrontMap: geo contains a ParaxialLensT, which adds no phase (TODO.md Code issues #41); OPD through it is wrong." maxlog=1
    end
    if isnothing(geoWavelength)
        @warn "wavefrontMap: geoWavelength not given; assuming geo's refractive indices match λ = $λ μm." maxlog=1
    elseif !isapprox(geoWavelength, λ; rtol = 1e-9)
        @warn "wavefrontMap: geo was built for $geoWavelength μm but λ = $λ μm. Indices are not changed, so chromatic effects are wrong; only phase hooks and the waves conversion use λ."
    end

    imgIdx = numsurfFromName(imageSurface, geo)
    imgIdx >= 2 || throw(ArgumentError("imageSurface must come after at least one other surface"))
    geoPre = geo[1:imgIdx-1]
    geoImg = geo[1:imgIdx]

    # chief ray, through to the image surface
    status, chiefRay = aimRay(aimer, geo, field, zero(T), zero(T))
    if status != 0
        @warn "wavefrontMap: could not aim the chief ray"
        return nothing
    end
    status, trcC = traceGeometry(chiefRay, geoImg)
    if status != 0
        @warn "wavefrontMap: chief ray did not reach the image surface (status $status)"
        return nothing
    end
    imagePoint = trcC[end].ray.base
    nImage = trcC[end].nIn
    lastC = trcC[end-1].ray
    bc, dc = lastC.base, lastC.dir
    oplC = pathLength(trcC, geoPre, λ; upto = length(trcC) - 1)

    # exit pupil
    status, sPupil, spread = findExitPupil(geoPre, aimer, chiefRay.dir, bc, dc, probeAngle)
    if status != 0
        @warn "wavefrontMap: exit pupil probe rays did not trace through the system"
        return nothing
    end
    sImage = sum(dc .* (imagePoint .- bc))
    if spread < pupilInfinityRatio * probeAngle
        R = T(Inf)
        σ = 1
        exitPupil = Point3{T}(NaN, NaN, NaN)
    else
        R = abs(sImage - sPupil)
        σ = sPupil < sImage ? 1 : -1
        exitPupil = Point3{T}(bc .+ sPupil .* dc)
        if R < 1e-9
            @warn "wavefrontMap: the exit pupil coincides with the image point; no reference sphere"
            return nothing
        end
    end

    oplRefC = oplToReference(bc, dc, oplC, nImage, imagePoint, R, σ)
    nanP = Point3{T}(NaN, NaN, NaN)
    nanV = Vec3{T}(NaN, NaN, NaN)
    us = collect(range(-one(T), one(T); length = N))
    vs = collect(range(-one(T), one(T); length = N))
    opd = fill(T(NaN), N, N)
    mask = fill(false, N, N)
    weight = zeros(T, N, N)
    pupilPoints = fill(nanP, N, N)
    lastBase = fill(nanP, N, N)
    lastDir = fill(nanV, N, N)
    oplToLast = fill(T(NaN), N, N)
    wavesPerLength = 1 / (λ * LENGTH_TO_WAVELENGTH)

    Threads.@threads for j in 1:N
        for i in 1:N
            u, v = us[i], vs[j]
            u^2 + v^2 <= 1 + 1e-12 || continue
            st, ray = aimRay(aimer, geo, field, u, v)
            st == 0 || continue
            st, trc = traceGeometry(ray, geoPre)
            st == 0 || continue
            r = trc[end].ray
            opl0 = pathLength(trc, geoPre, λ)
            oplRef = oplToReference(r.base, r.dir, opl0, nImage, imagePoint, R, σ)
            isnan(oplRef) && continue
            opd[i, j] = WAVEFRONT_OPD_SIGN * (oplRef - oplRefC) * wavesPerLength
            mask[i, j] = true
            weight[i, j] = pupilAmplitude(apodization, u, v, ray)
            lastBase[i, j] = r.base
            lastDir[i, j] = r.dir
            oplToLast[i, j] = opl0
            pupilPoints[i, j] = referencePoint(r.base, r.dir, imagePoint, R, σ, bc, dc)
        end
    end

    WavefrontResult{T}(us, vs, opd, mask, weight, pupilPoints, lastBase, lastDir, oplToLast,
        bc, dc, oplC, imagePoint, exitPupil, R, σ, nImage, T(λ),
        isnothing(geoWavelength) ? nothing : T(geoWavelength), field)
end


#=
Statistics, refocus, best focus, Strehl
=#

"""
    sampleWeights(r::WavefrontResult, opd, weighted)

Weights of the valid samples as a matrix: `r.weight` where the sample is
in the mask and `opd` is finite, `1` there if `weighted == false`, `0`
elsewhere.
"""
function sampleWeights(r::WavefrontResult{T}, opd, weighted::Bool) where T
    w = zeros(T, size(opd))
    for i in eachindex(opd)
        if r.mask[i] && isfinite(opd[i])
            w[i] = weighted ? r.weight[i] : one(T)
        end
    end
    w
end

"""
    weightedVariance(opd, w)

Weighted mean and variance `(m, v)` of `opd` over the samples with
nonzero weight `w`.
"""
function weightedVariance(opd, w)
    sw = sum(w)
    sw > 0 || return NaN, NaN
    m = sum(w[i] * opd[i] for i in eachindex(opd) if w[i] > 0) / sw
    v = sum(w[i] * (opd[i] - m)^2 for i in eachindex(opd) if w[i] > 0) / sw
    m, v
end

"""
    wavefrontRMS(r::WavefrontResult; weighted = true, opd = r.opd)

RMS wavefront error in waves with the mean (piston) removed. With
`weighted = true` each sample counts by its amplitude `r.weight`
(apodization); with `false` all valid samples count equally. `opd`
substitutes another OPD matrix on the same grid, e.g. from
[`refocusOPD`](@ref). Returns `NaN` if there are no valid samples.
"""
function wavefrontRMS(r::WavefrontResult; weighted::Bool = true, opd = r.opd)
    sqrt(weightedVariance(opd, sampleWeights(r, opd, weighted))[2])
end

"""
    wavefrontPV(r::WavefrontResult; minWeight = 0.0, opd = r.opd)

Peak-to-valley wavefront error in waves: `max - min` of the valid samples
whose weight is at least `minWeight` times the largest weight (so a
heavily apodized edge can be left out). Returns `NaN` if there are none.
PV depends on the sampling density; check it converges with `N`.
"""
function wavefrontPV(r::WavefrontResult; minWeight::Real = 0.0, opd = r.opd)
    wmax = maximum(r.weight[r.mask]; init = 0.0)
    vals = [opd[i] for i in eachindex(opd) if r.mask[i] && isfinite(opd[i]) && r.weight[i] >= minWeight * wmax]
    isempty(vals) ? NaN : maximum(vals) - minimum(vals)
end

"""
    strehlMarechal(rms)

Maréchal approximation to the Strehl ratio, `exp(-(2π rms)²)`, for an RMS
wavefront error `rms` in waves. Accurate only for `rms` below about 0.1.
"""
strehlMarechal(rms::Real) = exp(-(2π * rms)^2)

"""
    strehlRatio(r::WavefrontResult; weighted = true, opd = r.opd)

Strehl ratio at the reference point (the chief ray's image point, or the
shifted point if `opd` came from [`refocusOPD`](@ref)):
`|Σ w exp(i2πW)|² / (Σ w)²`. This is the on-axis intensity of the
scalar diffraction pattern relative to the aberration-free one. It is
the PSF peak only if the peak sits at the reference point; a wavefront
tilt (peak displaced) lowers it. The FFT-based peak Strehl comes with the
PSF computation.
"""
function strehlRatio(r::WavefrontResult; weighted::Bool = true, opd = r.opd)
    w = sampleWeights(r, opd, weighted)
    sw = sum(w)
    sw > 0 || return NaN
    abs2(sum(w[i] * cis(2π * opd[i]) for i in eachindex(opd) if w[i] > 0)) / sw^2
end

"""
    wavefrontStats(r::WavefrontResult; weighted = true, minWeight = 0.0)

NamedTuple of `rms`, `pv` (waves), `piston` (the removed weighted mean
OPD), `strehlMarechal`, `strehl` ([`strehlRatio`](@ref)), and `npoints`
(valid samples). See [`wavefrontRMS`](@ref), [`wavefrontPV`](@ref).
"""
function wavefrontStats(r::WavefrontResult; weighted::Bool = true, minWeight::Real = 0.0)
    w = sampleWeights(r, r.opd, weighted)
    m, v = weightedVariance(r.opd, w)
    rms = sqrt(v)
    (rms = rms, pv = wavefrontPV(r; minWeight), piston = m, strehlMarechal = strehlMarechal(rms),
     strehl = strehlRatio(r; weighted), npoints = count(>(0), w))
end

"""
    refocusOPD(r::WavefrontResult, δ)

OPD matrix (waves, chief-ray referenced) for the reference sphere
re-centered `δ` mm along the chief ray's image-space direction from the
image point (`δ > 0` is farther along the propagation direction). The
exit pupil is held fixed, so the sphere radius changes by `δ`. No
rays are re-traced: it re-evaluates the stored last segments.
`refocusOPD(r, 0)` equals `r.opd`. `NaN` where the shifted sphere is
missed.
"""
function refocusOPD(r::WavefrontResult{T}, δ::Real) where T
    dc, bc = r.chiefDir, r.chiefBase
    C = Point3{T}(r.imagePoint .+ δ .* dc)
    if isinf(r.refRadius)
        R, σ = T(Inf), 1
    else
        sImage = sum(dc .* (r.imagePoint .- bc))
        sPupil = sImage - r.refSide * r.refRadius
        sC = sImage + δ
        R = abs(sC - sPupil)
        σ = sPupil < sC ? 1 : -1
        R < 1e-9 && return fill(T(NaN), size(r.opd))
    end
    oplRefC = oplToReference(bc, dc, r.chiefOPL, r.nImage, C, R, σ)
    wavesPerLength = 1 / (r.λ * LENGTH_TO_WAVELENGTH)
    opd = fill(T(NaN), size(r.opd))
    for i in eachindex(opd)
        r.mask[i] || continue
        oplRef = oplToReference(r.lastBase[i], r.lastDir[i], r.oplToLast[i], r.nImage, C, R, σ)
        opd[i] = WAVEFRONT_OPD_SIGN * (oplRef - oplRefC) * wavesPerLength
    end
    opd
end

"""
    bestFocus(r::WavefrontResult; range = nothing, weighted = true)

Longitudinal shift `δ` (mm, see [`refocusOPD`](@ref)) of the reference
point along the chief ray that minimizes the RMS wavefront error, found
by a bounded 1D search (`Optim`, Brent) over `δ ∈ [-range, range]`.
`range` defaults to half the distance from the exit pupil to the image
point (or from the last surface to the image point if the pupil is at
infinity); give it explicitly if the best focus may lie farther away.
The minimum is found only inside the bracket and is a local one, so a
`δ` at the edge of `range` means the search should be widened.

Returns `(δ, rms, opd)`: the shift, the RMS (waves) there, and the OPD
matrix there. Tilt (lateral shift of the reference point) is not
optimized.
"""
function bestFocus(r::WavefrontResult; range = nothing, weighted::Bool = true)
    half = if range !== nothing
        range
    else
        sImage = sum(r.chiefDir .* (r.imagePoint .- r.chiefBase))
        isinf(r.refRadius) ? abs(sImage) / 2 : r.refRadius / 2
    end
    f(δ) = (x = wavefrontRMS(r; weighted, opd = refocusOPD(r, δ)); isfinite(x) ? x : Inf)
    res = Optim.optimize(f, -half, half)
    δ = Optim.minimizer(res)
    opd = refocusOPD(r, δ)
    δ, wavefrontRMS(r; weighted, opd), opd
end

#=
Legacy cross-check
=#

"""
    rmsWavefrontLegacy(geo, r::Point3, θmax, λ; surfview = "end", rings = 12)

RMS wavefront error (waves, piston removed) of the finite on-axis
point `r` (in the local frame of `geo[1]`) computed with the legacy
[`opdRel`](@ref) machinery of `plotOPD!`: a `refractSphere` placed at the
estimated exit pupil ([`computeExitPupilLoc`](@ref); `geo` needs a
surface named `"stop"`) and rays on `rings` equal-area rings (ring `k` at
radius `(k - 1/2)/rings` of the cone, with `4(2k - 1)` rays, `4rings²` in
total) out to the cone half-angle `θmax` (rad), equal weight per ray. Intended as an
independent cross-check for [`wavefrontMap`](@ref), which it should
match only inside the legacy assumptions: on-axis, finite object, air
(the sphere's index is `refIndexDefault`), exit pupil on the axis.
`λ` is in μm. Returns `(rms, opds)`; `opds` has `NaN` for failed rays.
"""
function rmsWavefrontLegacy(geo, r::Point3, θmax::Real, λ::Real; surfview = "end", rings::Integer = 12)
    surfnum = tracenumFromName(surfview, geo)
    localgeo = geo[1:surfnum-1]
    status, refTrace = traceGeometryRel(Ray(r, ZAXIS), localgeo)
    status == 0 || error("reference ray did not reach $surfview")
    status, zEP = computeExitPupilLoc(geo; epsilon = 0.0001, format = "quiet")
    status == 0 || error("could not compute the exit pupil location")
    if abs(zEP) > 1e5
        zEP = refTrace[end-1].ray.base[3]
    end
    baseImage = refTrace[end].ray.base
    dirImage = refTrace[end].ray.dir
    radius = (baseImage[3] - zEP) / dirImage[3]
    testgeo = [localgeo[1:end-1]
        [refractSphere("reference sphere", Point3(0., 0., zEP), dirImage, refIndexDefault, refIndexDefault, 1.0 / radius, 25., "testcoat")]]
    status, refTrace = traceGeometryRel(Ray(r, ZAXIS), testgeo)
    status == 0 || error("reference ray did not reach the reference sphere")
    sinMax = sin(θmax)
    # equal-area sampling: ring k is centered at radius (k - 1/2)/rings and has 4(2k - 1) rays
    angles = Tuple{Float64,Float64}[]
    for k in 1:rings, j in 0:(4(2k-1)-1)
        push!(angles, (sinMax * (k - 0.5) / rings, 2π * (j + 0.5) / (4(2k - 1))))
    end
    opds = map(angles) do (s, φ)
        dir = Vec3(s * cos(φ), s * sin(φ), sqrt(1 - s^2))
        opdRel(Ray(r, dir), refTrace, testgeo) / LENGTH_TO_WAVELENGTH / λ
    end
    good = filter(isfinite, opds)
    m = sum(good) / length(good)
    sqrt(sum(abs2, good .- m) / length(good)), opds
end


"""
    wavefrontRemove(r::WavefrontResult; remove = :piston, weighted = true, opd = r.opd)

OPD matrix (waves) with low-order terms fitted and subtracted by weighted
least squares over the valid samples, for display or for statistics that
should exclude them. `remove` is `:none`, `:piston` (the weighted mean),
`:pistonTilt` (piston plus a plane in the pupil coordinates `(u, v)`) or
`:pistonTiltFocus` (those plus `u² + v²`). `NaN` where there is no valid
sample.
"""
function wavefrontRemove(r::WavefrontResult{T}; remove::Symbol = :piston, weighted::Bool = true, opd = r.opd) where T
    remove in (:none, :piston, :pistonTilt, :pistonTiltFocus) ||
        throw(ArgumentError("remove must be :none, :piston, :pistonTilt or :pistonTiltFocus, got :$remove"))
    out = T.(opd)
    remove === :none && return out
    w = sampleWeights(r, opd, weighted)
    idx = findall(>(0), vec(w))
    isempty(idx) && return fill(T(NaN), size(opd))
    ci = CartesianIndices(w)
    uq = [r.u[ci[i][1]] for i in idx]
    vq = [r.v[ci[i][2]] for i in idx]
    cols = remove === :piston ? [ones(T, length(idx))] :
           remove === :pistonTilt ? [ones(T, length(idx)), uq, vq] :
           [ones(T, length(idx)), uq, vq, uq .^ 2 .+ vq .^ 2]
    X = hcat(cols...)
    sw = sqrt.(w[idx])
    coef = (sw .* X) \ (sw .* opd[idx])
    fit = X * coef
    for (k, i) in enumerate(idx)
        out[i] = opd[i] - fit[k]
    end
    for i in eachindex(out)
        w[i] > 0 || (out[i] = T(NaN))
    end
    out
end

#=
Point spread function and spot size / shape
=#

"""
    PSFResult{T}

Diffraction (wavefront-limited) point spread function from
[`wavefrontPSF`](@ref), on a square grid of lateral positions in the plane
through the reference point (the chief ray's image point) normal to the
chief ray's image-space direction.

- `x`, `y` -- positions (mm) along `axisX`, `axisY` (orthonormal global
  vectors in that plane, with `axisX × axisY` along the chief direction)
- `intensity[i, j]` -- at `(x[i], y[j])`, normalized so that the peak of the
  aberration-free PSF of the same pupil is `1` (so `maximum(intensity)` is
  the peak Strehl ratio, to the resolution of the grid)
- `pitch` -- grid spacing (mm)
- `totalEnergy` -- the integral of `intensity` over the whole plane (mm²),
  from Parseval's theorem for the sampled pupil; used to normalize encircled
  energy so that it does not depend on the window size
- `period` -- the PSF repeats every `period` mm (sampling the pupil on a
  grid aliases it); the window should be well inside `period / 2`
- `strehl` -- `maximum(intensity)`
- `pupilDistortion` -- largest deviation of the rays' direction cosines
  from a linear map of the pupil grid, relative to the largest direction
  cosine; `pupilAnisotropy` -- relative difference of that map's two
  scales. Both are `0` for a pupil that maps uniformly
- `λ` (μm), `nImage` -- as in [`WavefrontResult`](@ref)
"""
struct PSFResult{T<:Real}
    x::Vector{T}
    y::Vector{T}
    intensity::Matrix{T}
    axisX::Vec3{T}
    axisY::Vec3{T}
    pitch::T
    totalEnergy::T
    period::T
    strehl::T
    pupilDistortion::T
    pupilAnisotropy::T
    λ::T
    nImage::T
end

"""
    defaultImageAxes(dc)

Orthonormal axes `(e1, e2)` normal to the unit direction `dc`: `e1` is the
global x axis projected normal to `dc` (global y if x is nearly along
`dc`), and `e2 = dc × e1`.
"""
function defaultImageAxes(dc)
    e1 = Vec3(1., 0., 0.) .- sum(Vec3(1., 0., 0.) .* dc) .* dc
    if sum(abs2, e1) < 1e-12
        e1 = Vec3(0., 1., 0.) .- dc[2] .* dc
    end
    e1 = Vec3(normalize(e1))
    e1, Vec3(cross(dc, e1))
end

"""
    wavefrontPSF(r::WavefrontResult; opd = r.opd, weighted = true, pixels = 129,
                 pitch = nothing, halfWidth = nothing, axes = nothing)

Wavefront-limited PSF by direct (scalar Debye) summation over the traced
rays: `U(x) = Σ w exp(i2π W) exp(i2π s·x / λ')`, where `W` is the OPD
(physical sign, `WAVEFRONT_OPD_SIGN` undone), `s` the ray's image-space
direction cosines along the grid axes, `x` the lateral position relative to
the reference point, and `λ' = λ / nImage`. Each ray carries equal
quadrature weight, so a pupil whose grid maps non-uniformly onto the exit
pupil (see `pupilDistortion`) is only approximately weighted; obliquity
factors are ignored. Because the sum separates in `x` and `y` it is two
matrix products, and the output grid is independent of the pupil grid.

Output grid: `pixels` (odd) points per side. Give either `pitch` (mm) or
`halfWidth` (mm), not both; the default pitch is `λ' / (8 NA)` with `NA` the
largest direction cosine (twice finer than the intensity Nyquist spacing
`λ' / (4 NA)`). A warning is issued if the pitch is coarser than Nyquist or
the window reaches past half the aliasing `period`. `axes = (e1, e2)`
overrides [`defaultImageAxes`](@ref). `opd` and `weighted` are as in
[`wavefrontRMS`](@ref); pass `opd = refocusOPD(r, δ)` for the PSF at a
shifted reference point.
"""
function wavefrontPSF(r::WavefrontResult{T}; opd = r.opd, weighted::Bool = true, pixels::Integer = 129,
        pitch = nothing, halfWidth = nothing, axes = nothing) where T<:Real
    (pixels >= 3 && isodd(pixels)) || throw(ArgumentError("pixels must be odd and >= 3, got $pixels"))
    (pitch === nothing || halfWidth === nothing) || throw(ArgumentError("give pitch or halfWidth, not both"))
    w = sampleWeights(r, opd, weighted)
    idx = findall(>(0), vec(w))
    isempty(idx) && throw(ArgumentError("no valid samples"))
    e1, e2 = axes === nothing ? defaultImageAxes(r.chiefDir) : axes
    ci = CartesianIndices(w)
    uq = [r.u[ci[i][1]] for i in idx]
    vq = [r.v[ci[i][2]] for i in idx]
    sx = [sum(r.lastDir[i] .* e1) for i in idx]
    sy = [sum(r.lastDir[i] .* e2) for i in idx]
    NA = maximum(hypot.(sx, sy))
    NA > 1e-9 || throw(ArgumentError("the image-space beam has no angular extent; no PSF"))

    # how uniformly the pupil grid maps onto direction cosines
    Q = hcat(uq, vq)
    S = hcat(sx, sy)
    At = Q \ S
    resid = S .- Q * At
    distortion = maximum(sqrt.(sum(abs2, resid; dims = 2))) / NA
    sv = svdvals(At)
    anisotropy = (sv[1] - sv[2]) / ((sv[1] + sv[2]) / 2)
    dAs = abs(det(At)) * (r.u[2] - r.u[1]) * (r.v[2] - r.v[1])   # direction-cosine area per ray
    distortion > 0.05 && @warn "wavefrontPSF: the pupil grid maps non-uniformly onto the exit pupil (distortion $(round(distortion; digits = 3))); the PSF is only approximate." maxlog = 1

    λp = r.λ * LENGTH_TO_WAVELENGTH / r.nImage
    period = λp / sqrt(dAs)
    dx = if pitch !== nothing
        T(pitch)
    elseif halfWidth !== nothing
        T(halfWidth) / ((pixels - 1) / 2)
    else
        λp / (8NA)
    end
    dx > λp / (4NA) * (1 + 1e-9) && @warn "wavefrontPSF: pitch $dx mm is coarser than the intensity Nyquist spacing $(λp / (4NA)) mm"
    dx * (pixels - 1) / 2 > period / 2 && @warn "wavefrontPSF: the window reaches past half the aliasing period ($period mm); increase the pupil sampling N or reduce the window"
    half = (pixels - 1) ÷ 2
    xs = collect(T, (-half:half) .* dx)
    Ex = cis.((2π / λp) .* (xs * sx'))
    Ey = cis.((2π / λp) .* (xs * sy'))
    wv = w[idx]
    P = wv .* cis.((2π * WAVEFRONT_OPD_SIGN) .* opd[idx])
    U = (Ex .* transpose(P)) * transpose(Ey)
    sw = sum(wv)
    intensity = abs2.(U) ./ sw^2
    total = sum(abs2, wv) * λp^2 / (dAs * sw^2)
    PSFResult{T}(xs, copy(xs), intensity, e1, e2, dx, total, period, maximum(intensity), distortion, anisotropy,
        r.λ, r.nImage)
end

"""
    levelWidth(v, ipk, level, h)

Width (grid spacing `h`) of the part of the profile `v` around index `ipk`
that is at or above `level`, with linear interpolation of the two
crossings. `NaN` if `v` does not fall below `level` on both sides.
"""
function levelWidth(v, ipk, level, h)
    i = ipk
    while i > 1 && v[i-1] >= level
        i -= 1
    end
    j = ipk
    while j < length(v) && v[j+1] >= level
        j += 1
    end
    (i == 1 || j == length(v)) && return oftype(h, NaN)
    xl = i - (v[i] - level) / (v[i] - v[i-1])
    xr = j + (v[j] - level) / (v[j] - v[j+1])
    (xr - xl) * h
end

"""
    spotMetrics(psf::PSFResult; fractions = (0.5, 0.8, 0.9), center = :centroid)

Size and shape of the wavefront-limited spot, as a NamedTuple (lengths in
mm; positions are along `psf.axisX`/`psf.axisY`):

- `peak`, `peakX`, `peakY` -- largest intensity (the peak Strehl) and its position
- `centroidX`, `centroidY`, `rmsRadius` -- first and second intensity moments
  over the window (the window should contain essentially all the energy, as
  moments are sensitive to the tails), and `semiMajor`, `semiMinor`
  (`2σ` along the principal axes) with `angle` (rad, major axis from `axisX`)
- `fwhmX`, `fwhmY` -- full width at half maximum along the grid axes through
  the peak; `fwhmArea` -- diameter of the disc with the area where
  `intensity ≥ peak/2`; `diameterE2` -- the same for `peak/e²`
- `ee` -- encircled-energy radii for the energy `fractions`, about `center`
  (`:centroid`, `:peak` or `:chief`, the reference point at the origin),
  normalized by `psf.totalEnergy` so the window size does not matter;
  `NaN` where the window does not reach that energy
"""
function spotMetrics(psf::PSFResult{T}; fractions = (0.5, 0.8, 0.9), center::Symbol = :centroid) where T
    I = psf.intensity
    h = psf.pitch
    dA = h^2
    xs, ys = psf.x, psf.y
    ipk = argmax(I)
    peak = I[ipk]
    sumI = sum(I) * dA
    cx = sum(I[i, j] * xs[i] for i in eachindex(xs), j in eachindex(ys)) * dA / sumI
    cy = sum(I[i, j] * ys[j] for i in eachindex(xs), j in eachindex(ys)) * dA / sumI
    sxx = sum(I[i, j] * (xs[i] - cx)^2 for i in eachindex(xs), j in eachindex(ys)) * dA / sumI
    syy = sum(I[i, j] * (ys[j] - cy)^2 for i in eachindex(xs), j in eachindex(ys)) * dA / sumI
    sxy = sum(I[i, j] * (xs[i] - cx) * (ys[j] - cy) for i in eachindex(xs), j in eachindex(ys)) * dA / sumI
    ev = eigen(Symmetric([sxx sxy; sxy syy]))
    major = ev.vectors[:, 2]
    c0 = center === :centroid ? (cx, cy) : center === :peak ? (xs[ipk[1]], ys[ipk[2]]) :
         center === :chief ? (zero(T), zero(T)) : throw(ArgumentError("center must be :centroid, :peak or :chief"))
    d = [hypot(xs[i] - c0[1], ys[j] - c0[2]) for i in eachindex(xs), j in eachindex(ys)]
    order = sortperm(vec(d))
    dsort = vec(d)[order]
    cum = cumsum(vec(I)[order]) .* dA ./ psf.totalEnergy
    ee = map(fractions) do f
        k = findfirst(>=(f), cum)
        k === nothing && return T(NaN)
        k == 1 && return dsort[1]
        dsort[k-1] + (f - cum[k-1]) / (cum[k] - cum[k-1]) * (dsort[k] - dsort[k-1])
    end
    (peak = peak, peakX = xs[ipk[1]], peakY = ys[ipk[2]], centroidX = cx, centroidY = cy,
     rmsRadius = sqrt(sxx + syy), semiMajor = 2sqrt(ev.values[2]), semiMinor = 2sqrt(max(ev.values[1], 0)),
     angle = atan(major[2], major[1]),
     fwhmX = levelWidth(I[:, ipk[2]], ipk[1], peak / 2, h), fwhmY = levelWidth(I[ipk[1], :], ipk[2], peak / 2, h),
     fwhmArea = 2sqrt(count(>=(peak / 2), I) * dA / π), diameterE2 = 2sqrt(count(>=(peak / exp(2)), I) * dA / π),
     ee = collect(ee))
end

"""
    geometricSpot(r::WavefrontResult; opd = r.opd)

Geometric (ray) spot at the reference point, for comparison with the
wavefront-limited spot: where each valid ray crosses the plane through the
reference point normal to the chief ray's image-space direction. Returns
`(rmsRadius, offsets)`: the RMS distance (mm) of the crossings from their
centroid, and the `(x, y)` offsets of each valid ray from the reference
point along [`defaultImageAxes`](@ref). `opd` only selects which samples
are valid (those with finite OPD). Unweighted.
"""
function geometricSpot(r::WavefrontResult{T}; opd = r.opd) where T
    e1, e2 = defaultImageAxes(r.chiefDir)
    dc = r.chiefDir
    offsets = Tuple{T,T}[]
    for i in eachindex(opd)
        (r.mask[i] && isfinite(opd[i])) || continue
        b, d = r.lastBase[i], r.lastDir[i]
        t = sum(dc .* (r.imagePoint .- b)) / sum(dc .* d)
        p = b .+ t .* d .- r.imagePoint
        push!(offsets, (sum(p .* e1), sum(p .* e2)))
    end
    mx = sum(first, offsets) / length(offsets)
    my = sum(last, offsets) / length(offsets)
    sqrt(sum((o[1] - mx)^2 + (o[2] - my)^2 for o in offsets) / length(offsets)), offsets
end
