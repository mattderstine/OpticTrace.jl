# Plan: wavefront error analysis (`wavefrontMap` and friends)

Status: Step 0 done (Zemax sign convention still open); Steps 1-4 done. The
"Step N status" sections at the end record what was actually built and where it
differs from the design text above; where they disagree, the status sections
win. Supersedes the quick `opdRel`-based path in `src/plotting.jl` (which
stays as a legacy path).

## Goal

Compute the wavefront error of a system from an object point, relative to the
chief ray from that point, and derive from it: RMS, peak-to-valley, best focus
for RMS, plots (3D surface, heatmap, X/Y cuts), Strehl ratio, and the
wavefront-limited spot size and shape (PSF).

## Problems in the existing `opdRel`/`plotOPD!` path to fix

- No chief-ray piston removal / no mean-piston removal for statistics.
- Reference sphere vertex forced to `(0, 0, zExitPupil)`; radius from a
  z-difference only. Wrong for off-axis fields, tilted or decentred pupils.
- Exit pupil from a small-angle (1e-4 rad) ray, global-Z, assumes `dirb[2] != 0`;
  silently falls back to the prior surface when the pupil is at infinity.
- Reference sphere uses `refIndexDefault` instead of the image-space index.
- Only finite-object, only the `geo[1]` "stop" convention.

## Design decisions (agreed)

| Topic | Decision |
|---|---|
| Pupil sampling | Square Cartesian grid in normalized pupil coords, odd `N` (default 65), masked to the pupil. |
| Legacy cross-check | Separate function `rmsWavefrontLegacy(...)` using `opdRel` over area-uniform ring samples, piston mean removed. Test oracle only; valid only inside the legacy envelope (on-axis, finite object, air image space). |
| Extra phase | OPL = `Σ nᵢ·δᵢ + Σ surfacePhase(...)`; `surfacePhase` is a dispatch hook returning an OPL-equivalent length, default 0. Applied in the accumulation function, not inside `traceSurf` (leaves `Trace` and the hot path alone). |
| Result type | New `WavefrontResult{T}` struct from a new `wavefrontMap`; `opdRel` and the three `plotOPD!` methods stay. Stats and plots take the struct. |
| Reference | Real-ray exit pupil (4 probe rays from the pupil center tilted about the chief direction, least-squares point on the chief ray, 3D, global frame); sphere centred on the chief image point; closed-form signed ray-sphere intersection in image space using the actual image-space index; plane reference when the pupil is at infinity. |
| Wavelength | `geo` already has its indices instantiated at one wavelength. `λ` is used for the phase hooks, the waves conversion, and PSF scaling, not for refractive index. |

## Data structures

- `ObjectPoint{T}` (finite point source) and `CollimatedField{T}` (object at infinity), subtypes of `AbstractFieldPoint{T}`.
- `AbstractRayAimer` — chief-ray-aiming hook. `aimRay(aimer, geo, field, u, v) -> (status, Ray)`. Default `StopCenterAim` reproduces today's behavior; a later `IterativeAim` solves for launch direction. The same aimer launches the chief ray and every grid ray.
- `AbstractPupilAmplitude` — apodization / input-beam hook: `UniformPupil`, `GaussianPupil`, `FunctionPupil(f(u,v))`, plus an angle-dependent variant. Final weight = apodization × product of surface `AmpData.trans`.
- `WavefrontResult{T}` — grid `(u,v)`, mask, OPD in waves, weights, physical exit-pupil coordinates per ray, chief-ray data, reference-sphere data (centre, radius, image index), `λ`, `geoWavelength`, field.

Entry points (as implemented):
- `wavefrontMap(geo, field, λ; N=65, aimer, apodization, imageSurface, geoWavelength, ...) -> WavefrontResult` (or `nothing` after a warning)
- Not implemented: methods taking an `OpticalSystem`/`ExtendedGeometry`, and `wavefrontMap(system, fields, λs) -> Matrix{WavefrontResult}` (see Deferred and `TODO.md`).

Files: `src/wavefront.jl` (analysis, no Makie), included after `characterization.jl`; plot functions in `src/wavefront_plots.jl`, included after it.

## Core algorithm

1. Chief ray via the aimer at `(u,v) = (0,0)`.
2. Source wavefront: finite objects start in phase at the point source; collimated fields start in phase on a plane normal to the chief direction at the entrance pupil.
3. Per-ray `OPL = Σ nᵢ·δᵢ + Σ surfacePhase(...)`.
4. Exit pupil: four probe rays launched **from the pupil center**, tilted ±ε about the chief ray's launch direction in two perpendicular planes (not rays from the object point toward other pupil points, which converge at the image point, not the pupil); least-squares point on the chief ray's image-space line nearest the probes' image-space lines = image of the stop center. If the probes emerge parallel to the chief ray (RMS sine of the angle < `pupilInfinityRatio·ε`), the pupil is at infinity and the reference is a plane.
5. Reference sphere centred on the chief image point, radius = path distance to the exit pupil point; closed-form signed intersection of each ray's last segment, in the actual image-space index.
6. Subtract the chief ray's OPL (zero OPD at pupil centre). Mean-piston removal for statistics is separate; both are stored.

## Statistics, best focus, Strehl, spot

- RMS: weighted, mean-removed. PV: max − min over unmasked samples, optional weight threshold.
- Best focus: 1D minimization (`Optim`) of weighted variance over a shift δ of the sphere centre along the chief direction, re-evaluated analytically from stored last-ray data (no re-trace). Lateral shift (tilt removal) is not optimized.
- Strehl: Maréchal `exp(-(2πσ)²)`; the reference-point value `|Σ w e^{i2πW}|²/(Σ w)²`; and the PSF peak (`PSFResult.strehl`).
- Spot / PSF: a direct scalar-Debye ray sum (two matrix products, no FFT, no extra dependency) with the output pitch and window chosen freely; encircled-energy radii, x/y FWHM, second-moment ellipse; works unchanged for a plane reference. See the Step 3 status.
- Sampling sufficiency: `N = :auto` (start at 33, double until RMS changes < 0.1% and the largest neighbour step < ~0.25 wave) is **not implemented**; `N` is a plain keyword (default 65).

## Plots

`plotWavefront3D!`, `plotWavefrontHeatmap!`, `plotWavefrontCuts!`, and a combined `wavefrontFigure`; all take a `WavefrontResult`. Masked points are NaN. Display option `remove = :none | :piston | :pistonTilt | :pistonTiltFocus`.

## Wavelength handling

- `geo` carries indices fixed at build time. The `λ` argument to `wavefrontMap` is the evaluation wavelength: it converts OPD to waves, feeds `surfacePhase`, and sets the PSF scale. It does not change any index.
- Optional `geoWavelength` keyword (default from the system object when available); stored in the result; a warning is issued when it differs from `λ`. If unknown, the result says so rather than guessing.
- Multi-wavelength calls on one `geo` are physically meaningful only for phase-object contributions until the future task below exists.

## Tests (`test/wavefront.jl`, wired into `runtests.jl`)

- Paraboloid mirror, collimated: OPD ≈ 0.
- Defocus: RMS = a/(2√3), `a` the peak OPD.
- Pure W040: RMS = W040/(6√5) at best focus, best-focus shift = ρ² coefficient −W040.
- Rigid rotation/translation of the whole system leaves the result unchanged.
- Non-air image space.
- Grating phase hook; uniform apodization equals none.
- Maréchal and FFT Strehl agree for small σ.
- Zemax comparison by hand against local samples; constants added to tests only with approval (no sample files copied into the repo).
- Legacy cross-check inside its envelope.
- Per `CLAUDE.md`: per-method docstrings (verify with `@doc`), `TODO.md` entries without renumbering, update `CLAUDE.md` project structure.

## Zernike fit (requirements only)

- `WavefrontResult` stores normalized pupil coordinates, mask and weights, and a definition of the unit disk (the stop semi-diameter, not the extent of the unmasked points).
- Pick a convention: Noll, ANSI/OSA, or Zemax Fringe/Standard.
- Weighted linear least-squares of sampled OPD onto the basis up to a user-set order; samples must greatly outnumber terms.
- Outputs: coefficients in waves, residual map, RMS from coefficients, removal of chosen terms (feeds `remove=`, cross-checks best focus via the defocus term).
- Vignetted or elliptical pupils break orthogonality; needs a documented policy.

## Deferred

- Iterative ray aiming (hook exists).
- Polychromatic / multi-field aggregation: needs a reference-wavelength choice for lateral colour and the wavelength re-instantiation task below. (The fixed-image-pitch concern is gone: `wavefrontPSF` takes any output pitch.)
- `N = :auto` sampling sufficiency; `wavefrontMap` for `OpticalSystem`/`ExtendedGeometry`; a pupil surface other than `geo[1]`.
- Radiometric (Jacobian) amplitude weighting.
- Resampling to a uniform exit-pupil grid.
- **Re-instantiating `geo` at a different wavelength.** Indices are baked into surfaces as numbers. Zemax import and catalogs hold `wavelength -> index` functions but nothing connects them to built surfaces. Needs a design (e.g. `geoAtWavelength(system, λ)`, or surfaces that keep an index function), preserving type parameters for `ForwardDiff` and the `traceGeometry!` allocation behavior.

## Risks

- Pupil distortion: grid is uniform in the stop pupil but the FFT assumes uniform exit-pupil sampling. Version one stores exit-pupil coordinates, reports a distortion metric, and does not resample.
- Extending `AmpData` for complex coating phase touches a core struct; treat as its own change.
- Allocations: use `traceGeometryRel!` with per-thread buffers and `Threads.@threads` over grid rows.

## Step 0: verify conventions before writing code

1. `delta`/`nIn` pairing: the `Trace` docstring says both are "leaving surface"; `modFunc` says `nIn` is the index before the surface. Check with a plane-parallel plate whether `Σ nIn·delta` is the true OPL, including first and last segments.
2. Indices in `geo` are fixed numbers: confirm every index-bearing surface holds numeric indices and nothing in the trace path looks up an index from a wavelength.
3. Where the geometry's own wavelength is recorded (`ExtendedGeometry.wavelength`, `OpticalSystem.wavelengths`/`primaryWavelengthIndex`); decides how `geoWavelength` gets into the API.
4. Signed propagation: does the tracer accept a negative `delta` (virtual exit pupil behind the last surface)?
5. Sign convention for OPD: define explicitly; validate against a Zemax defocus case using local sample files.

### Step 0 findings

Checked 2026-10-17 by reading `src/tracing.jl`/`src/lens_definitions.jl` and
running a scratch trace (plane-parallel plate, tilted ray, sphere, mirror).

1. **`delta`/`nIn` pairing — resolved, consistent.** `Trace.delta` is the
   distance travelled from the previous point *to* this surface, and
   `Trace.nIn` is the index of the medium that segment travelled in
   (`modFunc` returns the pre-surface index). So `Σ nIn·delta` is the true OPL
   from the start ray to the last surface. Plate test: start on surface 1,
   glass 10 mm (n=1.5), air 30 mm gave 45.0 = 0 + 15 + 30; tilted ray gave
   46.60 as expected. The `Trace` docstring ("leaving surface" for both) is
   stale and should be corrected.
   - **Consequence for the design:** `traceGeometryRel` takes the start ray in
     the *local frame of `geo[1]`* and the first `Trace` entry has
     `delta = 0`, `nIn = 1`. The object-to-first-surface segment is only
     counted if the ray starts at the object point, in front of `geo[1]`. The
     new launcher must start rays at the object point (finite) or at the
     source-wavefront plane (collimated) so that segment is included.
     Note `trace[1].nIn = 1.0` is a placeholder, not the object-space index;
     the first real segment's index is the first surface's `refIndexIn`.
2. **Indices in `geo` are fixed numbers — confirmed.** `DielectricT`,
   `MirrorR`, `ParaxialLensT`, `CDiffuser`, `NoBendIndex` all store numeric
   `refIndexIn`/`refIndexOut`; `modFunc` reads them directly. Nothing in the
   trace path looks up an index from a wavelength.
3. **Geometry wavelength is not recorded — finding.**
   `OpticalSystem` has `wavelengths`/`primaryWavelengthIndex` (the system's
   design wavelengths) but **not** the wavelength `geo` was instantiated at:
   `readZemaxSystem(...; wavelength = 0.5)` bakes indices in at a keyword
   default (0.5 μm) and discards it. `ExtendedGeometry.wavelength` is the
   array handed to `funcSetup`, so it is a reasonable source for that type
   only. Decision: `geoWavelength` is an explicit keyword, defaulting from
   `ExtendedGeometry.wavelength[1]` when given one, and unknown (with a
   warning) for `OpticalSystem`/bare `Vector{AbstractSurface}` until an
   `instantiatedWavelength` field is added to `OpticalSystem` (small change in
   `readZemaxSystem`; to be done with the re-instantiation task).
4. **Signed propagation — confirmed.** `deltaToSurf` picks one root
   (`(g - sqrt(q))/C`) with no positivity check and `traceSurf` accepts it.
   A sphere start behind the vertex returned `delta = -29.8`, status 0.
   Negative `delta` is therefore legal in the tracer, though the plan's
   closed-form image-space sphere intersection still avoids relying on it.
   Caveat: the root choice is the near-side root of the surface; confirm the
   branch for rays that start beyond the sphere centre when testing the
   reference-sphere routine.
5. **OPD sign convention — open, owner: user.** Our convention (inherited from
   `opdRel`): `OPD = OPL_ray − OPL_chief`. There is no Zemax license
   available to check it, so the user will look for the correct convention
   elsewhere (documentation / literature). Until then the sign is a single
   documented constant in the code (easy to flip), tests use sign-independent
   quantities (RMS, PV, Strehl, |best-focus shift|) and sign-dependent tests
   are marked `@test_broken`/deferred. Zemax-comparison tests are limited to
   sign-independent quantities.

### Additional findings that affect the design

- **`ParaxialLensT` (ideal thin lens) adds no phase.** `modFunc` only bends
  the ray. A true ideal lens needs a phase of `-r²/(2f)` (in OPL terms), so a
  `surfacePhase` method for `ParaxialLensT` is required for correct OPD
  through paraxial/ideal-lens systems, including Zemax `PARAXIAL` imports.
- **`CDiffuser` is random** (samples a random direction), so it is
  incompatible with OPD calculations. **Decision:** `wavefrontMap` checks
  `geo` for any surface whose `mod` is a `CDiffuser` before tracing; if one is
  present it emits `@warn` naming the surface and exits, returning `nothing`
  (no partial result). The legacy `opdRel` path is left unchanged.
- **`ParaxialLensT` phase** is tracked as `TODO.md` Code issues #41 (a
  `surfacePhase` method giving `-r²/(2f)`); implement with the phase hook.
- **Unit conversion in the existing `plotOPD!` is consistent.**
  `LENGTH_TO_WAVELENGTH = 0.001` (μm→mm by multiplying), and both methods
  divide an mm OPD by it correctly. Reuse the constant in the new code.
- **Scratch script** used for these checks:
  `step0.jl` in the session scratchpad (not in the repo). A regression test
  for the plate OPL (45.0) and the negative-delta case should be added to
  `test/wavefront.jl` as the first tests.


## Step 1 status (implemented)

Files: `src/wavefront.jl`, `test/wavefront.jl` (51 tests, wired into
`runtests.jl`). Implemented: `ObjectPoint`/`CollimatedField`,
`StopCenterAim` (the `AbstractRayAimer` hook), pupil amplitude types
(`UniformPupil`, `GaussianPupil`, `FunctionPupil`, `AnglePupil`),
`bendPhase`/`surfacePhase`/`pathLength`, real-ray exit pupil
(`findExitPupil`), closed-form reference sphere (`oplToReference`,
`referencePoint`), `WavefrontResult`, `wavefrontMap`, the `CDiffuser`
warn-and-exit, the `ParaxialLensT` warning, and the `geoWavelength` check.
`WAVEFRONT_OPD_SIGN` is the single sign constant.

Not yet done (later steps): RMS/PV/best-focus/Strehl, PSF, plots, the
`opdRel`-based cross-check function, Zernike, an `ExtendedGeometry`/
`OpticalSystem` method of `wavefrontMap`, a `Base.show` for
`WavefrontResult`, surface-transmission (`AmpData.trans`) in `weight`.

Verified by tests: zero OPD for a perfect paraboloid (collimated) and a
sphere with the point source at its center of curvature; defocus RMS =
peak/(2√3) within 3% and the OPD flips sign with the shift; OPD scales by n
in a medium of index n (image-space index is used); identical results after
a rigid rotation + translation of the whole system; plane reference
(exit pupil at infinity) and its large-sphere limit; user phase hook adds
exactly `κ·x`; apodization weights; the warn/exit paths.

Findings and decisions made while implementing:
- **A plane reference is a limit of the sphere, not a "plane at the
  chief point".** For an exit pupil at infinity the reference sphere is
  centred on the image point with `R → ∞`; the constant `−σ·R` cancels and
  the OPD is the path to the plane through the image point normal to each
  ray (`q = 0`). `oplToReference` computes this stably for any `R` as
  `τ + σ·q`, `q = h²/(R + √(R²−h²))`, so large finite and infinite `R`
  agree. (The earlier idea of a plane normal to the chief ray was wrong:
  it makes a converging beam look like huge defocus.)
- A collimated beam through a bare plane to an image plane gives a large
  OPD against a finite-radius sphere. That is physically right (a flat
  wave against a converging reference), not a bug.
- `OptSurface` does not clip rays at its aperture (only `ModelSurface`
  does), so the mask only reflects misses, TIR and `ModelSurface`
  clipping. Real vignetting by lens edges is not modelled yet.
- The pupil surface must be `geo[1]` (probe rays start at its center).
  Generalising this needs aiming through any surfaces in front of it.
- The chief ray's last-segment index (`nIn` of the image-surface hit) is
  the image-space index; for an asymmetric image surface its
  `refIndexIn`/`refIndexOut` must be oriented to match the ray direction.
- Numerics: pure-Float64 arrays inside a `T` parameter; `ForwardDiff`
  compatibility of `wavefrontMap` is untested.


## Step 2 status (implemented)

In `src/wavefront.jl`, tests in `test/wavefront.jl` ("wavefront statistics").
`wavefrontRMS`, `wavefrontPV`, `strehlMarechal`, `strehlRatio`,
`wavefrontStats` (all weighted by `weight`, mean-removed, optional `opd`
keyword to evaluate another OPD matrix); `refocusOPD(r, δ)` re-evaluates the
stored last segments for a reference point shifted `δ` mm along the chief
ray's image-space direction (no re-tracing; `δ = 0` reproduces `r.opd`);
`bestFocus(r)` minimizes the weighted RMS over `δ` with Brent (`Optim`) and
returns `(δ, rms, opd)`; `rmsWavefrontLegacy` is the `opdRel`-based
cross-check.

Verified: perfect systems give RMS 0 and Strehl 1; a known 0.5 mm image-plane
shift is recovered as `bestFocus` `δ = +0.5` (the reference point moves
along the propagation direction toward the focus), RMS 0 there; spherical
aberration (sphere mirror): best-focus RMS is 1/4 of the paraxial-focus
RMS and the best-focus shift matches the `W040` relation within 15%;
Maréchal and reference-point Strehl agree within 3% for RMS ~ 0.06 waves;
PV of a pure defocus matches `2√3·RMS`; `rmsWavefrontLegacy` agrees with
`wavefrontMap` within 5% for an on-axis finite-object lens inside the
legacy envelope.

Findings:
- **The legacy fan sampling was biased, not the new code.** Hexapolar rings
  with an outer ring at the edge radius and equal weight per ray overweight
  the edge; RMS of the same lens read 4.96 / 4.61 / 4.42 for 6 / 12 / 24
  rings against a converged 4.205 from `wavefrontMap` (N = 129). The legacy
  exit pupil (z = 1.4035) matched the new real-ray one exactly. The
  cross-check function therefore uses equal-area ring midpoints.
- The "FFT" Strehl is still deferred to Step 3 (PSF); `strehlRatio` is the
  value at the reference point, which equals the PSF peak only if the peak
  is there (no tilt).
- `bestFocus` only optimizes the longitudinal shift inside a bracket (default
  half the pupil-to-image distance); tilt (lateral shift) is not optimized.


## Step 3 status (implemented)

`wavefrontPSF`, `PSFResult`, `spotMetrics`, `geometricSpot` in
`src/wavefront.jl`; tests in `test/wavefront.jl` ("wavefront PSF and spot").

Method: scalar Debye sum over the traced rays,
`U(x) = Σ w e^{i2πW} e^{i2π s·x/λ'}` with `s` the ray's image-space
direction cosines along the PSF axes and `λ' = λ/nImage`. Chosen over an FFT
because (a) it needs no new dependency, (b) it uses each ray's actual
direction cosines, so it is correct for a pupil that maps non-uniformly onto
the exit pupil and for a plane reference (exit pupil at infinity) with no
special case, and (c) because the sum separates in `x` and `y` it is two
matrix products, so the output grid (pitch, size, window) is independent of
the pupil grid -- which also removes the "fixed image-plane pitch across
wavelengths" obstacle listed under Deferred for polychromatic PSFs. The
physical sign of the OPD is used (`WAVEFRONT_OPD_SIGN` is undone), so the
PSF orientation does not depend on the reporting convention.

Normalization: intensity is relative to the aberration-free peak of the same
pupil, so `maximum(intensity)` is the peak Strehl; `totalEnergy` is the
whole-plane integral from Parseval, so encircled energy does not depend on
the window. The pupil sampling makes the PSF periodic with
`period = λ'/sqrt(area per ray in direction-cosine space)`; the window must
stay inside half of it (a warning is issued).

Verified: the perfect parabola gives an Airy pattern (FWHM = 1.029 λ'/(2NA)
within 3%; 83.8% energy inside 0.61 λ'/NA within 4%; Strehl 1); a 0.5 mm
defocus gives Strehl < 0.1, a PSF RMS radius within 15% of the geometric
ray-spot RMS, and refocusing restores Strehl 1; a wavefront tilt (phase hook)
gives reference-point Strehl < 0.1 but a PSF peak > 0.95 displaced by
focal length × tilt; the on-axis PSF value equals `strehlRatio`.

Limits and notes:
- Equal quadrature weight per ray (no Jacobian of the pupil mapping, no
  obliquity factor). `pupilDistortion`/`pupilAnisotropy` quantify how far the
  pupil grid is from a uniform map; a warning is issued above 5% distortion.
- PSF axes default to the global x axis projected normal to the chief ray and
  its cross product with the chief direction; override with `axes`.
- Spot moments (`rmsRadius`, ellipse) come from the window and are
  tail-sensitive; the encircled-energy radii are window-independent.
- `localRaysHexapolar` (used by `spotDiagramHex`) still has the legacy
  edge-weighted ring pattern; left unchanged.


## Step 4 status (implemented)

`src/wavefront_plots.jl` (included after `wavefront.jl`; every function takes
an analysis result, nothing is re-traced) and `wavefrontRemove` in
`src/wavefront.jl`; tests in `test/wavefront.jl` ("wavefrontRemove",
"wavefront plots"). Rendered images of the figures are written to the
gitignored `test/plots/` by `saveTestFigure` for visual inspection.

- `wavefrontRemove(r; remove = :none | :piston | :pistonTilt | :pistonTiltFocus)`
  -- weighted least-squares removal of low-order terms in `(u, v)`; used by
  every wavefront plot for its `remove` keyword and available for
  statistics (RMS/PV/Strehl accept the result through their `opd` keyword).
- `plotWavefrontHeatmap!`, `plotWavefront3D!`, `plotWavefrontCuts!` (Axis /
  Axis3-level, return the plots), `wavefrontFigure` (heatmap + colorbar, 3D
  surface, x/y cuts, a title line with RMS/PV/Strehl of what is shown).
  Color range is symmetric about zero with a floor of
  `MIN_WAVEFRONT_COLOR_RANGE` (0.001 wave) so rounding noise is not stretched
  to fill the range.
- `plotPSFHeatmap!`, `plotPSFCuts!`, `psfFigure` (linear or `log10` scale
  with a 1e-5 floor; positions in μm; title line with Strehl, FWHM, 80% EE
  radius).

Not done: the surface plot leaves a ragged edge where the circular mask meets
the square grid (Makie drops NaN cells); a PSF-vs-focus (through-focus) plot;
a polychromatic/multi-field figure (needs the deferred wavelength
re-instantiation work).
