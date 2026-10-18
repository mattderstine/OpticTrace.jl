# Plan: wavefront error analysis (`wavefrontMap` and friends)

Status: Step 0 in progress. Supersedes the quick `opdRel`-based path in
`src/plotting.jl` (which stays as a legacy path).

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
| Reference | Real-ray exit pupil (chief + 4 offset rays, least-squares point, 3D, local frame); sphere centred on the chief image point; closed-form signed ray-sphere intersection in image space using the actual image-space index; plane reference when the pupil is at infinity. |
| Wavelength | `geo` already has its indices instantiated at one wavelength. `λ` is used for the phase hooks, the waves conversion, and PSF scaling, not for refractive index. |

## Data structures

- `FieldPoint{T}` — object position (finite) or direction (collimated), mirroring `OpticalSystem.fieldType`.
- `AbstractRayAimer` — chief-ray-aiming hook. `aimRay(aimer, system, field, λ, u, v) -> (status, Ray)`. Default `StopCenterAim` reproduces today's behavior; a later `IterativeAim` solves for launch direction. The same aimer launches the chief ray and every grid ray.
- `AbstractPupilAmplitude` — apodization / input-beam hook: `UniformPupil`, `GaussianPupil`, `FunctionPupil(f(u,v))`, plus an angle-dependent variant. Final weight = apodization × product of surface `AmpData.trans`.
- `WavefrontResult{T}` — grid `(u,v)`, mask, OPD in waves, weights, physical exit-pupil coordinates per ray, chief-ray data, reference-sphere data (centre, radius, image index), `λ`, `geoWavelength`, field.

Entry points:
- `wavefrontMap(system, field, λ; N=65, aimer, apodization, geoWavelength, ...) -> WavefrontResult`
- `wavefrontMap(system, fields, λs) -> Matrix{WavefrontResult}`

Files: `src/wavefront.jl` (analysis, no Makie), included after `characterization.jl`; plot functions in `plotting.jl` or a new `wavefront_plots.jl` after it.

## Core algorithm

1. Chief ray via the aimer at `(u,v) = (0,0)`.
2. Source wavefront: finite objects start in phase at the point source; collimated fields start in phase on a plane normal to the chief direction at the entrance pupil.
3. Per-ray `OPL = Σ nᵢ·δᵢ + Σ surfacePhase(...)`.
4. Exit pupil: chief ray plus four offset rays (±x, ±y); least-squares point on the chief ray's image-space line. If beyond a threshold or the rays are nearly parallel, switch to a plane reference.
5. Reference sphere centred on the chief image point, radius = path distance to the exit pupil point; closed-form signed intersection of each ray's last segment, in the actual image-space index.
6. Subtract the chief ray's OPL (zero OPD at pupil centre). Mean-piston removal for statistics is separate; both are stored.

## Statistics, best focus, Strehl, spot

- RMS: weighted, mean-removed. PV: max − min over unmasked samples, optional weight threshold.
- Best focus: 1D minimization (`Optim`) of weighted variance over a shift δ of the sphere centre along the chief direction, re-evaluated analytically from stored last-ray data (no re-trace). Optional lateral shift (tilt removal).
- Strehl: Maréchal `exp(-(2πσ)²)` and the FFT value (PSF peak / diffraction-limited peak).
- Spot / PSF: pupil function `A·exp(i2πW)`, zero-padded (default 4×), FFT; image pitch `(λ/n)·R/(M·Δp)`; encircled-energy radii, x/y FWHM, second-moment ellipse. Plane reference: angle only.
- Sampling sufficiency: `N = :auto` starts at 33 and doubles until RMS changes < 0.1% and the largest neighbour step < ~0.25 wave; warns at a cap.

## Plots

`plotWavefront3D!`, `plotWavefrontHeatmap!`, `plotWavefrontCuts!`, and a combined `wavefrontFigure`; all take a `WavefrontResult`. Masked points are NaN. Display option `remove = :none | :piston | :pistonTilt`.

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
- Polychromatic / multi-field aggregation: needs a fixed image-plane pitch across wavelengths (resampling or chirp-z DFT), a reference-wavelength choice for lateral colour, and the wavelength re-instantiation task below.
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
