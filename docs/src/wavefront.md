# Wavefront error analysis

`wavefrontMap` computes the wavefront error (OPD) of a system for one source
and wavelength, relative to the chief ray from that source. Everything else
-- statistics, best focus, the PSF and the plots -- takes the
[`WavefrontResult`](@ref) it returns, so nothing is re-traced. The design,
its decisions and what is still open are in `plans/wavefront_analysis.md`.

## Example

A concave paraboloid mirror (focal length 100 mm) with the image plane 0.5 mm
short of the focus, lit by a collimated beam:

```julia
using OpticTrace, GeometryBasics

mirror = reflectConic("mirror", Point3(0., 0., 100.), Vec3(0., 0., 1.), 1., 1., -1 / 200, 0.0, 10., "c")
image  = refractSphere("img", Point3(0., 0., 0.5), Vec3(0., 0., 1.), 1., 1., 0.0, 25., "c")
geo = Vector{AbstractSurface}([mirror, image])

r = wavefrontMap(geo, CollimatedField(Vec3(0., 0., 1.)), 0.5; geoWavelength = 0.5)

st = wavefrontStats(r)          # st.rms ≈ 1.44, st.pv ≈ 5.0 waves, st.strehl ≈ 4e-7
δ, rms, opd = bestFocus(r)      # δ = 0.5 mm (move the reference point toward the focus), rms ≈ 0
psf = wavefrontPSF(r; pixels = 257, halfWidth = 0.07)
spotMetrics(psf)                # peak (Strehl), FWHM, encircled energy, RMS radius, ellipse

using GLMakie
wavefrontFigure(r)              # heatmap, 3D surface, x/y cuts
psfFigure(psf)
```

## What it does

For each point of an odd `N × N` grid over the unit pupil disk, a ray is
launched (`StopCenterAim`; the pupil surface must be `geo[1]`), traced through
the surfaces before `imageSurface`, and its optical path -- `Σ nIn·delta` plus
any [`surfacePhase`](@ref) -- to a reference sphere is compared with the
chief ray's. The reference sphere is centered on the point where the chief ray
meets the image surface and passes through the real-ray exit pupil
([`OpticTrace.findExitPupil`](@ref)), in the actual image-space refractive index; if the
exit pupil is at infinity it becomes a plane (`refRadius == Inf`).

- **Sources:** [`ObjectPoint`](@ref) (finite point) or
  [`CollimatedField`](@ref) (object at infinity).
- **Wavelength:** the refractive indices in `geo` are fixed numbers set when
  `geo` was built, so the `λ` argument (μm) only converts OPD to waves and
  feeds the phase hooks. Pass `geoWavelength` (the wavelength `geo` was built
  for); a warning is issued when it is missing or differs from `λ`.
- **Phase beyond `nIn·delta`:** add a [`bendPhase`](@ref) method for a new
  `AbstractBendType` (gratings, diffractive optics, ...). It returns an optical
  path length in mm.
- **Apodization / input beam:** `apodization = GaussianPupil(w)`,
  `FunctionPupil(f)` or `AnglePupil(f)`; the amplitude becomes the ray weight
  used by the statistics and the PSF.
- **Unsupported:** a `CDiffuser` in `geo` (random scatter) makes OPD
  meaningless, so `wavefrontMap` warns and returns `nothing`. A
  `ParaxialLensT` adds no phase (`TODO.md` Code issues #41), so OPD through it
  is wrong, and a warning says so.
- `OptSurface`s do not clip rays at their aperture (`TODO.md` Code issues
  #42), so the pupil mask only reflects misses, total internal reflection and
  `ModelSurface` clipping.

## Results and statistics

- [`wavefrontRMS`](@ref), [`wavefrontPV`](@ref), [`strehlMarechal`](@ref),
  [`strehlRatio`](@ref), [`wavefrontStats`](@ref) -- piston removed, weighted
  by the pupil amplitude.
- [`refocusOPD`](@ref) / [`bestFocus`](@ref) -- move the reference point along
  the chief ray without re-tracing; best focus minimizes the RMS.
- [`wavefrontRemove`](@ref) -- subtract piston, tilt and focus for display or
  for statistics that should exclude them.

## PSF and spot

[`wavefrontPSF`](@ref) sums the rays directly (scalar Debye, equal quadrature
weight per ray): `U(x) = Σ w e^{i2πW} e^{i2π s·x/λ'}`. The output pitch and
window are free; the PSF is normalized so its peak is the Strehl ratio. A
warning is issued when the pitch is coarser than Nyquist or the window reaches
past half the aliasing period set by the pupil sampling.
[`spotMetrics`](@ref) reports peak, centroid, RMS radius, the second-moment
ellipse, FWHM, and encircled-energy radii; [`geometricSpot`](@ref) gives the
ray spot for comparison.

## Plots

[`wavefrontFigure`](@ref) and [`psfFigure`](@ref) build complete figures;
[`plotWavefrontHeatmap!`](@ref), [`plotWavefront3D!`](@ref),
[`plotWavefrontCuts!`](@ref), [`plotPSFHeatmap!`](@ref) and
[`plotPSFCuts!`](@ref) draw into an axis you supply.

## Sign convention

`OPD = OPL_ray − OPL_chief`: a positive OPD means the ray's path is longer
than the chief ray's. Whether Zemax uses the same sign has not been verified
(`TODO.md` Code issues #45); the one constant `WAVEFRONT_OPD_SIGN` flips it.

## Cross-check

[`rmsWavefrontLegacy`](@ref) recomputes the RMS with the older `opdRel`
machinery on equal-area ring samples. It agrees with `wavefrontMap` only
inside the legacy assumptions (on-axis, finite object, air image space).
