# Plots of the wavefront-error analysis (src/wavefront.jl). Every function takes
# the analysis result (a WavefrontResult or PSFResult), never a geometry:
# nothing is re-traced here.

export plotWavefrontHeatmap!, plotWavefront3D!, plotWavefrontCuts!, wavefrontFigure
export plotPSFHeatmap!, plotPSFCuts!, psfFigure

"""
    MIN_WAVEFRONT_COLOR_RANGE

Smallest half-width (waves) of the symmetric color range used by the
wavefront plots, so that a wavefront flat to better than a milliwave (e.g.
rounding noise left after removing tilt) is not stretched to fill the
range.
"""
const MIN_WAVEFRONT_COLOR_RANGE = 1e-3

"""
    plotWavefrontHeatmap!(ax, r::WavefrontResult; remove = :piston, weighted = true,
                          opd = r.opd, colormap = :balance)

Heatmap of the wavefront error (waves) over the pupil `(u, v)` into the
`Axis` `ax`, with a color range symmetric about zero and at least
`MIN_WAVEFRONT_COLOR_RANGE` wide. Masked points are
left blank. `remove` is as in [`wavefrontRemove`](@ref). Returns the plot
(for a `Colorbar`).
"""
function plotWavefrontHeatmap!(ax, r::WavefrontResult; remove::Symbol = :piston, weighted::Bool = true,
        opd = r.opd, colormap = :balance)
    z = wavefrontRemove(r; remove, weighted, opd)
    lim = max(maximum(abs, filter(isfinite, z); init = 0.0), MIN_WAVEFRONT_COLOR_RANGE)
    Makie.heatmap!(ax, r.u, r.v, z; colormap, colorrange = (-lim, lim), nan_color = :transparent)
end

"""
    plotWavefront3D!(ax, r::WavefrontResult; remove = :piston, weighted = true,
                     opd = r.opd, colormap = :balance)

The wavefront error (waves) as a surface over the pupil `(u, v)`, into
the `Axis3` `ax`. See [`plotWavefrontHeatmap!`](@ref). Returns the plot.
"""
function plotWavefront3D!(ax, r::WavefrontResult; remove::Symbol = :piston, weighted::Bool = true,
        opd = r.opd, colormap = :balance)
    z = wavefrontRemove(r; remove, weighted, opd)
    lim = max(maximum(abs, filter(isfinite, z); init = 0.0), MIN_WAVEFRONT_COLOR_RANGE)
    Makie.surface!(ax, r.u, r.v, z; colormap, colorrange = (-lim, lim))
end

"""
    plotWavefrontCuts!(ax, r::WavefrontResult; remove = :piston, weighted = true, opd = r.opd)

The wavefront error (waves) along the two pupil cross-sections through the
center, into the `Axis` `ax`: `v = 0` against `u` (dashed, "x cut") and
`u = 0` against `v` (solid, "y cut"). The grid always has a center row and
column (`N` is odd). Returns the two line plots `(xcut, ycut)`.
"""
function plotWavefrontCuts!(ax, r::WavefrontResult; remove::Symbol = :piston, weighted::Bool = true, opd = r.opd)
    z = wavefrontRemove(r; remove, weighted, opd)
    mid = (length(r.u) + 1) ÷ 2
    xcut = Makie.lines!(ax, r.u, z[:, mid]; linestyle = :dash, label = "x cut (v = 0)")
    ycut = Makie.lines!(ax, r.v, z[mid, :]; label = "y cut (u = 0)")
    xcut, ycut
end

"""
    wavefrontFigure(r::WavefrontResult; remove = :piston, weighted = true, opd = r.opd,
                    size = (1000, 800))

A `Figure` with the wavefront heatmap (with colorbar), the 3D surface, the
x and y cuts, and a title line with the RMS, peak-to-valley and Strehl
ratio of what is shown (after `remove`). Returns the figure.
"""
function wavefrontFigure(r::WavefrontResult; remove::Symbol = :piston, weighted::Bool = true, opd = r.opd,
        size = (1000, 800))
    z = wavefrontRemove(r; remove, weighted, opd)
    fig = Figure(; size)
    rms = wavefrontRMS(r; weighted, opd = z)
    pv = wavefrontPV(r; opd = z)
    st = strehlRatio(r; weighted, opd = z)
    Label(fig[0, 1:3], "OPD at λ = $(r.λ) μm (waves, $remove removed):  RMS = $(round(rms; sigdigits = 3)),  PV = $(round(pv; sigdigits = 3)),  Strehl = $(round(st; sigdigits = 3))";
        tellwidth = false)
    ax1 = Axis(fig[1, 1]; xlabel = "u", ylabel = "v", aspect = DataAspect(), title = "Wavefront error")
    hm = plotWavefrontHeatmap!(ax1, r; remove, weighted, opd)
    Colorbar(fig[1, 2], hm)
    ax2 = Axis3(fig[1, 3]; xlabel = "u", ylabel = "v", zlabel = "waves", title = "Surface")
    plotWavefront3D!(ax2, r; remove, weighted, opd)
    ax3 = Axis(fig[2, 1:3]; xlabel = "normalized pupil coordinate", ylabel = "OPD (waves)", title = "Cross-sections")
    plotWavefrontCuts!(ax3, r; remove, weighted, opd)
    axislegend(ax3)
    colgap!(fig.layout, 40)
    fig
end

"""
    psfScaled(psf::PSFResult, scale)

The PSF intensity as plotted: `:linear` unchanged, `:log` as `log10` with a
floor of `1e-5` (relative to the aberration-free peak).
"""
function psfScaled(psf::PSFResult, scale::Symbol)
    scale === :linear && return psf.intensity
    scale === :log && return log10.(max.(psf.intensity, 1e-5))
    throw(ArgumentError("scale must be :linear or :log, got :$scale"))
end

"""
    plotPSFHeatmap!(ax, psf::PSFResult; scale = :linear, colormap = :viridis)

Heatmap of the PSF (see [`psfScaled`](@ref) for `scale`) into the `Axis`
`ax`, with positions in μm. Returns the plot.
"""
function plotPSFHeatmap!(ax, psf::PSFResult; scale::Symbol = :linear, colormap = :viridis)
    Makie.heatmap!(ax, psf.x ./ LENGTH_TO_WAVELENGTH, psf.y ./ LENGTH_TO_WAVELENGTH, psfScaled(psf, scale); colormap)
end

"""
    plotPSFCuts!(ax, psf::PSFResult; scale = :linear)

Cuts through the PSF's peak along its two axes, into the `Axis` `ax`
(positions in μm): along `x` dashed, along `y` solid. Returns
`(xcut, ycut)`.
"""
function plotPSFCuts!(ax, psf::PSFResult; scale::Symbol = :linear)
    z = psfScaled(psf, scale)
    ipk = argmax(psf.intensity)
    x = psf.x ./ LENGTH_TO_WAVELENGTH
    y = psf.y ./ LENGTH_TO_WAVELENGTH
    xcut = Makie.lines!(ax, x, z[:, ipk[2]]; linestyle = :dash, label = "along x")
    ycut = Makie.lines!(ax, y, z[ipk[1], :]; label = "along y")
    xcut, ycut
end

"""
    psfFigure(psf::PSFResult; scale = :linear, size = (1000, 500))

A `Figure` with the PSF heatmap (with colorbar), its peak cuts, and a title
line with the peak Strehl, the x and y FWHM and the radius enclosing 80% of
the energy (from [`spotMetrics`](@ref)). Returns the figure.
"""
function psfFigure(psf::PSFResult; scale::Symbol = :linear, size = (1000, 500))
    m = spotMetrics(psf; fractions = (0.8,))
    um(x) = round(x / LENGTH_TO_WAVELENGTH; sigdigits = 3)
    fig = Figure(; size)
    Label(fig[0, 1:3], "PSF at λ = $(psf.λ) μm:  Strehl = $(round(m.peak; sigdigits = 3)),  FWHM = $(um(m.fwhmX)) × $(um(m.fwhmY)) μm,  80% EE radius = $(um(m.ee[1])) μm";
        tellwidth = false)
    ax1 = Axis(fig[1, 1]; xlabel = "x (μm)", ylabel = "y (μm)", aspect = DataAspect(), title = "PSF ($scale)")
    hm = plotPSFHeatmap!(ax1, psf; scale)
    Colorbar(fig[1, 2], hm)
    ax2 = Axis(fig[1, 3]; xlabel = "position (μm)", ylabel = scale === :log ? "log10 intensity" : "intensity", title = "Cuts through the peak")
    plotPSFCuts!(ax2, psf; scale)
    axislegend(ax2)
    colgap!(fig.layout, 40)
    fig
end
