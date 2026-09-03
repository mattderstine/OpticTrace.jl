# Helper functions for testing

"""
define function to compute the normal vector using the gradient of the sag function
    This function takes a ray and a surface profile as input, and returns the normal vector at the intersection point of the ray with the surface.
The function first finds the intersection point of the ray with the surface using the deltaToSurf method, then finds the sag value at the intersection point, and finally finds the gradient of the sag value at the intersection point using an autodiff package. The normal vector is then obtained by normalizing the gradient vector.
"""
function normal_from_sag(ray::Ray, surfProfile)
    # find the intersection point of the ray with the surface using the deltaToSurf method
    delta = OpticTrace.deltaToSurf(ray, surfProfile)
    intersection_point = ray.base + delta * ray.dir
    # find the gradient of the sag value at the intersection point using an autodiff package
    grad_sag = ForwardDiff.gradient((x) -> sag(x[1], x[2], surfProfile), intersection_point[1:2])
    # the normal vector is then obtained by normalizing the gradient vector
    normal_vector = -normalize(Vec3(grad_sag[1], grad_sag[2], -1.0))
    return normal_vector
end

function normal_from_sag(intersection_point::Point3, surfProfile)
    # find the gradient of the sag value at the intersection point using an autodiff package
    grad_sag = ForwardDiff.gradient((x) -> sag(x[1], x[2], surfProfile), intersection_point[1:2])
    # the normal vector is then obtained by normalizing the gradient vector
    normal_vector = -normalize(Vec3(grad_sag[1], grad_sag[2], -1.0))
    return normal_vector
end

function normal_from_sag(x1::T, y1::T, surfProfile::S) where{T<:Real, S<:OpticTrace.AbstractSurfProfile{T}}
    f(x) = sag(x[1], x[2], surfProfile)

    x0 = [x1, y1]
    # find the gradient of the sag value at the intersection point using an autodiff package
    grad_sag = ForwardDiff.gradient(f, x0)
    # the normal vector is then obtained by normalizing the gradient vector
    normal_vector = -normalize(Vec3(grad_sag[1], grad_sag[2], -1.0))
    return normal_vector
end
function rprop(ray::Ray, delta::Float64)
    intersection_point = Point3(ray.base + delta * ray.dir)
    return intersection_point
end

"""
A trivial constant refractive-index function, for surface/lens-builder
tests that need a `riFunc` but shouldn't depend on any glass catalog.
"""
riFunc(λ) = 1.5

"""
Whether the real refractive-index glass catalog
(`OpticTrace.dirBaseRefractiveIndex`) is present on this machine. It's a
machine-local directory, not part of the repo (see `TODO.md`), so it
won't exist on a fresh checkout or in CI -- tests that need it should
check this and skip (not fail) when it's false, e.g.:

    if HAS_GLASS_CATALOG
        @testset "..." begin ... end
    else
        @info "Skipping ...: OpticTrace.dirBaseRefractiveIndex not found"
    end
"""
const HAS_GLASS_CATALOG = isdir(OpticTrace.dirBaseRefractiveIndex)

"""
Machine-local sample files for the `.zar`/`.zmf` reader tests in
`test/zemax.jl`. Both are real vendor/Zemax-derived files that live
outside the repo (possible IP concerns), not part of the repo -- same
situation as [`HAS_GLASS_CATALOG`](@ref) above. Tests that need them
should check the corresponding `HAS_...` constant and skip (not fail)
when it's false, e.g.:

    if HAS_ZAR_SAMPLE
        @testset "..." begin ... end
    else
        @info "Skipping ...: $ZAR_SAMPLE_PATH not found"
    end
"""
const ZAR_SAMPLE_PATH = "/Users/matt/Desktop/Zemax lenses/LF1988-Zemax.zar"
const HAS_ZAR_SAMPLE = isfile(ZAR_SAMPLE_PATH)

const ZMF_SAMPLE_PATH = "/Users/matt/Development/Software/Zemax/Stockcat/DiverseOptics.ZMF"
const HAS_ZMF_SAMPLE = isfile(ZMF_SAMPLE_PATH)

"""
Directory that `test/plotting.jl` saves its generated figures into (via
[`saveTestFigure`](@ref)), so they can be inspected after a test run.
Gitignored -- not part of the repo -- so it's created here rather than
relying on it being checked out.
"""
const TEST_PLOTS_DIR = joinpath(@__DIR__, "plots")
mkpath(TEST_PLOTS_DIR)

"""
Some plotting functions under test (e.g. `plotPerimeterRays`,
`rayHeatmap`) return a `Makie.FigureAxisPlot` rather than a bare
`Figure` when they create their own figure -- normalize to the
underlying `Figure` so callers like [`saveTestFigure`](@ref) can add a
title `Label` and save it uniformly.
"""
figureOf(fig) = fig
figureOf(fig::Makie.FigureAxisPlot) = fig.figure

"""
    saveTestFigure(fig, name)

Adds a title `Label` reading `name` to `fig` (or, if `fig` is a
`Makie.FigureAxisPlot`, to its underlying figure -- see
[`figureOf`](@ref)), then saves it into [`TEST_PLOTS_DIR`](@ref) via
`saveFigure`, so `test/plotting.jl` leaves a titled, identifiable image
behind for every figure it generates.
"""
function saveTestFigure(fig, name)
    f = figureOf(fig)
    Label(f[0, :], name, tellwidth = false)
    saveFigure(name, f; directory = TEST_PLOTS_DIR)
end