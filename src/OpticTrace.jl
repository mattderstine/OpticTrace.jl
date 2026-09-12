module OpticTrace


using StaticArrays
using FileIO # MeshIO should also be installed
using LinearAlgebra
using CoordinateTransformations
using IterTools
using Roots
using DataInterpolations
using StatsBase
using GeometryBasics
using GLMakie
# `import`, not `using` -- same collision-avoidance reason as the
# `import Bonito` note below (WGLMakie itself `using Bonito`s
# internally, and re-exports most of Makie, which would double up on
# GLMakie's already-`using`'d Button/Slider/Checkbox/Dropdown).
#
# Backend-activation policy: OpticTrace does not call `activate!()` on
# either backend itself, and has no `__init__` reclaiming one as
# default. `GLMakie.__init__` and `WGLMakie.__init__` each
# unconditionally call their own `activate!()` when loaded, so simply
# loading this package leaves whichever one's `__init__` happens to run
# last (an unspecified implementation detail of module load order) as
# the ambient default -- callers are responsible for calling
# `GLMakie.activate!()` or `WGLMakie.activate!()` themselves before
# using any OpticTrace plotting function that relies on the ambient
# backend (`saveFigure`, `multipleFigures`, any `display(fig)`-based
# usage), per normal Makie multi-backend convention. The one exception
# is `src/zemax_browser.jl`'s WGLMakie live geometry preview, which
# doesn't depend on the ambient backend at all -- see that file's
# `renderGeometryPreview` docstring.
import WGLMakie
using Printf
using Optim
using ForwardDiff
import YAML
# `import`, not `using`: Bonito and GLMakie (already `using`'d above) export
# several identical names bound to unrelated types (Button, Slider, Checkbox,
# Dropdown); `using Bonito` here would make those names ambiguous/undefined
# module-wide, breaking src/plotting.jl's unqualified use of GLMakie's Button
# in multipleFigures. src/zemax_browser.jl and src/UItools/filepicker.jl are
# the files that actually use this binding -- both qualify every reference
# (Bonito.App, Bonito.DOM, ...) rather than `using` it themselves, relying on
# this module-level import rather than their own local one.
import Bonito


include("constants.jl")
include("lens_definitions.jl")
include("lens_refractive_index.jl")
include("agf.jl")
include("mesh_primitives.jl")
include("plotting.jl")
include("printing.jl")
include("tracing.jl")
include("surfaces.jl")
include("aperture.jl")
include("extended_geo.jl")
include("lens_edmund.jl")
include("lens_thorlabs.jl")
include("characterization.jl")
include("zemax.jl")
include("surface_manipulation.jl")
include("UItools/filepicker.jl")
include("zemax_browser.jl")

end
