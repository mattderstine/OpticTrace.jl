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
