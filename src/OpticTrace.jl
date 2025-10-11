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
import YAML


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
end
