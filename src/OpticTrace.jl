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
include("lens definitions.jl")
include("lens refractive index.jl")
include("mesh primitives.jl")
include("plotting.jl")
include("tracing.jl")
include("extended geo.jl")
end
