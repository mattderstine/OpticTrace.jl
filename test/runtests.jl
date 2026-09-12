using OpticTrace
using BenchmarkTools
using Test
using GeometryBasics
using ForwardDiff
using LinearAlgebra
using StaticArrays
using GLMakie
using StatsBase


include("helper.jl")


@testset "OpticTrace.jl" begin
    # Write your tests here.

    #include("testing.jl")
    include("allocations.jl")
    include("foundations.jl")
    include("optics.jl")
    include("surface_builders.jl")
    include("mesh_primitives.jl")
    include("trace_geometry.jl")
    include("surface_manipulation.jl")
    include("characterization.jl")
    include("refractive_index.jl")
    include("agf.jl")
    include("lens_catalogs.jl")
    include("zemax.jl")
    include("zemax_browser.jl")
    include("filepicker.jl")
    include("printing.jl")
    include("plotting.jl")



end
