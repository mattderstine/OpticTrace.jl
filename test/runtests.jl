using OpticTrace
using BenchmarkTools
using Test
using GeometryBasics
using ForwardDiff
using LinearAlgebra


include("helper.jl")


@testset "OpticTrace.jl" begin
    # Write your tests here.

    #include("testing.jl")
    include("allocations.jl")
    include("optics.jl")



end
