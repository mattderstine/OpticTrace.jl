using OpticTrace
using BenchmarkTools
using Test
using GeometryBasics




@testset "OpticTrace.jl" begin
    # Write your tests here.

    #include("testing.jl")
    include("allocations.jl")
    include("optics.jl")



end
