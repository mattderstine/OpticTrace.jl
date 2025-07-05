@testset "Allocation tests" begin
    @test 0 == @allocations sag(0., 1., OpticTrace.SurfProfileConic(0.01, 0.)) 
end

