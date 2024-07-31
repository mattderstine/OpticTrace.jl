@testset "Allocation tests" begin
    @test 0 == @allocations pp = sag(0., 1., OpticTrace.SurfProfileConic(0.01, 0.)) 
    @test pp = 0.005
    surfprofile = OpticTrace.SurfProfileAsphere(0.01, 0, SVector(0., 0.1, -0.1))
    @test 0 = @allocations pp = sag(0. , 1., surfprofile)
end

