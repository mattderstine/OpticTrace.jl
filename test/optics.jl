@testset "optics.jl" begin
    # Write your tests here.
    #sag tests
    @testset "sag tests" begin
        spc = OpticTrace.SurfProfileConic(1.0, 0.0)
        sag_value = sag(1.0, 1.0, spc)
        @test sag_value ≈ 1.0

        spc = OpticTrace.SurfProfileConic(sqrt(0.5), 1.0 )
        sag_value = sag(0.9999999999999999, 1.0, spc)
        @test sag_value ≈ sqrt(2)

        spc = OpticTrace.SurfProfileConic(0.0, 1.0)
        sag_value = sag(1.0, 1.0, spc)
        @test sag_value == 0.0  

        sps = OpticTrace.SurfProfileSphere(sqrt(0.5))
        sag_value = sag(0.9999999999999999, 1.0, sps)
        @test sag_value ≈ 1.0  
    end

end

