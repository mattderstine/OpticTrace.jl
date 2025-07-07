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

        sps = OpticTrace.SurfProfileAsphere(0.5, 0.0, [0.0, 0.01, 0.1])
        sag_value = sag(1.0, 1.0, sps)
        @test sag_value ≈ 1.1056854249492383 #this assumes code is correct as of 7 july 2025    

        #Tests not implemented for SurfProfileCyl & SurfProfileToroid
    end
    @testset "deltaToSurf tests" begin

        #=
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileSphere(1.0)) == 0.0
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileConic(1.0, 0.0)) == 0.0
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileAsphere(1.0, 0.0, [0.0, 0.01, 0.1])) == 0.0
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileCyl(1.0, 0.0, [0.0, 0.01, 0.1])) == 0.0
        =#
    end
end
    #@testset "modFunc tests" begin
        normal =  Vec3( 0.0, 0.0, 1.0)
        ray = Ray(Point(0.0, 0.0, 0.0), Vec(0.0, 0.4, sqrt(1.0 - 0.4^2)))
        dT= OpticTrace.DielectricT(1.0, 1.5)
        status,dir, nIn = OpticTrace.modFunc(ray, normal, dT ) 
        @test status == true
        @test nIn == dT.refIndexIn
        @test dT.refIndexOut * dir[2] ≈ dT.refIndexIn * ray.dir[2]

    #end
#end
