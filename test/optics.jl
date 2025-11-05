
# test objects
spc1 = OpticTrace.SurfProfileConic(1.0, 0.0)
spc2 = OpticTrace.SurfProfileConic(sqrt(0.5), 1.0 )
spc3 = OpticTrace.SurfProfileConic(0.0, 1.0)
sps = OpticTrace.SurfProfileSphere(sqrt(0.5))
spsAsphere = OpticTrace.SurfProfileAsphere(0.0, 0.0, [0.00, 0.1, 0.0, 0.01])
spsEAsphere = OpticTrace.SurfProfileEvenAsphere(0.0, 0.0, [0.1, 0.01])

@testset "optics.jl" begin
    # Write your tests here.
    #sag tests
    @testset "sag tests" begin

        sag_value = sag(1.0, 1.0, spc1)
        @test sag_value ≈ 1.0

        sag_value = sag(0.9999999999999999, 1.0, spc2)
        @test sag_value ≈ sqrt(2)

        sag_value = sag(1.0, 1.0, spc3)
        @test sag_value == 0.0  

        sag_value = sag(0.9999999999999999, 1.0, sps)
        @test sag_value ≈ 1.0  


        asphere_coeff1 = 0.1
        asphere_coeff2 = 0.01
        factor1 = 4.0
        factor2 = 8.0

        #r = sqrt(2) r2=2, r4 = 4, r6 = 8

        sag_value = sag(1.0, 1.0, spsAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2

        sag_value = sag(1.0, 1.0, spsEAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2
        sag_value = sag(1.0, 1.0, spsAsphere)
        @test sag_value ≈ asphere_factor1 * asphere_coeff1 + asphere_factor2 * asphere_coeff2

        sag_value = sag(1.0, 1.0, spsEAsphere)
        @test sag_value ≈ asphere_factor1 * asphere_coeff1 + asphere_factor2 * asphere_coeff2

        #=
        spsC = OpticTrace.SurfProfileCyl(1.0, 0.0, [0.00, 0.1, 0.0, 0.01])
        sag_value = sag(1.0, 1.0, spsC)
        @test sag_value ≈ 1.0

        spsT = OpticTrace.SurfProfileToroid(1.0, sqrt(0.5), [0.00, 0.1, 0.0, 0.01])
        sag_value = sag(1.0, 1.0, spsT)
        @test sag_value ≈ 1.4142135623730951   

        #Tests not implemented for SurfProfileCyl & SurfProfileToroid
        =#
    end
    @testset "deltaToSurf tests" begin

        #=
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileSphere(1.0)) == 0.0
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileConic(1.0, 0.0)) == 0.0
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileAsphere(1.0, 0.0, [0.0, 0.01, 0.1])) == 0.0
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileCyl(1.0, 0.0, [0.0, 0.01, 0.1])) == 0.0
        =#
    end

    @testset "modFunc tests" begin
        normal =  Vec3( 0.0, 0.0, 1.0)
        ray = Ray(Point(0.0, 0.0, 0.0), Vec(0.0, 0.4, sqrt(1.0 - 0.4^2)))
        dT= OpticTrace.DielectricT(1.0, 1.5)
        status,dir, nIn = OpticTrace.modFunc(ray, normal, dT ) 
        @test status == true
        @test nIn == dT.refIndexIn
        @test dT.refIndexOut * dir[2] ≈ dT.refIndexIn * ray.dir[2]

    end

    @testset "surfNormal tests" begin
        normal = OpticTrace.surfNormal(ORIGIN, spsAsphere)
        @test normal ≈ [0.0, 0.0, 1.0]
    end

    @testset "Lens Tests" begin
        
        g = [referencePlane(ORIGIN, ZAXIS, 10.0),
            lens_TLAC254_060( Point(0.0, 0.0, 10.0), ZAXIS,  0.45; order = "forward", lensname = "TL_AC254-060"),
            referencePlane(Point(0.0, 0.0, ),   ZAXIS, 10.0)
        ]
    end

end
