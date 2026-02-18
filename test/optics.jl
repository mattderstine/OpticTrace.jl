


@testset "optics.jl" begin
    # Write your tests here.

    # test objects

    asphere_coeff1 = 0.1
    asphere_coeff2 = 0.01
    factor1 = 4.0
    factor2 = 8.0
    curve1 = 1.0
    curve2 = sqrt(0.5)
    curve3 = 0.0
    curve = curve2
    spc1 = OpticTrace.SurfProfileConic(curve1, 0.0)
    spc2 = OpticTrace.SurfProfileConic(curve2, 1.0)
    spc3 = OpticTrace.SurfProfileConic(curve3, 1.0)
    sps = OpticTrace.SurfProfileSphere(curve)
    spsAsphere = OpticTrace.SurfProfileAsphere(0.0, 0.0, [0.00,  asphere_coeff1,0.0, asphere_coeff2])
    spsEAsphere = OpticTrace.SurfProfileEvenAsphere(0.0, 0.0, [asphere_coeff1, asphere_coeff2])

    simplesystem = [[referencePlane("object", ORIGIN, ZAXIS, 1.0, 10.0, "test coating")]
        lens_TLAC254_060(Point(0.0, 0.0, 10.0), ZAXIS, 0.45; order="forward", lensname="TL_AC254-060");
        [referencePlane("image", Point(0.0, 0.0, 70.0), ZAXIS, 1.0, 10.0, "test coating")]
    ]

    #sag tests
    @testset "sag tests" begin

        sag_value = sag(1.0, 1.0, spc1)
        @test sag_value ≈ 1.0

        sag_value = sag(0.9999999999999999, 1.0, spc2)
        @test sag_value ≈ sqrt(2)

        sag_value = sag(1.0, 1.0, spc3)
        @test sag_value == 0.0

        sag_value = sag(0.9999999999999999, 1.0, sps)
        @test sag_value ≈ sqrt(2)

        sag_value = sag(100.0, 1.0, sps)
        @test isnan(sag_value)

        sag_value = sag(100.0, 1.0, spc2)
        @test isnan(sag_value)



        #r = sqrt(2) r2=2, r4 = 4, r6 = 8

        sag_value = sag(1.0, 1.0, spsAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2

        sag_value = sag(1.0, 1.0, spsEAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2
        sag_value = sag(1.0, 1.0, spsAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2

        sag_value = sag(1.0, 1.0, spsEAsphere)
        @test sag_value ≈ factor1 * asphere_coeff1 + factor2 * asphere_coeff2

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

    @testset "surfNormal & deltaToSurf tests" begin



        offset = 0.2
        r2 = offset^2
        ray = Ray(Point(0.0, offset, -1.0), Vec(0.0, 0.0, 1.0))

        #sphere test
        sag_value = curve * r2 / (1.0 + sqrt(1 - curve^2 * r2))
        delta = OpticTrace.deltaToSurf(ray, sps)
        @test delta ≈ 1.0+sag_value
        point = rprop(ray, delta)
        @test point ≈ Point3(0.0, offset, sag_value)
        normal = OpticTrace.surfNormal(point, sps)
        @test normal_from_sag(point, sps) ≈ normal

        #paraboloid test
        sag_value = curve1 * offset^2 * 0.5 
        delta = OpticTrace.deltaToSurf(ray, spc1)
        @test delta ≈ 1.0+sag_value
        point = rprop(ray, delta)
        @test point ≈ Point3(0.0, offset, sag_value)
        normal = OpticTrace.surfNormal(point, spc1)
        @test normal_from_sag(point, spc1) ≈ normal

        #asphere test
        sag_value = offset^4 * asphere_coeff1 + offset^6 * asphere_coeff2
        delta = OpticTrace.deltaToSurf(ray, spsAsphere)
        @test delta ≈ 1.0+sag_value
        point = rprop(ray, delta)
        @test point ≈ Point3(0.0, offset, sag_value)
        normal = OpticTrace.surfNormal(point, spsAsphere)
        @test normal_from_sag(point, spsAsphere) ≈ normalize(normal)

        #even asphere test
        sag_value = offset^4 * asphere_coeff1 + offset^6 * asphere_coeff2
        delta = OpticTrace.deltaToSurf(ray, spsEAsphere)
        @test delta ≈ 1.0+sag_value
        point = rprop(ray, delta)
        @test point ≈ Point3(0.0, offset, sag_value)
        normal = OpticTrace.surfNormal(point, spsEAsphere)
        @test normal_from_sag(point, spsEAsphere) ≈ normalize(normal)


    end



    @testset "modFunc tests" begin
        normal = Vec3(0.0, 0.0, 1.0)
        ray = Ray(Point(0.0, 0.0, 0.0), Vec(0.0, 0.4, sqrt(1.0 - 0.4^2)))
        dT = OpticTrace.DielectricT(1.0, 1.5)
        status, dir, nIn = OpticTrace.modFunc(ray, normal, dT)
        @test status == true
        @test nIn == dT.refIndexIn
        @test dT.refIndexOut * dir[2] ≈ dT.refIndexIn * ray.dir[2]

    end
    #=
        a test set for surfNormal methods
        The tests checks that the normal vector of different surfaces and different intersection points are correct.

        This is done by comparing the normal vector obtained from the surfNormal method
        with the expected normal vector found using the gradient of the sag methods found using and autodiff package.
        The first step for this is to find the intersection point of the ray with the surface, then find the sag value at that point,
             and then find the gradient of the sag value at that point. The normal vector is then obtained by normalizing the gradient vector.
    =#




    
    @testset "Lens Tests" begin


    end

end
 