
# test objects
spc1 = OpticTrace.SurfProfileConic(1.0, 0.0)
spc2 = OpticTrace.SurfProfileConic(sqrt(0.5), 1.0 )
spc3 = OpticTrace.SurfProfileConic(0.0, 1.0)
sps = OpticTrace.SurfProfileSphere(sqrt(0.5))
spsAsphere = OpticTrace.SurfProfileAsphere(0.0, 0.0, [0.00, 0.1, 0.0, 0.01])
spsEAsphere = OpticTrace.SurfProfileEvenAsphere(0.0, 0.0, [0.1, 0.01])

simplesystem = [[referencePlane("object",ORIGIN, ZAXIS, 1.0, 10.0, "test coating")]
            lens_TLAC254_060( Point(0.0, 0.0, 10.0), ZAXIS,  0.45; order = "forward", lensname = "TL_AC254-060")
            [referencePlane("image",Point(0.0, 0.0, 70.0),   ZAXIS, 1.0, 10.0, "test coating")]
            ]
        

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
    @testset "deltaToSurf tests" begin

        @test OpticTrace.deltaToSurf(Ray(Point(0.0, 0.0, -1.0), Vec(0.0, 0.0, 1.0)), spsAsphere) == 1.0

        #=
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileSphere(1.0)) == 0.0
        @test OpticTrace.deltaToSurf(0.0, OpticTrace.SurfProfileConic(1.0, 0.0)) == 0.0
   
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
    #=
        a test set for surfNormal methods
        The tests checks that the normal vector of different surfaces and different intersection points are correct.
        
        This is done by comparing the normal vector obtained from the surfNormal method
        with the expected normal vector found using the gradient of the sag methods found using and autodiff package.
        The first step for this is to find the intersection point of the ray with the surface, then find the sag value at that point,
             and then find the gradient of the sag value at that point. The normal vector is then obtained by normalizing the gradient vector.
    =#

    """
    define function to compute the normal vector using the gradient of the sag function
     This function takes a ray and a surface profile as input, and returns the normal vector at the intersection point of the ray with the surface.
    The function first finds the intersection point of the ray with the surface using the deltaToSurf method, then finds the sag value at the intersection point, and finally finds the gradient of the sag value at the intersection point using an autodiff package. The normal vector is then obtained by normalizing the gradient vector.
    """
    function normal_from_sag(ray::Ray, surfProfile)
        # find the intersection point of the ray with the surface using the deltaToSurf method
        delta = OpticTrace.deltaToSurf(ray, surfProfile)
        intersection_point = ray.origin + delta * ray.dir
        # find the gradient of the sag value at the intersection point using an autodiff package
        grad_sag = ForwardDiff.gradient((x) -> sag(x[1], x[2], surfProfile), intersection_point[1:2])
        # the normal vector is then obtained by normalizing the gradient vector
        normal_vector = normalize(Vec3(grad_sag[1], grad_sag[2], -1.0))
        return normal_vector
    end

  
    @testset "surfNormal tests" begin
        normal = OpticTrace.surfNormal(ORIGIN, spsAsphere)
        @test normal ≈ [0.0, 0.0, 1.0]



        

    end

    @testset "Lens Tests" begin
        

    end

end
