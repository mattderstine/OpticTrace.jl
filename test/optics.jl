


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

    offsetOA = Vec3(0.3, -0.2, 0.05)
    spOA = OpticTrace.SurfProfileOAConic(curve1, 0.0, offsetOA)
    baseConicOA = OpticTrace.SurfProfileConic(curve1, 0.0)

    spsCyl1 = OpticTrace.SurfProfileCyl(curve1, 0.0, Float64[])
    spsCyl2 = OpticTrace.SurfProfileCyl(curve2, 1.0, Float64[])

    curvYToroid = 0.3
    curvXToroid = 0.2
    spsToroid = OpticTrace.SurfProfileToroid(curvYToroid, curvXToroid)

    # `simplesystem` is unused elsewhere in this file, but building it
    # calls lens_TLAC254_060, which looks up real glass files from
    # OpticTrace.dirBaseRefractiveIndex -- a machine-local directory,
    # not part of the repo (see TODO.md). Guarded so this file still
    # runs (skipping only this dead assignment) on a checkout without
    # that directory, e.g. CI.
    if HAS_GLASS_CATALOG
        simplesystem = [[referencePlane("object", ORIGIN, ZAXIS, 1.0, 10.0, "test coating")]
            lens_TLAC254_060(Point(0.0, 0.0, 10.0), ZAXIS, 0.45; order="forward", lensname="TL_AC254-060");
            [referencePlane("image", Point(0.0, 0.0, 70.0), ZAXIS, 1.0, 10.0, "test coating")]
        ]
    end

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

        # SurfProfileCyl: x does not appear in the formula (only y does)
        sag_value = sag(1.0, 1.0, spsCyl1) # ϵ=0 branch: z = curv*y^2/2
        @test sag_value ≈ curve1 * 1.0^2 * 0.5
        @test sag(5.0, 1.0, spsCyl1) == sag_value # x is ignored

        sag_value = sag(1.0, 1.0, spsCyl2) # ϵ=1 branch (circular cross-section)
        @test sag_value ≈ sag(0.0, 1.0, sps) # reduces to SurfProfileSphere's formula along y
        @test sag(-3.0, 1.0, spsCyl2) == sag_value # x is ignored

        # SurfProfileToroid: author-flagged "likely incorrect" (see TODO.md).
        # A toroid's x=0 cross-section should reduce to a plain circular sag
        # along y with curvature curvY, matching SurfProfileSphere -- the
        # current formula is dimensionally inconsistent with that (it's
        # missing the division by curvY that every other sag formula in
        # this codebase has), so this is expected to fail.
        expectedToroidSag = sag(0.0, 0.5, OpticTrace.SurfProfileSphere(curvYToroid))
        actualToroidSag = sag(0.0, 0.5, spsToroid)
        @test_broken actualToroidSag ≈ expectedToroidSag
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

        #off-axis conic test (sag & surfNormal correct; deltaToSurf has a
        #known bug -- author's own comment flags "logic is flawed in this
        #one", see TODO.md)
        x0, y0 = 0.1, 0.15
        z0 = sag(x0, y0, spOA)
        @test z0 ≈ sag(x0 - offsetOA[1], y0 - offsetOA[2], baseConicOA) + offsetOA[3]

        normalOA = OpticTrace.surfNormal(Point3(x0, y0, z0), spOA)
        expectedNormalOA = OpticTrace.surfNormal(Point3(x0, y0, z0) .- offsetOA, baseConicOA)
        @test normalOA ≈ expectedNormalOA
        @test normal_from_sag(x0, y0, spOA) ≈ normalize(normalOA)

        oaRay = Ray(Point3(x0, y0, -1.0), Vec3(0.0, 0.0, 1.0))
        oaDelta = OpticTrace.deltaToSurf(oaRay, spOA)
        oaPoint = rprop(oaRay, oaDelta)
        @test_broken oaPoint ≈ Point3(x0, y0, z0)

        #cylinder test
        xC, yC = 0.3, 0.4
        zC = sag(xC, yC, spsCyl2)
        pointC = Point3(xC, yC, zC)
        normalC = OpticTrace.surfNormal(pointC, spsCyl2)
        @test normal_from_sag(xC, yC, spsCyl2) ≈ normalize(normalC)

        # deltaToSurf, ϵ=1 branch: reduces to the same circular
        # cross-section formula as the sphere test above, since
        # curve2 == curve and ray's x is 0
        sag_valueCyl = curve2 * r2 / (1.0 + sqrt(1 - curve2^2 * r2))
        deltaCyl = OpticTrace.deltaToSurf(ray, spsCyl2)
        @test deltaCyl ≈ 1.0 + sag_valueCyl
        pointCyl = rprop(ray, deltaCyl)
        @test pointCyl ≈ Point3(0.0, offset, sag_valueCyl)
        normalCyl = OpticTrace.surfNormal(pointCyl, spsCyl2)
        @test normal_from_sag(pointCyl, spsCyl2) ≈ normalize(normalCyl)
    end



    @testset "modFunc tests" begin


        normal = Vec3(0.0, 0.0, 1.0)
        ray = Ray(Point(0.0, 0.0, 0.0), Vec(0.0, 0.4, sqrt(1.0 - 0.4^2)))
        dT = OpticTrace.DielectricT(1.0, 1.5)
        status, dir, nIn = OpticTrace.modFunc(ray, normal, dT)
        @test status == true
        @test nIn == dT.refIndexIn
        @test dT.refIndexOut * dir[2] ≈ dT.refIndexIn * ray.dir[2]

        dR = OpticTrace.MirrorR(1.0, 1.0)
        status, dir, nIn = OpticTrace.modFunc(ray, normal, dR)
        @test status == true
        @test nIn ≈ dR.refIndexIn
        @test dir ≈ Vec3(0.0, 0.4, -sqrt(1.0 - 0.4^2))

        @testset "CDiffuser" begin
            a = normalize(Vec3(0.1, 0.2, 0.97))
            θ = 0.15
            cd = OpticTrace.CDiffuser(tan(θ), 1.0, 1.5)
            diffRay = Ray(Point3(0.0, 0.0, 0.0), a)

            dirs = Vector{Vec3{Float64}}()
            for _ in 1:20
                status, dir, nIn = OpticTrace.modFunc(diffRay, normal, cd)
                @test status == true
                @test nIn == cd.refIndexIn # a ⋅ normal > 0 -> hitting "forward"
                @test norm(dir) ≈ 1.0
                @test dot(dir, a) >= cos(θ) - 1e-9 # stays within the diffuser cone
                push!(dirs, dir)
            end
            @test length(unique(dirs)) > 1 # scattering is randomized
        end

        @testset "NoBendIndex" begin
            nbRay = Ray(Point3(0.0, 0.0, 0.0), Vec3(0.1, 0.2, sqrt(1.0 - 0.1^2 - 0.2^2)))
            nb = OpticTrace.NoBendIndex(1.33)
            status, dir, nIn = OpticTrace.modFunc(nbRay, normal, nb)
            @test status == true
            @test dir == nbRay.dir
            @test nIn == nb.refIndexIn
        end

    end
    #=
        a test set for surfNormal methods
        The tests checks that the normal vector of different surfaces and different intersection points are correct.

        This is done by comparing the normal vector obtained from the surfNormal method
        with the expected normal vector found using the gradient of the sag methods found using and autodiff package.
        The first step for this is to find the intersection point of the ray with the surface, then find the sag value at that point,
             and then find the gradient of the sag value at that point. The normal vector is then obtained by normalizing the gradient vector.
    =#

end
